#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <algorithm>
#include <fstream>

// Helper function to find the maximum value in a vector
// Inputs: a vector of integers
// Outputs: the maximum value in the vector, or the minimum possible integer if the vector is empty
int find_max(const std::vector<int>& values) {
    if (values.empty()) {
        return std::numeric_limits<int>::min();
    }

    return *std::max_element(values.begin(), values.end());
}

// Helper function to reverse a string in place
// Inputs: a pointer to a character array (string) and its length
// Outputs: the string is reversed in place
void reverse_string(char* str, unsigned int length) {
    for (unsigned int i = 0, j = length - 1; i < length / 2; ++i, --j) {
        std::swap(str[i], str[j]);
    }
}

// Enum for traceback direction in the alignment matrix
enum class Direction {
    MATCH,       // Represents a match or mismatch in the alignment
    INSERTION,   // Represents an insertion in the alignment
    DELETION,    // Represents a deletion in the alignment
    NONE         // Represents no alignment direction (used for initialization)
};

// Enum for alignment mode (local or global)
enum class AlignmentMode {
    LOCAL,       // Local alignment mode
    GLOBAL       // Global alignment mode
};

// Struct to represent a block of the alignment matrix to be processed by a worker thread
struct MatrixBlock {
    unsigned int startRow; // Starting row of the block
    unsigned int startCol; // Starting column of the block
};

// Class for parallel sequence alignment
class ParallelSequenceAligner {
public:
    char* seq1; // Sequence 1
    char* seq2; // Sequence 2
    unsigned int len1, len2; // Lengths of sequences
    int gap_penalty, match_score, mismatch_score; // Scoring parameters
    unsigned int num_threads, block_height, block_width; // Thread and block parameters

    std::vector<std::thread> workers; // Worker threads
    std::mutex* mtx; // Mutex array for thread synchronization
    std::condition_variable* cv; // Condition variable array for thread synchronization

    std::atomic_int completed_threads; // Counter for completed threads
    std::atomic_int current_phase; // Current phase of computation
    int total_phases; // Total number of phases

    std::vector<std::vector<int>> score_matrix; // Alignment score matrix
    std::vector<std::vector<Direction>> traceback_matrix; // Traceback matrix
    char* aligned_seq1; // Aligned sequence 1
    char* aligned_seq2; // Aligned sequence 2
    unsigned int aligned_len1, aligned_len2; // Lengths of aligned sequences
    AlignmentMode mode; // Alignment mode (local or global)

    // Constructor to initialize the aligner with sequences and parameters
    // Inputs: sequences (seq1, seq2), their lengths (len1, len2), number of threads, scoring parameters (gap_penalty, match_score, mismatch_score), and alignment mode
    ParallelSequenceAligner(
        char* seq1, char* seq2,
        unsigned int len1, unsigned int len2,
        unsigned int num_threads,
        int gap_penalty, int match_score, int mismatch_score,
        AlignmentMode mode = AlignmentMode::LOCAL)
        : seq1(seq1), seq2(seq2), len1(len1), len2(len2), num_threads(num_threads),
          gap_penalty(gap_penalty), match_score(match_score), mismatch_score(mismatch_score), mode(mode) {
        mtx = new std::mutex[num_threads];
        cv = new std::condition_variable[num_threads];

        score_matrix.resize(len1 + 1, std::vector<int>(len2 + 1, 0));
        traceback_matrix.resize(len1 + 1, std::vector<Direction>(len2 + 1, Direction::NONE));

        completed_threads.store(0);
        current_phase.store(1);

        block_height = (len1 + num_threads - 1) / num_threads;
        block_width = (len2 + num_threads - 1) / num_threads;
        total_phases = (len1 + block_height - 1) / block_height + (len2 + block_width - 1) / block_width - 1;

        aligned_seq1 = new char[len1 + len2 + 2];
        aligned_seq2 = new char[len1 + len2 + 2];
        aligned_len1 = 0;
        aligned_len2 = 0;
    }

    // Destructor to clean up allocated resources
    ~ParallelSequenceAligner() {
        delete[] mtx;
        delete[] cv;
        delete[] aligned_seq1;
        delete[] aligned_seq2;
    }

    // Function to start the alignment process using worker threads
    // This function initializes and starts the worker threads, then waits for them to complete, and finally performs the traceback to generate the aligned sequences
    void align_sequences() {
        for (unsigned int i = 0; i < num_threads; ++i) {
            workers.emplace_back(&ParallelSequenceAligner::worker_task, this, i + 1);
        }

        for (auto& worker : workers) {
            if (worker.joinable())
                worker.join();
        }

        trace_back();
    }

    // Function executed by each worker thread
    // This function is responsible for computing the alignment matrix in parallel, working in phases
    // Inputs: worker_id - the ID of the worker thread
    void worker_task(unsigned int worker_id) {
        int local_phase = 1;
        while (local_phase <= total_phases) {
            std::vector<MatrixBlock> blocks_to_compute;
            assign_blocks(worker_id, blocks_to_compute, local_phase);

            for (const auto& block : blocks_to_compute) {
                for (unsigned int i = block.startRow; i < std::min(block.startRow + block_height, len1 + 1); ++i) {
                    for (unsigned int j = block.startCol; j < std::min(block.startCol + block_width, len2 + 1); ++j) {
                        compute_cell(i, j);
                    }
                }
            }

            local_phase++;
            unsigned int fetched = completed_threads.fetch_add(1);
            if (fetched == num_threads - 1) {
                int expected = fetched + 1;
                while (!completed_threads.compare_exchange_weak(expected, 0)) {
                    expected = fetched + 1;
                }
                current_phase.fetch_add(1);

                for (unsigned int i = 0; i < num_threads; ++i) {
                    cv[i].notify_one();
                }
                continue;
            }
            std::unique_lock<std::mutex> lock(mtx[worker_id - 1]);
            while (current_phase.load() != local_phase)
                cv[worker_id - 1].wait(lock);
        }
    }

    // Function to assign blocks of the alignment matrix to a worker thread for a given phase
    // Inputs: worker_id - the ID of the worker thread, blocks - reference to a vector to store assigned blocks, phase - the current phase of computation
    // Outputs: the assigned blocks are added to the blocks vector
    void assign_blocks(unsigned int worker_id, std::vector<MatrixBlock>& blocks, unsigned int phase) {
        if (worker_id > num_threads) {
            std::cerr << "Invalid: worker_id is greater than the number of threads." << std::endl;
            return;
        }

        int i = 1 + (worker_id - 1) * block_height;
        int j = 1 + (phase - worker_id) * block_width;

        while (i < static_cast<int>(len1) + 1 && j >= 1) {
            blocks.push_back({static_cast<unsigned int>(i), static_cast<unsigned int>(j)});
            i += num_threads * block_height;
            j -= num_threads * block_width;
        }

        if (i < static_cast<int>(len1) + 1 && j >= 1) {
            blocks.push_back({static_cast<unsigned int>(i), static_cast<unsigned int>(j)});
        }
    }

    // Function to compute the score for a single cell in the alignment matrix
    // Inputs: i - row index, j - column index
    // This function calculates the match/mismatch, insertion, and deletion scores, then selects the best score and updates the score matrix and traceback matrix accordingly
    void compute_cell(unsigned int i, unsigned int j) {
        int match_mismatch = score_matrix[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match_score : mismatch_score);
        int insertion = score_matrix[i - 1][j] - gap_penalty;
        int deletion = score_matrix[i][j - 1] - gap_penalty;

        if (mode == AlignmentMode::LOCAL) {
            score_matrix[i][j] = std::max({0, match_mismatch, insertion, deletion});
        } else {
            score_matrix[i][j] = std::max({match_mismatch, insertion, deletion});
        }

        if (score_matrix[i][j] == match_mismatch) {
            traceback_matrix[i][j] = Direction::MATCH;
        } else if (score_matrix[i][j] == insertion) {
            traceback_matrix[i][j] = Direction::INSERTION;
        } else if (score_matrix[i][j] == deletion) {
            traceback_matrix[i][j] = Direction::DELETION;
        } else {
            traceback_matrix[i][j] = Direction::NONE;
        }
    }

    // Function to perform traceback and generate the aligned sequences
    // This function traces back through the traceback matrix to build the aligned sequences based on the best scores calculated during matrix computation
    void trace_back() {
        unsigned int i = len1, j = len2;
        if (mode == AlignmentMode::LOCAL) {
            int max_score = 0;
            for (unsigned int x = 0; x <= len1; ++x) {
                for (unsigned int y = 0; y <= len2; ++y) {
                    if (score_matrix[x][y] > max_score) {
                        max_score = score_matrix[x][y];
                        i = x;
                        j = y;
                    }
                }
            }
        }

        while (i > 0 || j > 0) {
            if (traceback_matrix[i][j] == Direction::MATCH) {
                aligned_seq1[aligned_len1++] = seq1[i - 1];
                aligned_seq2[aligned_len2++] = seq2[j - 1];
                i--;
                j--;
            } else if (traceback_matrix[i][j] == Direction::INSERTION) {
                aligned_seq1[aligned_len1++] = seq1[i - 1];
                aligned_seq2[aligned_len2++] = '-';
                i--;
            } else if (traceback_matrix[i][j] == Direction::DELETION) {
                aligned_seq1[aligned_len1++] = '-';
                aligned_seq2[aligned_len2++] = seq2[j - 1];
                j--;
            } else {
                break;
            }
        }

        reverse_string(aligned_seq1, aligned_len1);
        reverse_string(aligned_seq2, aligned_len2);
    }
};

// Function to run a single test and log the results
// Inputs: seq1 and seq2 are the sequences to align, num_threads is the number of threads to use, scoring parameters (gap_penalty, match_score, mismatch_score), alignment mode, and a reference to the log file
// Outputs: logs the original and aligned sequences to the log file
void run_test(const std::string& seq1, const std::string& seq2, unsigned int num_threads, int gap_penalty, int match_score, int mismatch_score, AlignmentMode mode, std::ofstream& log_file) {
    unsigned int len1 = seq1.size();
    unsigned int len2 = seq2.size();
    ParallelSequenceAligner aligner(const_cast<char*>(seq1.c_str()), const_cast<char*>(seq2.c_str()), len1, len2, num_threads, gap_penalty, match_score, mismatch_score, mode);
    aligner.align_sequences();

    log_file << "Test with sequences:\n";
    log_file << "Sequence 1: " << seq1 << "\n";
    log_file << "Sequence 2: " << seq2 << "\n";
    log_file << "Aligned Sequence 1: ";
    for (unsigned int i = 0; i < aligner.aligned_len1; ++i) {
        log_file << aligner.aligned_seq1[i];
    }
    log_file << "\nAligned Sequence 2: ";
    for (unsigned int i = 0; i < aligner.aligned_len2; ++i) {
        log_file << aligner.aligned_seq2[i];
    }
    log_file << "\n\n";
}

// Function to generate a random DNA sequence of a given length
// Inputs: length - the length of the sequence to generate
// Outputs: a random DNA sequence of the specified length
std::string generate_random_sequence(size_t length) {
    const char nucleotides[] = {'A', 'C', 'G', 'T'};
    std::string seq;
    seq.reserve(length);
    for (size_t i = 0; i < length; ++i) {
        seq += nucleotides[rand() % 4];
    }
    return seq;
}

int main() {
    std::ofstream log_file("alignment_results.log");

    // Predefined test cases
    run_test("AGCTGAC", "GCTGAT", 4, 2, 1, -1, AlignmentMode::GLOBAL, log_file);
    run_test("GATTACA", "GCATGCU", 4, 1, 2, -1, AlignmentMode::GLOBAL, log_file);
    run_test("ACGT", "ACGT", 4, 1, 2, -1, AlignmentMode::GLOBAL, log_file);
    run_test("AAAA", "TTTT", 4, 1, 2, -1, AlignmentMode::GLOBAL, log_file);
    run_test("A", "G", 4, 1, 2, -1, AlignmentMode::GLOBAL, log_file);
    run_test("ACGTACGTACGT", "ACGTACGT", 4, 1, 2, -1, AlignmentMode::GLOBAL, log_file);
    run_test("ACGTTGCACGTA", "ACGTCACGTG", 4, 1, 2, -1, AlignmentMode::GLOBAL, log_file);

    // Randomized tests with various lengths
    for (int i = 0; i < 100; ++i) {
        std::string seq1 = generate_random_sequence(10 + rand() % 100);
        std::string seq2 = generate_random_sequence(10 + rand() % 100);
        run_test(seq1, seq2, 4, 2, 1, -1, AlignmentMode::GLOBAL, log_file);
    }

    log_file.close();
    std::cout << "All tests completed. Results written to alignment_results.log" << std::endl;

    return 0;
}
