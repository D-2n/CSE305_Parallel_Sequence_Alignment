#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <algorithm>

// Enum for traceback direction in the alignment matrix
enum class Direction {
    MATCH,       // Represents a match or mismatch in the alignment
    INSERTION,   // Represents an insertion in the alignment
    DELETION,    // Represents a deletion in the alignment
    NONE         // Represents no alignment direction (used for initialization)
};

// Struct to represent a block of the alignment matrix to be processed by a worker thread
struct MatrixBlock {
    unsigned int startRow; // Starting row of the block
    unsigned int startCol; // Starting column of the block
};

// Class for parallel sequence alignment focusing on partitioning (Section 5)
class ParallelSequenceAligner {
public:
    char* seq1; // Sequence 1
    char* seq2; // Sequence 2
    unsigned int len1, len2; // Lengths of sequences
    unsigned int num_threads; // Number of threads

    std::vector<std::thread> workers; // Worker threads
    std::mutex* mtx; // Mutex array for thread synchronization
    std::condition_variable* cv; // Condition variable array for thread synchronization

    std::atomic_int completed_threads; // Counter for completed threads
    std::atomic_int current_phase; // Current phase of computation
    int total_phases; // Total number of phases

    std::vector<std::vector<int>> score_matrix; // Alignment score matrix
    std::vector<std::vector<Direction>> traceback_matrix; // Traceback matrix

    // Constructor to initialize the aligner with sequences and parameters
    ParallelSequenceAligner(
        char* seq1, char* seq2,
        unsigned int len1, unsigned int len2,
        unsigned int num_threads)
        : seq1(seq1), seq2(seq2), len1(len1), len2(len2), num_threads(num_threads) {
        mtx = new std::mutex[num_threads];
        cv = new std::condition_variable[num_threads];

        score_matrix.resize(len1 + 1, std::vector<int>(len2 + 1, 0));
        traceback_matrix.resize(len1 + 1, std::vector<Direction>(len2 + 1, Direction::NONE));

        completed_threads.store(0);
        current_phase.store(1);

        total_phases = (len1 + len2) - 1;
    }

    // Destructor to clean up allocated resources
    ~ParallelSequenceAligner() {
        delete[] mtx;
        delete[] cv;
    }

    // Function to start the partitioning process using worker threads
    void find_partition() {
        for (unsigned int i = 0; i < num_threads; ++i) {
            workers.emplace_back(&ParallelSequenceAligner::worker_task, this, i + 1);
        }

        for (auto& worker : workers) {
            if (worker.joinable())
                worker.join();
        }

        // Output the partition points, traceback matrix, and score matrix for use in final alignment
        std::cout << "Partition Points and Traceback Matrix:\n";
        for (unsigned int i = 0; i <= len1; ++i) {
            for (unsigned int j = 0; j <= len2; ++j) {
                if (traceback_matrix[i][j] != Direction::NONE) {
                    std::cout << "Partition Point: (" << i << ", " << j << ") - "
                              << "Score: " << score_matrix[i][j] << ", "
                              << "Direction: " << static_cast<int>(traceback_matrix[i][j]) << "\n";
                }
            }
        }
    }

    // Function executed by each worker thread
    void worker_task(unsigned int worker_id) {
        int local_phase = 1;
        while (local_phase <= total_phases) {
            std::vector<MatrixBlock> blocks_to_compute;
            assign_blocks(worker_id, blocks_to_compute, local_phase);

            for (const auto& block : blocks_to_compute) {
                for (unsigned int i = block.startRow; i < std::min(block.startRow + 1, len1 + 1); ++i) {
                    for (unsigned int j = block.startCol; j < std::min(block.startCol + 1, len2 + 1); ++j) {
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
    void assign_blocks(unsigned int worker_id, std::vector<MatrixBlock>& blocks, unsigned int phase) {
        if (worker_id > num_threads) {
            std::cerr << "Invalid: worker_id is greater than the number of threads." << std::endl;
            return;
        }

        int i = 1 + (worker_id - 1);
        int j = 1 + (phase - worker_id);

        while (i < static_cast<int>(len1) + 1 && j >= 1) {
            blocks.push_back({static_cast<unsigned int>(i), static_cast<unsigned int>(j)});
            i += num_threads;
            j -= num_threads;
        }

        if (i < static_cast<int>(len1) + 1 && j >= 1) {
            blocks.push_back({static_cast<unsigned int>(i), static_cast<unsigned int>(j)});
        }
    }

    // Function to compute the score for a single cell in the alignment matrix
    void compute_cell(unsigned int i, unsigned int j) {
        int match_mismatch = score_matrix[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? 1 : -1);
        int insertion = score_matrix[i - 1][j] - 1;
        int deletion = score_matrix[i][j - 1] - 1;

        score_matrix[i][j] = std::max({0, match_mismatch, insertion, deletion});

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
};

// Runs multiple tests with different sequences.
// For each sequence pair, it initializes the ParallelSequenceAligner and finds the partition.
void run_tests() {
    std::vector<std::pair<std::string, std::string>> test_sequences = {
        {"AGCTGAC", "GCTGAT"},
        {"GATTACA", "GCATGCU"},
        {"ACGT", "ACGT"},
        {"AAAA", "TTTT"},
        {"A", "G"},
        {"ACGTACGTACGT", "ACGTACGT"},
        {"ACGTTGCACGTA", "ACGTCACGTG"}
    };

    for (const auto& seq_pair : test_sequences) {
        unsigned int len1 = seq_pair.first.size();
        unsigned int len2 = seq_pair.second.size();
        unsigned int num_threads = 4;

        ParallelSequenceAligner aligner(const_cast<char*>(seq_pair.first.c_str()), const_cast<char*>(seq_pair.second.c_str()), len1, len2, num_threads);
        aligner.find_partition();
    }
}

int main() {
    run_tests();
    return 0;
}
