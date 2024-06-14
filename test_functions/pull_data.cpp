#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_set>
#include "read_test_data.h"
#include <thread>
#include <atomic>
#include <algorithm>
/**
* Read gene sequences from a file and store them in a string vector, and store their names in another one.
* @param names reference to a list that stores names and identifiers of gene sequences.
* @param sequences reference to a list that stores gene sequences.
* @param filename path to the data file.
* @return 0 if executed sucessfully, 1 otherwise.
* @note Canis Lupus familiaris a.k.a "dog" genes taken from the NCBI database: https://www.ncbi.nlm.nih.gov/datasets/gene/taxon/9615/
*/
int read_and_store_sequences(std::vector<std::string>& names, std::vector<std::string>& sequences, std::string& filename){
    
    std::cout << "Opening data file: " << filename <<"\n";
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file! Check the file path/name!\n";
        return 1;
    }
    std::cout << "File opened successfully!\n";

    std::string line;
    std::string current_sequence;

    std::cout << "Storing sequences...\n";
    while (std::getline(file, line)) {
        if (line[0] == '>') { //lines of new sequences begin with >NC_006587.....
            if (!current_sequence.empty()) {
                sequences.push_back(current_sequence);
                current_sequence.clear();
            }
            names.push_back(line);
        } else {
            current_sequence += line;
        }
    }
    if (!current_sequence.empty()) {
        sequences.push_back(current_sequence);
    }

    file.close();
    if (sequences.size() != names.size()){
        std::cout << "Error: mismatch in sequences and names list sizes\n";
        return 1;
    }
    //Sanity check for duplicate sequences

    std::cout << "Checking for duplicate sequences...\n";
    std::unordered_set<std::string> sequence_set;
    bool has_duplicates = false;
    for (const auto& seq : sequences) {
        if (!sequence_set.insert(seq).second) {
            has_duplicates = true;
            std::cout << "Duplicate sequence found!\n";
        }
    }

    if (!has_duplicates) {
        std::cout << "No duplicate sequences found.\n";
    }else{
        std::cout << "There is at least one duplicate sequence found. Please check your data file.\n";
    }
    std::cout << "Dataset read successfully!\n";
    return 0;
}

/**
Compares characters in a chunk of two sequences.
* @param sequence1 First gene sequence for comparison
* @param sequence2 Second gene sequence for comparison
* @param start starting index of comparison
* @param end ending index of comparison
* @param score reference to the overall (atomic) score
*/
void compare_chars(const std::string &sequence1, const std::string &sequence2, size_t start, size_t end, std::atomic<int>& score) {
    int local_score = 0;
    for (size_t i = start; i < end; ++i){
        if (sequence1[i] == sequence2[i]) {
            local_score +=1;
        }
    }
    score.fetch_add(local_score,std::memory_order_relaxed);
}
/**
Computes similarity between two sequences: each matching pair of values awards one point, and the total amount is divided by the minimum length
of the two sequences. Threaded version.
* @param sequence1 First gene sequence for comparison
* @param sequence2 Second gene sequence for comparison
* @return similarity of two sequences, normalized between 0 and 1
*/
double sequence_similarity(const std::string &sequence1, const std::string &sequence2) {
    size_t l1 = sequence1.size();
    size_t l2 = sequence2.size();
    size_t iteration_size = std::min(l1, l2);
    size_t max_len = std::max(l1, l2);

    std::atomic<int> score(0);
    std::vector<std::thread> threads;

    size_t num_threads = std::thread::hardware_concurrency();;
    size_t chunk_size = iteration_size / num_threads;
    size_t n_chunks = iteration_size / chunk_size;
    size_t remainder = iteration_size % num_threads;
    
    for (size_t i = 0; i < n_chunks; ++i) {
        size_t start = i * chunk_size;
        size_t end = (i == num_threads - 1) ? (start + chunk_size + remainder) : (start + chunk_size); 
        threads.emplace_back(compare_chars, std::ref(sequence1), std::ref(sequence2), start, end, std::ref(score));
    }

    for (auto& t : threads) {
        if (t.joinable()) {
            t.join();
        }else{
            std::cout << "ERROR: Failed to join thread!";
        }
    }

    double similarity = static_cast<double>(score.load()) / max_len;
    return similarity;
}