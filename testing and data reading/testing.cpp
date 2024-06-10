#include "read_test_data.h"
#include <iostream>
#include <random>
#include <chrono>
#include <thread>
#include <fstream>
/*
Testing section:
1) Test on sequences of different sizes.
2) Test depending on number of cores.
3) Test based on the similarity of input sequences.
*/

/**
test_input_size - Performs a set number of test batches, each of them iterating over sequences pre-determined variable size.
* @param names vector of gene sequence names
* @param sequences vector of gene sequences
* @return 0 if executed correctly
*/
int test_input_size( std::vector<std::string> &names, std::vector<std::string> &sequences){
    std::cout << "Testing with different input sizes\n\n";
    int batches = 10;
    int input_size_increment = 1000;

    int sequence_one_index;
    int sequence_two_index;

    size_t dataset_size = sequences.size();
    int range = dataset_size - 1;
    for (int t = 0; t < batches; t++){
        std::cout << "Testing batch " << t << "\n\n";

        //get two random sequences
        sequence_one_index = rand() % range;
        sequence_two_index = rand() % range;

        while (sequence_two_index == sequence_one_index){ //avoid duplicates here
            sequence_two_index = rand() % range;
        }
        size_t sequence_size = std::min(sequences[sequence_one_index].size(),sequences[sequence_two_index].size());
        
        std::cout << "Testing sequences \n" << names[sequence_one_index] << "\n and \n" << names[sequence_two_index] << "\n\n";
        for (size_t i = 1; i * input_size_increment < sequence_size; i++){
            std::cout << "Testing with input size: " << i * input_size_increment << "\n";
            auto start = std::chrono::high_resolution_clock::now();

            std::this_thread::sleep_for(std::chrono::seconds(2)); //replace with main function call when
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;

            std::cout << "Execution time: " << elapsed.count() << " seconds\n";
        }
    }
    return 0;
}

/**
test_n_cores - Performs a set number of test batches, each of them iterating over sequences while changing the number of available cores.
* @param names vector of gene sequence names
* @param sequences vector of gene sequences
* @return 0 if executed correctly
*/
int test_n_cores( std::vector<std::string> &names, std::vector<std::string> &sequences){
    std::cout << "Testing with different number of cores\n\n";
    int batches = 10;
    int n_tests = 5; //number of tests per batch, each with diff number of cores
    int sequence_one_index;
    int sequence_two_index;
    int n_cores;
    size_t dataset_size = sequences.size();
    int range = dataset_size - 1;
    for (int t = 0; t < batches; t++){
        std::cout << "Testing batch " << t << "\n\n";

        //get two random sequences
        sequence_one_index = rand() % range;
        sequence_two_index = rand() % range;

        while (sequence_two_index == sequence_one_index){ //avoid duplicates here
            sequence_two_index = rand() % range;
        }
        
        std::cout << "Testing sequences \n" << names[sequence_one_index] << "\n and \n" << names[sequence_two_index] << "\n\n";

        for (int i = 1; i <= n_tests; i++){
            n_cores = rand(); //change to appropriate when ready
            std::cout <<"("<< i <<"/" << n_tests << ") " <<"Testing with number of cores: " << n_cores << "\n";
            auto start = std::chrono::high_resolution_clock::now();

            std::this_thread::sleep_for(std::chrono::seconds(2)); //replace with main function call when
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;

            std::cout << "Execution time: " << elapsed.count() << " seconds\n";
        }
    }
    return 0;
}

/**
test_similarity - Performs a set number of test batches, each of them picking two sequences from the dataset, computing their similarity,
and then performing alignment.
* @param names vector of gene sequence names
* @param sequences vector of gene sequences
* @return 0 if executed correctly
*/
int test_similarity( std::vector<std::string> &names, std::vector<std::string> &sequences){
    std::cout << "Testing with similarity computation\n\n";
    int test_pairs = 20;
    int sequence_one_index;
    int sequence_two_index;
    size_t dataset_size = sequences.size();
    int range = dataset_size - 1;

    
    std::vector<std::thread> threads;
    std::vector<double> similarity_results(test_pairs);
    std::vector<std::chrono::duration<double>> execution_times(test_pairs);



    std::ofstream similarity_file;
    similarity_file.open ("similarity_testing.csv");

    similarity_file << "Testing with similarity computation\n";
    similarity_file << "Test number,Similarity,Execution time\n";

    
    auto sequence_similarity_threaded = [&](const std::string &seq1, const std::string &seq2, size_t id) {
        auto start = std::chrono::high_resolution_clock::now();
        
        std::cout << similarity_results[id] << "\n";
        similarity_results[id] = sequence_similarity(seq1, seq2);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        execution_times[id] = elapsed;
    };

    for (int t = 0; t < test_pairs; t++){
        std::cout << "Testing pair " << t+1 << "\n\n";
        //get two random sequences
        sequence_one_index = rand() % range;
        sequence_two_index = rand() % range;

        while (sequence_two_index == sequence_one_index){ //avoid duplicates here
            sequence_two_index = rand() % range;
        }
        threads.emplace_back(sequence_similarity_threaded, std::ref(sequences[sequence_one_index]), std::ref(sequences[sequence_two_index]), t);
        
    }
    
    for (auto& t : threads) {
        t.join();
    }
    //write results into .csv
    for (size_t j = 0; j < test_pairs; j++){
        std::cout << similarity_results[j] << "\n";
        similarity_file << j << "," << similarity_results[j] << "," << execution_times[j].count() << "\n";
    }
    similarity_file.close();
    return 0;
}
int main(){
    std::string filename = "./gene_sequences"; 

    //initialize data structs
    std::vector<std::string> names;
    std::vector<std::string> sequences;

    //import data
    read_and_store_sequences(names,sequences,filename);

    //test_input_size(names,sequences);
    
    //test_n_cores(names,sequences);
    test_similarity(names, sequences);

}
