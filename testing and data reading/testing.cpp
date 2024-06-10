#include "read_test_data.h"
#include <iostream>
#include <random>
#include <chrono>
#include <thread>
#include <fstream>
#include <list>
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
int test_input_size_thread(  std::vector<std::string> &names, std::vector<std::string> &sequences){
    //get the size of the smallest list
    std::cout << "Testing with different input sizes\n\n";
    size_t available_threads =  std::thread::hardware_concurrency();
    int test_pairs = 2000;
    int input_size_inc = 2000;
    int sequence_one_index;
    int sequence_two_index;

    size_t chunk_size = (test_pairs + available_threads - 1) / available_threads; // Correct chunk size calculation
    size_t n_chunks = available_threads;
    size_t remainder = test_pairs % chunk_size;


    size_t dataset_size = sequences.size();
    int range = dataset_size - 1;

    
    std::vector<std::thread> threads;
    std::vector<double> input_sizes(test_pairs);
    std::vector<std::chrono::duration<double>> execution_times(test_pairs);



    std::ofstream test_res_file;
    test_res_file.open ("input_size_testing.csv");

    test_res_file << "Testing with different input sizes\n";
    test_res_file << "Test number,Input size,Execution time\n";

    
    auto test_input_threaded = [&](const std::vector<std::string> &sequences, size_t start, size_t end, size_t input_size) {
        for (size_t i = start; i < end && i < test_pairs; ++i){
            
            sequence_one_index = rand() % range;
            sequence_two_index = rand() % range;
            std::string seq1 = sequences[sequence_one_index];
            std::string seq2 = sequences[sequence_two_index];

            auto start = std::chrono::high_resolution_clock::now();
            //perform call, also do input size = max(minimal length of the two sequences, provided input size)
            input_sizes[i] = input_size;
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            
            execution_times[i] = elapsed;
        }
    };
    size_t input_size;
    for (size_t t = 0; t < n_chunks; t++){
        //get two random sequences
        size_t start = t * chunk_size;
        size_t end = (t == available_threads - 1) ? (test_pairs) : (start + chunk_size); 
        input_size = (((t+1) * chunk_size) / 150) * input_size_inc;
        threads.emplace_back(test_input_threaded, std::ref(sequences), start, end, input_size);
        
    }
    
    for (auto& t : threads) {
        t.join();
    }
    //write results into .csv
    for (size_t j = 0; j < test_pairs; j++){
        test_res_file << j << "," << input_sizes[j] << "," << execution_times[j].count() << "\n";
    }
    test_res_file.close();
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
    size_t available_threads =  std::thread::hardware_concurrency();
    int test_pairs = 2000;
    int sequence_one_index;
    int sequence_two_index;

    size_t chunk_size = (test_pairs + available_threads - 1) / available_threads; // Correct chunk size calculation
    size_t n_chunks = available_threads;
    size_t remainder = test_pairs % chunk_size;


    size_t dataset_size = sequences.size();
    int range = dataset_size - 1;

    
    std::vector<std::thread> threads;
    std::vector<double> similarity_results(test_pairs);
    std::vector<std::chrono::duration<double>> execution_times(test_pairs);



    std::ofstream similarity_file;
    similarity_file.open ("similarity_testing.csv");

    similarity_file << "Testing with similarity computation\n";
    similarity_file << "Test number,Similarity,Execution time\n";

    
    auto sequence_similarity_threaded = [&](const std::vector<std::string> &sequences, size_t start, size_t end) {
        for (size_t i = start; i < end && i < test_pairs; ++i){
            
            sequence_one_index = rand() % range;
            sequence_two_index = rand() % range;
            std::string seq1 = sequences[sequence_one_index];
            std::string seq2 = sequences[sequence_two_index];

            auto start = std::chrono::high_resolution_clock::now();
            
            similarity_results[i] = sequence_similarity(std::ref(seq1), std::ref(seq2));
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;
            
            execution_times[i] = elapsed;
        }
    };
    for (size_t t = 0; t < n_chunks; t++){
        //get two random sequences
        size_t start = t * chunk_size;
        size_t end = (t == available_threads - 1) ? (test_pairs) : (start + chunk_size); 
        threads.emplace_back(sequence_similarity_threaded, std::ref(sequences), start, end);
        
    }
    
    for (auto& t : threads) {
        t.join();
    }
    //write results into .csv
    for (size_t j = 0; j < test_pairs; j++){
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

    test_input_size_thread(names,sequences);
    
    //test_n_cores(names,sequences);
    //test_similarity(names, sequences);

}
