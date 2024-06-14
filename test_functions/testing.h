#pragma once
#ifndef testing_H
#define testing_H

#include <string>
#include <vector>
int test_input_size( std::vector<std::string> &names, std::vector<std::string> &sequences);

int test_input_size_thread(  std::vector<std::string> &names, std::vector<std::string> &sequences);

/**
test_n_cores - Performs a set number of test batches, each of them iterating over sequences while changing the number of available cores.
* @param names vector of gene sequence names
* @param sequences vector of gene sequences
* @return 0 if executed correctly
*/
int test_n_cores(std::vector<std::string> &names, std::vector<std::string> &sequences);
int test_n_cores_thread(  std::vector<std::string> &names, std::vector<std::string> &sequences);
/**
test_similarity - Performs a set number of test batches, each of them picking two sequences from the dataset, computing their similarity,
and then performing alignment.
* @param names vector of gene sequence names
* @param sequences vector of gene sequences
* @return 0 if executed correctly
*/
int test_similarity( std::vector<std::string> &names, std::vector<std::string> &sequences);
#endif