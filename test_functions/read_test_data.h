
#ifndef read_test_data_H
#define read_test_data_H

#include <vector>
#include <string>

int read_and_store_sequences( std::vector<std::string>& names, std::vector<std::string>& sequences, std::string& filename);

double sequence_similarity(const std::string &sequence1, const std::string &sequence2);
#endif 