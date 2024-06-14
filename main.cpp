#include "alignment_algorithm/main_alignment.h"
#include "test_functions/testing.h"
#include "test_functions/read_test_data.h"
#include <vector>
#include <string>
int main(){
    std::string filename = "./gene_sequences_test"; 

    //initialize data structs
    std::vector<std::string> names;
    std::vector<std::string> sequences;

    //import data
    read_and_store_sequences(names,sequences,filename);
    
    test_input_size_thread(names,sequences);
    
   // test_n_cores_thread(names,sequences);
   // test_similarity(names, sequences);
    return 0;
}