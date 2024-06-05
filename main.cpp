#include <stdio.h>
#include <thread>
#include <vector>
#include <math.h>

typedef struct alignment_point {
    int i;
    int j;
    int t;
    struct alignment_point* next = NULL;
} align;

void OptimalAlignmentMapThread(){
    // to complete
}

void optimal_alignment(std::vector<align> partial_bp, int n, int m, int p) {
    std::vector<int> omega(0);
    std::vector<int> partial_sums(0);
    std::vector<align>::iterator begin = partial_bp.begin();
    std::vector<align>::iterator end = partial_bp.end();
    int i = 0;
    while (begin + 1 != end) {
        int a = std::ceil(((begin+1)->i - begin->i) * p * 1.0 / m);
        int b = std::ceil(((begin+1)->j - begin->j) * p * 1.0 / n);
        omega.push_back(std::max(a, b));
        begin++;
        if (i == 0){
            partial_sums.push_back(omega[0]);
        }
        else{
            partial_sums.push_back(omega[i] + partial_sums[i-1]);
        }
        i++;
    }
    // add parallel prefix for partial sum

    /*for(int j=0; j<omega.size();j++){
        printf("%d ", omega[j]);
    }
    printf("\n");
    for(int j=0; j<partial_sums.size();j++){
        printf("%d ", partial_sums[j]);
    }
    printf("\n");*/
    
    // solve subproblems 0, 3, ...
    // solve subproblems 1, 4, ...
    // solve subproblems 2, 5, ...
    int num_suproblems = omega.size();
    std::vector<int> num_proc_per_sub(0); // to be computed after parallel-prefix
    std::vector<std::thread> workers(num_subproblems - 1);
    std::vector<align*> pointers (num_subproblems * 2); // beginning-end

    for (int i=0; i<num_subproblems - 1; i+=3) {
        workers[i] = std::thread(&OptimalAlignmentMapThread, ...)
    }
    OptimalAlignmentMapThread(...);
    for(int i=0; i<num_subproblems - 1; i++) {
        workers[i].join();
    }
    for (int i=1; i<num_subproblems - 1; i+=3) {
        workers[i] = std::thread(&OptimalAlignmentMapThread, ...)
    }    
    OptimalAlignmentMapThread(...);
    for(int i=0; i<num_subproblems - 1; i++) {
        workers[i].join();
    }
    for (int i=2; i<num_subproblems - 1; i+=3) {
        workers[i] = std::thread(&OptimalAlignmentMapThread, ...)
    }
    OptimalAlignmentMapThread(...);
    for(int i=0; i<num_subproblems - 1; i++) {
        workers[i].join();
    }
    // concatenate the results
    for(int i = 1; i<num_subproblems-1; i++) {
        align* prev = pointers[2*(i-1)+1];
        align* next = pointers[2*i];
        prev -> next = next;
    }
}

int main(){
    std::vector<align> pbp(0);
    align s1, s2, s3, s4, s5;
    s1.i = 0; s1.j = 0; s1.t = -1;
    s2.i = 8; s2.j = 7; s2.t = -1;
    s3.i = 16; s3.j = 33; s3.t = -1;
    s4.i = 40; s4.j = 36; s4.t = -2;
    s5.i = 48; s5.j = 48; s5.t = 1;
    pbp.push_back(s1);
    pbp.push_back(s2);
    pbp.push_back(s3);
    pbp.push_back(s4);
    pbp.push_back(s5);
    optimal_alignment(pbp, 48, 48, 6);

    return 0;
}