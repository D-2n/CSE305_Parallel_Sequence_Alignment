#include <stdio.h>
#include <thread>
#include <vector>
#include <math.h>
#include <deque>
#include "subproblem_alignment.h"
#include "../sequence_alignment/partial.h"



void OptimalAlignmentMapThread(char *A, char *B, size_t m, size_t n, size_t ida, size_t idb, size_t p, int start_type, int end_type, double g, double h, align* &begin, align* &end){
    printf("bp1\n");
    Subproblem subp = Subproblem(A, B, m, n, ida, idb, p, start_type, end_type, g, h);
    printf("bp1.2\n");
    subp.compute_tables();
    printf("bp2\n");
    subp.find_alignment();
    printf("bp3\n");
    begin = subp.alignment_begin;
    end = subp.alignment_end;
    printf("bp4\n");
}



void print_align(align *begin) {
    while (begin != NULL) {
        printf("(%ld, %ld, %d)\n", begin->i, begin->j, begin->t);
        begin = begin -> next;
    }
}
void print_seq(char *A, char *B, align *begin) {
    align *begin_a = begin;
    align *begin_b = begin;
    while (begin_a != NULL) {
        if (begin_a->t == 1 || begin_a->t == 3) {
            printf("%c", A[begin_a->i]);
        }
        else {
            printf("-");
        }
        begin_a = begin_a -> next;
    }
    printf("\n");
    while (begin_b != NULL) {
        if (begin_b->t == 1 ||  begin_b->t == 2) {
            printf("%c", B[begin->j]);
        }
        else {
            printf("-");
        }
        begin_b = begin_b -> next;
    }
    printf("\n");
}


typedef struct parallel_prefix_queue_element {
    long int value;
    size_t begin_id;
    size_t end_id;
    struct parallel_prefix_queue_element* next;
} queue_indices;

void PrefixSumMapThread(std::vector<long int> &sums, long int value, queue_indices *curr) {
    if (curr != NULL) {
        for (size_t i = curr->begin_id; i < curr->end_id; i++) {
            sums[i] += value;
        }
    }
}

void PrefixInitMapThread(std::vector<long int> &values, std::vector<long int> &sums, queue_indices &q) {
    sums[q.begin_id] = values[q.begin_id];
    for (size_t i = q.begin_id + 1; i < q.end_id; i++) {
        sums[i] = sums[i-1] + values[i];
    }
    q.value = sums[q.end_id - 1];
}

void ParallelPrefix(size_t p, std::vector<long int> &values, std::vector<long int> &partial_sums) {

    std::vector<queue_indices> units;
    std::deque<queue_indices*> processes;

    // initialization of the elements
    size_t n = values.size();
    size_t block_size = n / p;
    size_t num_partitions = p;
    if (p > n) {
        block_size = 1;
        num_partitions = n;
    }
    queue_indices q;
    q.value = 0;
    q.begin_id = 0;
    q.next = NULL;
    q.end_id = block_size;
    units.push_back(q);
    for (size_t i = 1; i < num_partitions; i++) {
        queue_indices q;
        q.value = 0;
        q.begin_id = block_size * i;
        q.next = NULL;
        units[i-1].next = &q;
        q.end_id = block_size * (i + 1);
        if (q.end_id >= n || i == num_partitions - 1) {
            q.end_id = n;
        }
        units.push_back(q);
        units[i-1].next = &(units[i]);
    }


    // suming the elements of p partitions

    std::vector<std::thread> workers(num_partitions-1);
    for (size_t i = 1; i < num_partitions; i++) {
        units[i-1].next = &(units[i]);
        workers[i-1] = std::thread(&PrefixInitMapThread, std::ref(values), std::ref(partial_sums), std::ref(units[i]));
    }
    PrefixInitMapThread(values, partial_sums, units[0]);
    for (size_t i = 0; i < num_partitions - 1; i++) {
        processes.push_back(&(units[i]));
        workers[i].join();
    }
    processes.push_back(&(units[num_partitions - 1]));

    // pointer-jumping
    while (!processes.empty()) {
        size_t m = processes.size();
        std::vector<std::thread> workers(0);
        std::deque<queue_indices*> aux(m-1);
        for (size_t i = 0; i < m-1; i++) {
            queue_indices *q = processes.front();
            processes.pop_front();
            if(q->next != NULL) {
                queue_indices* r = q->next;
                q->next = r->next;
                q->value = partial_sums[q->end_id - 1];
                workers.push_back(std::thread(&PrefixSumMapThread, std::ref(partial_sums), q->value, r));
                processes.push_back(q);
            }
        }
        
        queue_indices *q = processes.front();
        processes.pop_front();
        PrefixSumMapThread(partial_sums, q->value, q->next);

        for(size_t i = 0; i < workers.size(); i++) {
            workers[i].join();
        }
        
    }

}

void ComputeOmegaMapThread(std::vector<align>::iterator begin, std::vector<align>::iterator end, size_t m, size_t n, size_t p, std::vector<long int> &omega, size_t offset) {
    size_t i = 0;
    while (begin + 1 != end) {
        long int a = std::ceil(1.0 * ((begin+1)->i - begin->i) / ((1.0 * m) / p));
        long int b = std::ceil(1.0 * ((begin+1)->j - begin->j) / ((1.0 * n) / p));
        omega[i + offset] = std::max(a, b);
        begin++;
        i++;
    }
}

void compute_omega_parallel(std::vector<align> &partial_bp, size_t m, size_t n, size_t p, size_t len, std::vector<long int> &omega) {
    size_t block_size = len / p;
    size_t num_threads = p;
    if (p > len) {
        block_size = 1;
        num_threads = len;
    }
    std::vector<std::thread> workers(num_threads-1);
    std::vector<align>::iterator begin = partial_bp.begin();
    std::vector<align>::iterator end = partial_bp.end();
    std::vector<align>::iterator start_block = begin;
    for(size_t i = 0; i < num_threads - 1; i++) {
        std::vector<align>::iterator end_block = start_block + block_size + 1;
        workers[i] = std::thread(&ComputeOmegaMapThread, start_block, end_block, m, n, p, std::ref(omega), i*block_size);
        start_block += block_size;
    }
    ComputeOmegaMapThread(start_block, end, m, n, p, omega, (num_threads-1) * block_size);
    
    for(size_t i = 0; i < num_threads - 1; i++) {
        workers[i].join();
    }
}

size_t assign_processors(long int sum_prev, long int curr_subproblem) {
    if (sum_prev % 3 == 0) {
        return (curr_subproblem + 2) / 3;
    }
    else if (sum_prev % 3 == 1) {
        return 1 + curr_subproblem / 3;
    }
    return 1 + (curr_subproblem + 1) / 3;
}

void optimal_alignment(char *A, char *B, std::vector<align> partial_bp, size_t m, size_t n, size_t p, double g, double h) {
    size_t num_subproblems = partial_bp.size() - 1;
    std::vector<long int> omega(num_subproblems);
    std::vector<long int> partial_sums(num_subproblems);

    //compute omega
    compute_omega_parallel(partial_bp, m, n, p, num_subproblems, omega);

    // add parallel prefix for partial sum
    ParallelPrefix(p, omega, partial_sums);

    /* testing */
    /*
    for(size_t j=0; j<omega.size();j++){
        printf("%ld ", omega[j]);
    }
    printf("\n");
    for(size_t j=0; j<partial_sums.size();j++){
        printf("%ld ", partial_sums[j]);
    }
    printf("\n");

    printf("%d ", assign_processors(0, omega[0]));
    for(size_t j=1; j<partial_sums.size();j++){
        printf("%d ", assign_processors(partial_sums[j-1], omega[j]));
    }
    printf("\n");
    */
    
    /* solving subproblems*/
    std::vector<std::thread> workers(num_subproblems - 1);
    std::vector<align*> pointers (num_subproblems * 2); // beginning-end

    size_t i = 0;
    // solve subproblems 0, 3, ...
    while (i < num_subproblems - 3) {
        size_t lenA = partial_bp[i+1].i - partial_bp[i].i;
        size_t lenB = partial_bp[i+1].j - partial_bp[i].j;
        size_t ida = partial_bp[i].i;
        size_t idb = partial_bp[i].j;
        size_t p_prime;
        if (i == 0) {
            p_prime = assign_processors(0, omega[0]);
        }
        else {
            p_prime = assign_processors(partial_sums[i-1], omega[i]);
        }
        int start_type = partial_bp[i].t;
        int end_type = -partial_bp[i+1].t;
        
        workers[i] = std::thread(&OptimalAlignmentMapThread, A, B, lenA, lenB, ida, idb, p_prime, start_type, end_type, g, h, std::ref(pointers[2*i]), std::ref(pointers[2*i+1]));
        
        i += 3;
    }
    size_t lenA, lenB, ida, idb, p_prime;
    int start_type, end_type;
    if (i < num_subproblems) {
        lenA= partial_bp[i+1].i - partial_bp[i].i;
        lenB = partial_bp[i+1].j - partial_bp[i].j;
        ida = partial_bp[i].i;
        idb = partial_bp[i].j;
        if (i == 0) {
            p_prime = assign_processors(0, omega[0]);
        }
        else {
            p_prime = assign_processors(partial_sums[i-1], omega[i]);
        }
        start_type = partial_bp[i].t;
        end_type = -partial_bp[i+1].t;
        OptimalAlignmentMapThread(A, B, lenA, lenB, ida, idb, p_prime, start_type, end_type, g, h, pointers[2*i], pointers[2*i+1]);
    
    }
    for(size_t i=0; i<num_subproblems - 3; i += 3) {
        workers[i].join();
    }
    

    // solve subproblems 1, 4, ...
    i = 1;
    while (i < num_subproblems - 3) {
        size_t lenA = partial_bp[i+1].i - partial_bp[i].i;
        size_t lenB = partial_bp[i+1].j - partial_bp[i].j;
        size_t ida = partial_bp[i].i;
        size_t idb = partial_bp[i].j;
        size_t p_prime = assign_processors(partial_sums[i-1], omega[i]);
        int start_type = partial_bp[i].t;
        int end_type = -partial_bp[i+1].t;
        workers[i] = std::thread(&OptimalAlignmentMapThread, A, B, lenA, lenB, ida, idb, p_prime, start_type, end_type, g, h, std::ref(pointers[2*i]), std::ref(pointers[2*i+1]));
        
        i += 3;
    }
    if (i < num_subproblems) {
        lenA= partial_bp[i+1].i - partial_bp[i].i;
        lenB = partial_bp[i+1].j - partial_bp[i].j;
        ida = partial_bp[i].i;
        idb = partial_bp[i].j;
        p_prime = assign_processors(partial_sums[i-1], omega[i]);
        start_type = partial_bp[i].t;
        end_type = -partial_bp[i+1].t;
        OptimalAlignmentMapThread(A, B, lenA, lenB, ida, idb, p_prime, start_type, end_type, g, h, pointers[2*i], pointers[2*i+1]);
    
    }
    for(size_t i=1; i<num_subproblems - 3; i += 3) {
        workers[i].join();
    }

    // solve subproblems 2, 5, ...
    i = 2;
    while (i < num_subproblems - 3) {
        size_t lenA = partial_bp[i+1].i - partial_bp[i].i;
        size_t lenB = partial_bp[i+1].j - partial_bp[i].j;
        size_t ida = partial_bp[i].i;
        size_t idb = partial_bp[i].j;
        size_t p_prime = assign_processors(partial_sums[i-1], omega[i]);
        int start_type = partial_bp[i].t;
        int end_type = -partial_bp[i+1].t;
        workers[i] = std::thread(&OptimalAlignmentMapThread, A, B, lenA, lenB, ida, idb, p_prime, start_type, end_type, g, h, std::ref(pointers[2*i]), std::ref(pointers[2*i+1]));
        i += 3;
    }
    if (i < num_subproblems) {
        lenA= partial_bp[i+1].i - partial_bp[i].i;
        lenB = partial_bp[i+1].j - partial_bp[i].j;
        ida = partial_bp[i].i;
        idb = partial_bp[i].j;
        p_prime = assign_processors(partial_sums[i-1], omega[i]);
        start_type = partial_bp[i].t;
        end_type = -partial_bp[i+1].t;
        OptimalAlignmentMapThread(A, B, lenA, lenB, ida, idb, p_prime, start_type, end_type, g, h, pointers[2*i], pointers[2*i+1]);
    
    }
    for(size_t i=2; i<num_subproblems - 3; i += 3) {
        workers[i].join();
    }

    // concatenate the results
    for(size_t i = 1; i<num_subproblems-1; i++) {
        align* prev = pointers[2*(i-1)+1];
        align* next = pointers[2*i];
        prev -> next = next;
    }
  // print_align(pointers[0]);
    print_seq(A,B,pointers[0]);
}

int main_alignment_function(char* A, char* B, size_t m, size_t n, size_t p, double g, double h){
    std::cout << m << "\n";

    A = "-AGGA";
    B = "-AGTGC";
    m = 4;
    n = 5;
    p = 32;
    g = 1;
    h = 2;
    
    std::vector<align> partitions;
    findPartialBalancedPartitionParallel(A, B, m, n, p, g, h, partitions);
    printf("Found partitions\n");
    std::cout << partitions.size() << "\n";
    // to comment when pbp is done
    /*
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
    */
    // uncomment when pbp is done
    // optimal_partition(A, B, m, n, pbp, p);
    
    /*
    align s1, s2, s3, s4;
    s1.i = 0; s1.j = 0; s1.t = -1;
    s2.i = 0; s2.j = 2; s2.t = 2;
    s3.i = 3; s3.j = 4; s3.t = 1;
    s4.i = m; s4.j = n; s4.t = 1;
    pbp.push_back(s1);
    pbp.push_back(s2);
    pbp.push_back(s3);
    pbp.push_back(s4);
    */

    optimal_alignment(A, B, partitions, m, n, p, g, h);

    return 0;
}
