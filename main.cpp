#include <stdio.h>
#include <thread>
#include <vector>
#include <math.h>
#include <deque>

typedef struct alignment_point {
    size_t i;
    size_t j;
    size_t t;
    struct alignment_point* next = NULL;
} align;

void OptimalAlignmentMapThread(){
    // to compete
}

typedef struct parallel_prefix_queue_element {
    size_t value;
    size_t begin_id;
    size_t end_id;
    struct parallel_prefix_queue_element* next;
} queue_indices;

void PrefixSumMapThread(std::vector<size_t> &sums, size_t value, queue_indices *curr) {
    if (curr != NULL) {
        for (size_t i = curr->begin_id; i < curr->end_id; i++) {
            sums[i] += value;
        }
    }
}

void PrefixInitMapThread(std::vector<size_t> &values, std::vector<size_t> &sums, queue_indices &q) {
    sums[q.begin_id] = values[q.begin_id];
    for (size_t i = q.begin_id + 1; i < q.end_id; i++) {
        sums[i] = sums[i-1] + values[i];
    }
    q.value = sums[q.end_id - 1];
}

void ParallelPrefix(size_t p, std::vector<size_t> &values, std::vector<size_t> &partial_sums) {

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

void ComputeOmegaMapThread(std::vector<align>::iterator begin, std::vector<align>::iterator end, size_t n, size_t m, size_t p, std::vector<size_t> &omega, size_t offset) {
    size_t i = 0;
    while (begin + 1 != end) {
        size_t a = std::ceil(((begin+1)->i - begin->i) * p * 1.0 / m);
        size_t b = std::ceil(((begin+1)->j - begin->j) * p * 1.0 / n);
        omega[i + offset] = std::max(a, b);
        begin++;
        i++;
    }
}

void compute_omega_parallel(std::vector<align> &partial_bp, size_t n, size_t m, size_t p, size_t len, std::vector<size_t> &omega) {
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
        workers[i] = std::thread(&ComputeOmegaMapThread, start_block, end_block, n, m, p, std::ref(omega), i*block_size);
        start_block += block_size;
    }
    ComputeOmegaMapThread(start_block, end, n, m, p, omega, (num_threads-1) * block_size);
    
    for(size_t i = 0; i < num_threads - 1; i++) {
        workers[i].join();
    }
}

void optimal_alignment(std::vector<align> partial_bp, size_t n, size_t m, size_t p) {
    size_t num_subproblems = partial_bp.size() - 1;
    std::vector<size_t> omega(num_subproblems);
    std::vector<size_t> partial_sums(num_subproblems);

    //compute omega
    compute_omega_parallel(partial_bp, n, m, p, num_subproblems, omega);

    // add parallel prefix for partial sum
    ParallelPrefix(p, omega, partial_sums);

    /* testing */
    for(size_t j=0; j<omega.size();j++){
        printf("%d ", omega[j]);
    }
    printf("\n");
    for(size_t j=0; j<partial_sums.size();j++){
        printf("%d ", partial_sums[j]);
    }
    printf("\n");
    
    
    /* solving subproblems
    // solve subproblems 0, 3, ...
    // solve subproblems 1, 4, ...
    // solve subproblems 2, 5, ...
    size_t num_subproblems = omega.size();
    std::vector<size_t> num_proc_per_sub(0); // to be computed after parallel-prefix
    std::vector<std::thread> workers(num_subproblems - 1);
    std::vector<align*> pointers (num_subproblems * 2); // beginning-end

    for (size_t i=0; i<num_subproblems - 1; i+=3) {
        workers[i] = std::thread(&OptimalAlignmentMapThread, ...)
    }
    OptimalAlignmentMapThread(...);
    for(size_t i=0; i<num_subproblems - 1; i++) {
        workers[i].join();
    }
    for (size_t i=1; i<num_subproblems - 1; i+=3) {
        workers[i] = std::thread(&OptimalAlignmentMapThread, ...)
    }    
    OptimalAlignmentMapThread(...);
    for(size_t i=0; i<num_subproblems - 1; i++) {
        workers[i].join();
    }
    for (size_t i=2; i<num_subproblems - 1; i+=3) {
        workers[i] = std::thread(&OptimalAlignmentMapThread, ...)
    }
    OptimalAlignmentMapThread(...);
    for(size_t i=0; i<num_subproblems - 1; i++) {
        workers[i].join();
    }
    // concatenate the results
    for(size_t i = 1; i<num_subproblems-1; i++) {
        align* prev = pointers[2*(i-1)+1];
        align* next = pointers[2*i];
        prev -> next = next;
    }
    */
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