#include "subproblem_alignment.h"
#include <deque>

typedef struct parallel_prefix_max_queue_element {
    long int value;
    size_t begin_id;
    size_t end_id;
    struct parallel_prefix_max_queue_element* next;
} queue_indices_max;

void PrefixMaxMapThread(std::vector<long int> &sums, long int value, queue_indices_max *curr) {
    if (curr != NULL) {
        for (size_t i = curr->begin_id; i < curr->end_id; i++) {
            sums[i] = std::max(sums[i], value);
        }
    }
}

void PrefixMaxInitMapThread(std::vector<long int> &values, std::vector<long int> &sums, queue_indices_max &q) {
    sums[q.begin_id] = values[q.begin_id];
    for (size_t i = q.begin_id + 1; i < q.end_id; i++) {
        sums[i] = std::max(sums[i-1], values[i]);
    }
    q.value = sums[q.end_id - 1];
}

void ParallelPrefixMax(size_t p, std::vector<long int> &values, std::vector<long int> &partial_sums) {

    std::vector<queue_indices_max> units;
    std::deque<queue_indices_max*> processes;

    // initialization of the elements
    size_t n = values.size();
    size_t block_size = n / p;
    size_t num_partitions = p;
    if (p > n) {
        block_size = 1;
        num_partitions = n;
    }
    queue_indices_max q;
    q.value = 0;
    q.begin_id = 0;
    q.next = NULL;
    q.end_id = block_size;
    units.push_back(q);
    for (size_t i = 1; i < num_partitions; i++) {
        queue_indices_max q;
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
        workers[i-1] = std::thread(&PrefixMaxInitMapThread, std::ref(values), std::ref(partial_sums), std::ref(units[i]));
    }
    PrefixMaxInitMapThread(values, partial_sums, units[0]);
    for (size_t i = 0; i < num_partitions - 1; i++) {
        processes.push_back(&(units[i]));
        workers[i].join();
    }
    processes.push_back(&(units[num_partitions - 1]));

    // pointer-jumping
    while (!processes.empty()) {
        size_t m = processes.size();
        std::vector<std::thread> workers(0);
        std::deque<queue_indices_max*> aux(m-1);
        for (size_t i = 0; i < m-1; i++) {
            queue_indices_max *q = processes.front();
            processes.pop_front();
            if(q->next != NULL) {
                queue_indices_max* r = q->next;
                q->next = r->next;
                q->value = partial_sums[q->end_id - 1];
                workers.push_back(std::thread(&PrefixMaxMapThread, std::ref(partial_sums), q->value, r));
                processes.push_back(q);
            }
        }
        
        queue_indices_max *q = processes.front();
        processes.pop_front();
        PrefixMaxMapThread(partial_sums, q->value, q->next);

        for(size_t i = 0; i < workers.size(); i++) {
            workers[i].join();
        }
        
    }

}

int main() {
    std::vector<long int> values = {-1, 3, 5, 1, 2, 7};
    std::vector<long int> partial_sums(6);
    ParallelPrefixMax(7, values, partial_sums);
    for(auto &p: partial_sums) {
        printf("%ld ", p);
    }
    printf("\n");
    return 0;
}

void Subproblem::ComputeFirstRowMapThread(Subproblem *subp, size_t start, size_t end) {
    while (start < end) {
        subp->T1[0][start] = std::numeric_limits<long int>::min();
        subp->T3[0][start] = std::numeric_limits<long int>::min();
        if (subp->start_type == -2) {
            subp->T2[0][start] = -subp->g * start;
        }
        else if (subp->start_type == 1 || subp->start_type == 3) {
            subp->T2[0][start] = std::numeric_limits<long int>::min();
        }
        else {
            subp->T2[0][start] = -subp->h - subp->g * start;
        }
        start++;
    }
}

void Subproblem::ComputeRowMapThread13(Subproblem *subp, size_t i, size_t start, size_t end) {
    while (start < end) {
        subp->T1[i][start] = subp->f(i, start) + std::max(std::max(subp->T1[i-1][start-1], subp->T2[i-1][start-1]), subp->T3[i-1][start-1]);
        subp->T3[i][start] = std::max(std::max(subp->T1[i-1][start] - subp->g - subp->h, subp->T2[i-1][start] - subp->g - subp->h), subp->T3[i-1][start] - subp->g);
    }
    start++;
}

void Subproblem::ComputeOmegaMapThread(Subproblem *subp, size_t i, size_t start, size_t end, std::vector<long int> &omega) {
    while (start < end) {
        omega[start] = start * subp->g + std::max(subp->T1[i][start-1] - subp->g - subp->h, subp->T3[i][start-1] - subp->g - subp->h);
        start++;
    }
}

void Subproblem::ComputeRowMapThread2(Subproblem *subp, size_t i, size_t start, size_t end, std::vector<long int> &partial) {
    while (start < end) {
        subp->T2[i][start] = partial[start] - start * subp->g;
        start++;
    }
}

void Subproblem::compute_row(size_t i) {
    size_t block_size = n / p;
    size_t num_threads = p;
    if (p > n) {
        num_threads = n;
        block_size = 1;
    }
    std::vector<std::thread> workers(num_threads-1);
    if (i == 0) {
        // top row
        T1[0][0] = std::numeric_limits<long int>::min();
        T2[0][0] = std::numeric_limits<long int>::min();
        T3[0][0] = std::numeric_limits<long int>::min();
        if (start_type == 1 || start_type == -1) {
            T1[0][0] = 0;
        }
        else if (start_type == -2) {
            T2[0][0] = 0;
        }
        else if (start_type == -3) {
            T3[0][0] = 0;
        }
        for (int j = 0; j < num_threads - 1; j++) {
            workers[j] = std::thread(&ComputeFirstRowMapThread, this, 1 + j*block_size, 1 + (j + 1) * block_size);
        }
        ComputeFirstRowMapThread(this, 1 + (num_threads-1)*block_size, n+1);
        for (int j = 0; j < num_threads - 1; j++) {
            workers[j].join();
        }
    }
    else {
        T1[i][0] = std::numeric_limits<long int>::min();
        T2[i][0] = std::numeric_limits<long int>::min();
        if (start_type == -3) {
            T3[i][0] = -g * i;
        }
        else if (start_type == 1 || start_type == 2) {
            T3[i][0] = std::numeric_limits<long int>::min();
        }
        else {
            T3[i][0] = -h - g * i;
        }
        // compute T1 and T3
        for (int j = 0; j < num_threads - 1; j++) {
            workers[j] = std::thread(&ComputeRowMapThread13, this, i, 1 + j*block_size, 1 + (j + 1) * block_size);
        }
        ComputeRowMapThread13(this, i, 1 + (num_threads-1)*block_size, n+1);
        for (int j = 0; j < num_threads - 1; j++) {
            workers[j].join();
        }
        // compute T2 using parallel prefix

        // compute omega +jg values
        std::vector<long int> omega(n+1);
        omega[0] = T2[i][0];
        for (int j = 0; j < num_threads - 1; j++) {
            workers[j] = std::thread(&ComputeOmegaMapThread, this, i, 1 + j*block_size, 1 + (j + 1) * block_size, std::ref(omega));
        }
        ComputeOmegaMapThread(this, i, 1 + (num_threads-1)*block_size, n+1, omega);
        for (int j = 0; j < num_threads - 1; j++) {
            workers[j].join();
        }

        // parallel prefix
        std::vector<long int> partials(n+1);
        ParallelPrefixMax(this->p, omega, partials);

        // deduce T2
        for (int j = 0; j < num_threads - 1; j++) {
            workers[j] = std::thread(&ComputeRowMapThread2, this, i, 1 + j*block_size, 1 + (j + 1) * block_size, std::ref(partials));
        }
        ComputeRowMapThread2(this, i, 1 + (num_threads-1)*block_size, n+1, partials);
        for (int j = 0; j < num_threads - 1; j++) {
            workers[j].join();
        }
    }
}