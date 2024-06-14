#include "subproblem_alignment.h"
#include <deque>
#include <cstdlib>


typedef struct parallel_prefix_max_queue_element {
    double value;
    size_t begin_id;
    size_t end_id;
    struct parallel_prefix_max_queue_element* next;
} queue_indices_max;

void PrefixMaxMapThread(std::vector<double> &sums, double value, queue_indices_max *curr) {
    if (curr != NULL) {
        for (size_t i = curr->begin_id; i < curr->end_id; i++) {
            sums[i] = std::max(sums[i], value);
        }
    }
}

void PrefixMaxInitMapThread(std::vector<double> &values, std::vector<double> &sums, queue_indices_max &q) {
    sums[q.begin_id] = values[q.begin_id];
    for (size_t i = q.begin_id + 1; i < q.end_id; i++) {
        sums[i] = std::max(sums[i-1], values[i]);
    }
    q.value = sums[q.end_id - 1];
}

void ParallelPrefixMax(size_t p, std::vector<double> &values, std::vector<double> &partial_sums) {

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

void Subproblem::find_alignment() {
    align *curr_point, *new_point;
    curr_point = (align*)std::malloc(sizeof(align));
    size_t i = m;
    size_t j = n;
    curr_point->next = NULL;
    alignment_end = curr_point;
    if (end_type > 0) {
        curr_point->t = end_type;
        if (end_type == 1) {
            curr_point->i = i + id_A;
            curr_point->j = j + id_B;
        }
        else if (end_type == 2) {
            curr_point->i = 0;
            curr_point->j = j + id_B;
        }
        else {
            curr_point->i = i + id_A;
            curr_point->j = 0;
        }
    }
    else {
        double t1 = T1[m][n];
        double t2 = T2[m][n] + h_prime(-2);
        double t3 = T3[m][n] + h_prime(-3);
        if (t1 >= t2 && t1 >= t3) {
            curr_point->t = 1;
            curr_point->i = i + id_A;
            curr_point->j = j + id_B;
        }
        else if (t2 >= t1 && t2 >= t3) {
            curr_point->t = 2;
            curr_point->i = 0;
            curr_point->j = j + id_B;
        }
        else {
            curr_point->t = 3;
            curr_point->i = i + id_A;
            curr_point->j = 0;
        }
    }
    while (i > 0 && j > 0) {
        //printf("(i, j, t): %ld, %ld, %d\n", curr_point->i, curr_point->j, curr_point->t);
        new_point = (align*)std::malloc(sizeof(align));
        if (curr_point->t == 1) {
            if (T1[i][j] == f(i, j) + T1[i-1][j-1]) {new_point->t = 1; new_point->i = i-1 + id_A; new_point->j = j-1 + id_A; i--; j--;}
            else if (T1[i][j] == f(i, j) + T2[i-1][j-1]) {new_point->t = 2; new_point->i = 0; new_point->j = j-1 + id_B; i--; j--;}
            else if (T1[i][j] == f(i, j) + T3[i-1][j-1]) {new_point->t = 3; new_point->i = i-1 + id_A; new_point->j = 0; i--; j--;}

        }
        else if (curr_point->t == 2) {
            if (T2[i][j] == -g - h + T1[i][j-1]) {new_point->t = 1; new_point->i = i + id_A; new_point->j = j-1 + id_B; j--;}
            else if (T2[i][j] == -g + T2[i][j-1]) {new_point->t = 2; new_point->i = 0; new_point->j = j-1 + id_B; j--;}
            else if (T2[i][j] == -g - h + T3[i][j-1]) {new_point->t = 3; new_point->i = i + id_A; new_point->j = 0; j--;}

        }
        else {
            if (T3[i][j] == -g - h + T1[i-1][j]) {new_point->t = 1; new_point->i = i-1 + id_A; new_point->j = j + id_B; i--;}
            else if (T3[i][j] == -g - h + T2[i-1][j]) {new_point->t = 2; new_point->i = 0; new_point->j = j + id_B; i--;}
            else if (T3[i][j] == -g + T3[i-1][j]) {new_point->t = 3; new_point->i = i-1 + id_A; new_point->j = 0; i--;}
        }
        new_point->next = curr_point;
        curr_point = new_point;
    }
    alignment_begin = curr_point->next;
    //print_alignment();
}

void Subproblem::print_alignment() {
    align* begin = alignment_begin;
    while (begin != NULL) {
        printf("(%ld, %ld, %d)\n", begin->i, begin->j, begin->t);
        begin = begin -> next;
    }
}
/*TEST
int main() {
    /*std::vector<double> values = {-1, 3, 5, 1, 2, 7};
    std::vector<double> partial_sums(6);
    ParallelPrefixMax(7, values, partial_sums);
    for(auto &p: partial_sums) {
        printf("%lf ", p);
    }
    printf("\n"); 
    char* A = "-AGGA";
    char* B = "-ATGTC";
    size_t m = 4;
    size_t n = 5;
    size_t ida = 0;
    size_t idb = 0;
    size_t p = 3;
    int start_type = -1;
    int end_type = -1;
    double g = 2;
    double h = 1;
    Subproblem subp = Subproblem(A, B, m, n, ida, idb, p, start_type, end_type, g, h);
    //printf("parallel:\n");
    subp.compute_tables();
    /*
    printf("non-parallel:\n");
    subp.non_parallel_tables();
    subp.find_alignment();

    return 0;
}
*/
void Subproblem::ComputeFirstRowMapThread(Subproblem *subp, size_t start, size_t end) {
    while (start < end) {
        subp->T1[0][start] = -std::numeric_limits<double>::infinity();
        subp->T3[0][start] = -std::numeric_limits<double>::infinity();
        if (subp->start_type == -2) {
            subp->T2[0][start] = -subp->g * start;
        }
        else if (subp->start_type == 1 || subp->start_type == 3) {
            subp->T2[0][start] = -std::numeric_limits<double>::infinity();
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
        start++;
    }
}

void Subproblem::ComputeOmegaMapThread(Subproblem *subp, size_t i, size_t start, size_t end, std::vector<double> &omega) {
    while (start < end) {
        omega[start] = start * subp->g + std::max(subp->T1[i][start-1] - subp->g - subp->h, subp->T3[i][start-1] - subp->g - subp->h);
        start++;
    }
}

void Subproblem::ComputeRowMapThread2(Subproblem *subp, size_t i, size_t start, size_t end, std::vector<double> &partial) {
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
        T1[0][0] = -std::numeric_limits<double>::infinity();
        T2[0][0] = -std::numeric_limits<double>::infinity();
        T3[0][0] = -std::numeric_limits<double>::infinity();
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
        T1[i][0] = -std::numeric_limits<double>::infinity();
        T2[i][0] = -std::numeric_limits<double>::infinity();
        if (start_type == -3) {
            T3[i][0] = -g * i;
        }
        else if (start_type == 1 || start_type == 2) {
            T3[i][0] = -std::numeric_limits<double>::infinity();
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
        std::vector<double> omega(n+1);
        omega[0] = T2[i][0];
        for (int j = 0; j < num_threads - 1; j++) {
            workers[j] = std::thread(&ComputeOmegaMapThread, this, i, 1 + j*block_size, 1 + (j + 1) * block_size, std::ref(omega));
        }
        ComputeOmegaMapThread(this, i, 1 + (num_threads-1)*block_size, n+1, omega);
        for (int j = 0; j < num_threads - 1; j++) {
            workers[j].join();
        }

        // parallel prefix
        std::vector<double> partials(n+1);
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

void Subproblem::compute_tables() {
    for (size_t i = 0; i <= m; i++) {
        compute_row(i);
    }
    /* printf("T1:\n");
    for (size_t i = 0; i<= m; i++) {
        for (size_t j = 0; j<= n; j++) {
            printf("%lf ", T1[i][j]);
        }
        printf("\n");
    }
    printf("T2:\n");
    for (size_t i = 0; i<= m; i++) {
        for (size_t j = 0; j<= n; j++) {
            printf("%lf ", T2[i][j]);
        }
        printf("\n");
    }
    printf("T3:\n");
    for (size_t i = 0; i<= m; i++) {
        for (size_t j = 0; j<= n; j++) {
            printf("%lf ", T3[i][j]);
        }
        printf("\n");
    } */

}

void Subproblem::non_parallel_tables() {
    T1[0][0] = -std::numeric_limits<double>::infinity();
    T2[0][0] = -std::numeric_limits<double>::infinity();
    T3[0][0] = -std::numeric_limits<double>::infinity();
    if (start_type == 1 || start_type == -1) {
        T1[0][0] = 0;
    }
    else if (start_type == -2) {
        T2[0][0] = 0;
    }
    else if (start_type == -3) {
        T3[0][0] = 0;
    }
    for (size_t j = 1; j <= n; j++) {
        T1[0][j] = -std::numeric_limits<double>::infinity();
        T3[0][j] = -std::numeric_limits<double>::infinity();
        if (start_type == -2) {
            T2[0][j] = -g * j;
        }
        else if (start_type == 1 || start_type == 3) {
            T2[0][j] = -std::numeric_limits<double>::infinity();
        }
        else {
            T2[0][j] = -h - g * j;
        }
    }
    for (size_t i = 1; i <= m; i++) {
        T1[i][0] = -std::numeric_limits<double>::infinity();
        T2[i][0] = -std::numeric_limits<double>::infinity();
        if (start_type == -3) {
            T3[i][0] = -g * i;
        }
        else if (start_type == 1 || start_type == 2) {
            T3[i][0] = -std::numeric_limits<double>::infinity();
        }
        else {
            T3[i][0] = -h - g * i;
        }
        for (size_t j = 1; j <= n; j++) {
            T1[i][j] = f(i, j) + std::max(std::max(T1[i-1][j-1], T2[i-1][j-1]), T3[i-1][j-1]);
            T3[i][j] = std::max(std::max(T1[i-1][j] - g - h, T2[i-1][j] - g - h), T3[i-1][j] - g);
            T2[i][j] = std::max(std::max(T1[i][j-1] -g-h, T2[i][j-1] - g), T3[i][j-1] - g-h); 
        }
    }
    printf("T1:\n");
    for (size_t i = 0; i<= m; i++) {
        for (size_t j = 0; j<= n; j++) {
            printf("%lf ", T1[i][j]);
        }
        printf("\n");
    }
    printf("T2:\n");
    for (size_t i = 0; i<= m; i++) {
        for (size_t j = 0; j<= n; j++) {
            printf("%lf ", T2[i][j]);
        }
        printf("\n");
    }
    printf("T3:\n");
    for (size_t i = 0; i<= m; i++) {
        for (size_t j = 0; j<= n; j++) {
            printf("%lf ", T3[i][j]);
        }
        printf("\n");
    }
}