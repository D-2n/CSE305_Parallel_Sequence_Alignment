#pragma once
#ifndef main_alignment_H
#define main_alignment_H

#include <stdio.h>
#include <thread>
#include <vector>
#include <math.h>
#include <deque>
#include "subproblem_alignment.h"
typedef struct parallel_prefix_queue_element {
    long int value;
    size_t begin_id;
    size_t end_id;
    struct parallel_prefix_queue_element* next;
} queue_indices;
void OptimalAlignmentMapThread(char *A, char *B, size_t m, size_t n, size_t ida, size_t idb, size_t p, int start_type, int end_type, double g, double h, align* &begin, align* &end);

void print_align(align *begin);


void PrefixSumMapThread(std::vector<long int> &sums, long int value, queue_indices *curr);


void PrefixInitMapThread(std::vector<long int> &values, std::vector<long int> &sums, queue_indices &q);

void ParallelPrefix(size_t p, std::vector<long int> &values, std::vector<long int> &partial_sums);

void ComputeOmegaMapThread(std::vector<align>::iterator begin, std::vector<align>::iterator end, size_t n, size_t m, size_t p, std::vector<long int> &omega, long int offset);


void compute_omega_parallel(std::vector<align> &partial_bp, size_t n, size_t m, size_t p, size_t len, std::vector<long int> &omega);

size_t assign_processors(long int sum_prev, long int curr_subproblem);

void optimal_alignment(char *A, char *B, std::vector<align> partial_bp, size_t n, size_t m, size_t p, double g, double h);

int main_alignment_function(char* A, char* B, size_t m, size_t n, size_t p, double g, double h);

#endif
