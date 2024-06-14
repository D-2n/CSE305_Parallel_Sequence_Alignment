#pragma once

#ifndef PARTIAL_H
#define PARTIAL_H

#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <algorithm>
#include <cmath>
#include <climits>

// Define the alignment_point struct
typedef struct alignment_point {
    size_t i;
    size_t j;
    int t;
    struct alignment_point* next = nullptr;
} align;

// Function prototypes
int score(char a, char b);

void initializeTables(std::vector<std::vector<int>>& T1, std::vector<std::vector<int>>& T2, std::vector<std::vector<int>>& T3, size_t m, size_t n, double g, double h, int start_type);

void initializeReverseTables(std::vector<std::vector<int>>& TR1, std::vector<std::vector<int>>& TR2, std::vector<std::vector<int>>& TR3, size_t m, size_t n, double g, double h, int end_type);

void fillTablesParallel(const char* A, const char* B, size_t m, size_t n, 
    std::vector<std::vector<int>>& T1, std::vector<std::vector<int>>& T2, 
    std::vector<std::vector<int>>& T3, double g, double h, size_t p);

void fillReverseTablesParallel(const char* A, const char* B, size_t m, size_t n, 
    std::vector<std::vector<int>>& TR1, std::vector<std::vector<int>>& TR2, 
    std::vector<std::vector<int>>& TR3, double g, double h, size_t p);

std::vector<align> findPartitionParallel(const std::vector<std::vector<int>>& T1, const std::vector<std::vector<int>>& T2, const std::vector<std::vector<int>>& T3, 
    const std::vector<std::vector<int>>& TR1, const std::vector<std::vector<int>>& TR2, const std::vector<std::vector<int>>& TR3, 
    size_t m, size_t n, size_t p, double h);

void findPartialBalancedPartitionParallel(const char* A, const char* B, size_t m, size_t n, size_t p, double g, double h, int start_type, int end_type, std::vector<align>& partition);

void extractPartitions(const std::vector<align>& partition, const char* A, const char* B);

#endif // PARTIAL_H
