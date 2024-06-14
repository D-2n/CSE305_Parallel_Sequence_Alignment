#pragma once

#ifndef partial_H
#define partial_H

#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <algorithm>
#include <cmath>
#include <climits>

using namespace std;
int score(char a, char b);
void initializeTables(std::vector<std::vector<int>>& T, size_t m, size_t n);

void fillTablesParallel(const char* A, const char* B, size_t m, size_t n, 
    std::vector<std::vector<int>>& T1, std::vector<std::vector<int>>& T2, 
    std::vector<std::vector<int>>& T3, double g, double h, size_t p) ;

void fillReverseTablesParallel(const char* A, const char* B, size_t m, size_t n, 
    std::vector<std::vector<int>>& TR1, std::vector<std::vector<int>>& TR2, 
    std::vector<std::vector<int>>& TR3, double g, double h, size_t p);

std::vector<align> findPartitionParallel(const std::vector<std::vector<int>>& T1, const std::vector<std::vector<int>>& T2, const std::vector<std::vector<int>>& T3, 
    const std::vector<std::vector<int>>& TR1, const std::vector<std::vector<int>>& TR2, const std::vector<std::vector<int>>& TR3, 
    size_t m, size_t n, size_t p, double h);

void findPartialBalancedPartitionParallel(const char* A, const char* B, size_t m, size_t n, size_t p, double g, double h, std::vector<align>& partition);


#endif