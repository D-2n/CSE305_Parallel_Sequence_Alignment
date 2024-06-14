#include <vector>
#include <algorithm>
#include <omp.h>
#include <ostream>
#include <climits>
#include <cstring>
#include <iostream>

typedef struct alignment_point {
    size_t i;
    size_t j;
    int t;
    struct alignment_point* next = NULL;
} align;

int score(char a, char b) {
    return (a == b) ? 0 : 1;
}

void initializeTables(std::vector<std::vector<int>>& T, size_t m, size_t n) {
    for (size_t i = 0; i <= m; ++i) {
        T[i][0] = 0;
    }
    for (size_t j = 0; j <= n; ++j) {
        T[0][j] = 0;
    }
}

void fillTablesParallel(const char* A, const char* B, size_t m, size_t n, 
    std::vector<std::vector<int>>& T1, std::vector<std::vector<int>>& T2, 
    std::vector<std::vector<int>>& T3, double g, double h, size_t p) {
    
    #pragma omp parallel for num_threads(p)
    for (size_t i = 1; i <= m; ++i) {
        for (size_t j = 1; j <= n; ++j) {
            T1[i][j] = std::max({T1[i-1][j-1] + score(A[i-1], B[j-1]), T2[i-1][j-1] + score(A[i-1], B[j-1]), T3[i-1][j-1] + score(A[i-1], B[j-1])});
            T2[i][j] = std::max({T1[i][j-1] - (int)(g + h), T2[i][j-1] - (int)g, T3[i][j-1] - (int)(g + h)});
            T3[i][j] = std::max({T1[i-1][j] - (int)(g + h), T2[i-1][j] - (int)(g + h), T3[i-1][j] - (int)g});
        }
    }
}

void fillReverseTablesParallel(const char* A, const char* B, size_t m, size_t n, 
    std::vector<std::vector<int>>& TR1, std::vector<std::vector<int>>& TR2, 
    std::vector<std::vector<int>>& TR3, double g, double h, size_t p) {
    
    #pragma omp parallel for num_threads(p)
    for (size_t i = m; i > 0; --i) {
        for (size_t j = n; j > 0; --j) {
            TR1[i][j] = std::max({TR1[i+1][j+1] + score(A[i-1], B[j-1]), TR2[i+1][j+1] + score(A[i-1], B[j-1]), TR3[i+1][j+1] + score(A[i-1], B[j-1])});
            TR2[i][j] = std::max({TR1[i][j+1] - (int)(g + h), TR2[i][j+1] - (int)g, TR3[i][j+1] - (int)(g + h)});
            TR3[i][j] = std::max({TR1[i+1][j] - (int)(g + h), TR2[i+1][j] - (int)(g + h), TR3[i+1][j] - (int)g});
        }
    }
}

std::vector<align> findPartitionParallel(const std::vector<std::vector<int>>& T1, const std::vector<std::vector<int>>& T2, const std::vector<std::vector<int>>& T3, 
    const std::vector<std::vector<int>>& TR1, const std::vector<std::vector<int>>& TR2, const std::vector<std::vector<int>>& TR3, 
    size_t m, size_t n, size_t p, double h) {
    
    std::vector<align> partition;
    size_t block_size_m = m / p;
    size_t block_size_n = n / p;

    partition.push_back({0, 0, -1, nullptr});

    #pragma omp parallel for num_threads(p)
    for (size_t k = 1; k < p; ++k) {
        int max_val = INT_MIN;
        align best_align = {0, 0, 0, nullptr};

        // Find max value in the special row
        for (size_t i = k * block_size_m; i < (k + 1) * block_size_m && i <= m; ++i) {
            for (size_t j = 1; j <= n; ++j) {
                int val = std::max({
                    T1[i][j] + TR1[i][j],
                    T2[i][j] + TR2[i][j] + (int)h,
                    T3[i][j] + TR3[i][j] + (int)h
                });
                if (val > max_val) {
                    max_val = val;
                    best_align = {i, j, (T1[i][j] >= T2[i][j] && T1[i][j] >= T3[i][j]) ? 1 : (T2[i][j] >= T3[i][j]) ? 2 : 3, nullptr};
                }
            }
        }
        #pragma omp critical
        partition.push_back(best_align);

        max_val = INT_MIN;
        best_align = {0, 0, 0, nullptr};

        // Find max value in the special column
        for (size_t j = k * block_size_n; j < (k + 1) * block_size_n && j <= n; ++j) {
            for (size_t i = 1; i <= m; ++i) {
                int val = std::max({
                    T1[i][j] + TR1[i][j],
                    T2[i][j] + TR2[i][j] + (int)h,
                    T3[i][j] + TR3[i][j] + (int)h
                });
                if (val > max_val) {
                    max_val = val;
                    best_align = {i, j, (T1[i][j] >= T2[i][j] && T1[i][j] >= T3[i][j]) ? 1 : (T2[i][j] >= T3[i][j]) ? 2 : 3, nullptr};
                }
            }
        }
        #pragma omp critical
        partition.push_back(best_align);
    }

    partition.push_back({m, n, 1, nullptr});

    std::sort(partition.begin(), partition.end(), [](const align& a, const align& b) {
        return (a.i < b.i) || (a.i == b.i && a.j < b.j);
    });

    return partition;
}

void findPartialBalancedPartitionParallel(const char* A, const char* B, size_t m, size_t n, size_t p, double g, double h, std::vector<align>& partition) {
    std::vector<std::vector<int>> T1(m + 1, std::vector<int>(n + 1));
    std::vector<std::vector<int>> T2(m + 1, std::vector<int>(n + 1));
    std::vector<std::vector<int>> T3(m + 1, std::vector<int>(n + 1));
    std::vector<std::vector<int>> TR1(m + 2, std::vector<int>(n + 2, 0));
    std::vector<std::vector<int>> TR2(m + 2, std::vector<int>(n + 2, 0));
    std::vector<std::vector<int>> TR3(m + 2, std::vector<int>(n + 2, 0));

    initializeTables(T1, m, n);
    initializeTables(T2, m, n);
    initializeTables(T3, m, n);
    initializeTables(TR1, m+1, n+1);
    initializeTables(TR2, m+1, n+1);
    initializeTables(TR3, m+1, n+1);

    fillTablesParallel(A, B, m, n, T1, T2, T3, g, h, p);
    fillReverseTablesParallel(A, B, m, n, TR1, TR2, TR3, g, h, p);
    partition = findPartitionParallel(T1, T2, T3, TR1, TR2, TR3, m, n, p, h);
}

