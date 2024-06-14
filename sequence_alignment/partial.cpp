#include <vector>
#include <algorithm>
#include <omp.h>
#include <iostream>
#include <climits>
#include <cstring>
#include "partial.h"

int score(char a, char b) {
    return (a == b) ? 0 : 1;
}

void initializeTables(std::vector<std::vector<int>>& T1, std::vector<std::vector<int>>& T2, std::vector<std::vector<int>>& T3, size_t m, size_t n, double g, double h, int start_type) {
    for (size_t i = 0; i <= m; ++i) {
        for (size_t j = 0; j <= n; ++j) {
            T1[i][j] = T2[i][j] = T3[i][j] = INT_MIN;  // Initialize all cells to -∞
        }
    }

    if (start_type == 1) {
        T1[0][0] = 0;
    } else if (start_type == 2) {
        for (size_t j = 1; j <= n; ++j) {
            T2[0][j] = -(int)(g + h) * j;
        }
    } else if (start_type == 3) {
        for (size_t i = 1; i <= m; ++i) {
            T3[i][0] = -(int)(g + h) * i;
        }
    }
}

void initializeReverseTables(std::vector<std::vector<int>>& TR1, std::vector<std::vector<int>>& TR2, std::vector<std::vector<int>>& TR3, size_t m, size_t n, double g, double h, int end_type) {
    for (size_t i = 0; i <= m + 1; ++i) {
        for (size_t j = 0; j <= n + 1; ++j) {
            TR1[i][j] = TR2[i][j] = TR3[i][j] = INT_MIN;  // Initialize all cells to -∞
        }
    }

    if (end_type == 1) {
        TR1[m+1][n+1] = 0;
    } else if (end_type == 2) {
        for (size_t j = n; j > 0; --j) {
            TR2[m+1][j] = -(int)(g + h) * (n - j + 1);
        }
    } else if (end_type == 3) {
        for (size_t i = m; i > 0; --i) {
            TR3[i][n+1] = -(int)(g + h) * (m - i + 1);
        }
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
    for (size_t i = m; i >= 1; --i) {
        for (size_t j = n; j >= 1; --j) {
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
        int max_val_row = INT_MIN;
        int max_val_col = INT_MIN;
        align best_align_row = {0, 0, 0, nullptr};
        align best_align_col = {0, 0, 0, nullptr};

        // Find max value in the special row
        for (size_t i = k * block_size_m; i < (k + 1) * block_size_m && i <= m; ++i) {
            for (size_t j = 1; j <= n; ++j) {
                int val = std::max({
                    T1[i][j] + TR1[i][j],
                    T2[i][j] + TR2[i][j] + (int)h,
                    T3[i][j] + TR3[i][j] + (int)h
                });
                if (val > max_val_row) {
                    max_val_row = val;
                    best_align_row = {i, j, (T1[i][j] >= T2[i][j] && T1[i][j] >= T3[i][j]) ? 1 : (T2[i][j] >= T3[i][j]) ? 2 : 3, nullptr};
                }
            }
        }

        // Find max value in the special column
        for (size_t j = k * block_size_n; j < (k + 1) * block_size_n && j <= n; ++j) {
            for (size_t i = 1; i <= m; ++i) {
                int val = std::max({
                    T1[i][j] + TR1[i][j],
                    T2[i][j] + TR2[i][j] + (int)h,
                    T3[i][j] + TR3[i][j] + (int)h
                });
                if (val > max_val_col) {
                    max_val_col = val;
                    best_align_col = {i, j, (T1[i][j] >= T2[i][j] && T1[i][j] >= T3[i][j]) ? 1 : (T2[i][j] >= T3[i][j]) ? 2 : 3, nullptr};
                }
            }
        }

        // Add the best align point to the partition list
        #pragma omp critical
        {
            if (max_val_row > max_val_col) {
                partition.push_back(best_align_row);
            } else {
                partition.push_back(best_align_col);
            }
        }
    }

    partition.push_back({m, n, 1, nullptr});

    std::sort(partition.begin(), partition.end(), [](const align& a, const align& b) {
        return (a.i < b.i) || (a.i == b.i && a.j < b.j);
    });

    return partition;
}


void findPartialBalancedPartitionParallel(const char* A, const char* B, size_t m, size_t n, size_t p, double g, double h, int start_type, int end_type, std::vector<align>& partition) {
    std::vector<std::vector<int>> T1(m + 1, std::vector<int>(n + 1));
    std::vector<std::vector<int>> T2(m + 1, std::vector<int>(n + 1));
    std::vector<std::vector<int>> T3(m + 1, std::vector<int>(n + 1));
    std::vector<std::vector<int>> TR1(m + 2, std::vector<int>(n + 2, INT_MIN));
    std::vector<std::vector<int>> TR2(m + 2, std::vector<int>(n + 2, INT_MIN));
    std::vector<std::vector<int>> TR3(m + 2, std::vector<int>(n + 2, INT_MIN));

    initializeTables(T1, T2, T3, m, n, g, h, start_type);
    initializeReverseTables(TR1, TR2, TR3, m, n, g, h, end_type);

    fillTablesParallel(A, B, m, n, T1, T2, T3, g, h, p);
    fillReverseTablesParallel(A, B, m, n, TR1, TR2, TR3, g, h, p);
    partition = findPartitionParallel(T1, T2, T3, TR1, TR2, TR3, m, n, p, h);
}

/*void extractPartitions(const std::vector<align>& partition, const char* A, const char* B) {
    for (size_t k = 0; k < partition.size() - 1; ++k) {
        size_t start_i = partition[k].i;
        size_t start_j = partition[k].j;
        size_t end_i = partition[k + 1].i;
        size_t end_j = partition[k + 1].j;

        std::cout << "Partition " << k + 1 << ":\n";
        std::cout << "A: ";
        for (size_t i = start_i; i < end_i; ++i) {
            std::cout << A[i];
        }
        std::cout << "\n";

        std::cout << "B: ";
        for (size_t j = start_j; j < end_j; ++j) {
            std::cout << B[j];
        }
        std::cout << "\n\n";
    }
}*/
