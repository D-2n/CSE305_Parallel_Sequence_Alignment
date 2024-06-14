#include <iostream>
#include <vector>
#include <thread>
#include <cmath>
#include <algorithm>

using namespace std;

typedef struct alignment_point {
    size_t i;
    size_t j;
    int t;
    struct alignment_point* next = NULL;
} align;

void fill_tables_parallel(char *A, char *B, size_t m, size_t n, vector<vector<int>> &T1, vector<vector<int>> &T2, vector<vector<int>> &T3, double g, double h, size_t p) {
    auto fill_chunk = [&](size_t start_row, size_t end_row) {
        for (size_t i = start_row; i <= end_row; i++) {
            for (size_t j = 1; j <= n; j++) {
                T1[i][j] = max({T1[i-1][j-1], T2[i-1][j-1], T3[i-1][j-1]}) + (A[i-1] == B[j-1] ? 0 : 1);
                T2[i][j] = max({T1[i][j-1] - g - h, T2[i][j-1] - g, T3[i][j-1] - g - h});
                T3[i][j] = max({T1[i-1][j] - g - h, T2[i-1][j] - g - h, T3[i-1][j] - g});
            }
        }
    };

    vector<thread> threads;
    size_t chunk_size = m / p;
    for (size_t i = 0; i < p; i++) {
        size_t start_row = i * chunk_size + 1;
        size_t end_row = (i == p - 1) ? m : (i + 1) * chunk_size;
        threads.emplace_back(fill_chunk, start_row, end_row);
    }

    for (auto &t : threads) {
        t.join();
    }
}

vector<align> find_partition_parallel(vector<vector<int>> &T1, vector<vector<int>> &T2, vector<vector<int>> &T3, size_t m, size_t n, size_t p) {
    vector<align> partition;

    size_t block_size_m = m / p;
    size_t block_size_n = n / p;

    partition.push_back({0, 0, -1, nullptr});

    auto find_max_in_block = [&](size_t i, align &max_align) {
        size_t row = i * block_size_m;
        size_t col = i * block_size_n;

        int max_score = -1;

        for (size_t j = 1; j < n; j++) {
            if (T1[row][j] > max_score) {
                max_score = T1[row][j];
                max_align = {row, j, 1, nullptr};
            }
            if (T2[row][j] > max_score) {
                max_score = T2[row][j];
                max_align = {row, j, 2, nullptr};
            }
            if (T3[row][j] > max_score) {
                max_score = T3[row][j];
                max_align = {row, j, 3, nullptr};
            }
        }

        for (size_t k = 1; k < m; k++) {
            if (T1[k][col] > max_score) {
                max_score = T1[k][col];
                max_align = {k, col, 1, nullptr};
            }
            if (T2[k][col] > max_score) {
                max_score = T2[k][col];
                max_align = {k, col, 2, nullptr};
            }
            if (T3[k][col] > max_score) {
                max_score = T3[k][col];
                max_align = {k, col, 3, nullptr};
            }
        }
    };

    vector<thread> threads;
    vector<align> max_aligns(p);

    for (size_t i = 1; i < p; i++) {
        threads.emplace_back(find_max_in_block, i, ref(max_aligns[i]));
    }

    for (auto &t : threads) {
        t.join();
    }

    for (const auto &max_align : max_aligns) {
        partition.push_back(max_align);
    }

    partition.push_back({m, n, 1, nullptr});
    return partition;
}

void find_partial_balanced_partition_parallel(char *A, char *B, size_t m, size_t n, size_t p, double g, double h, vector<align> &partition) {
    vector<vector<int>> T1(m + 1, vector<int>(n + 1));
    vector<vector<int>> T2(m + 1, vector<int>(n + 1));
    vector<vector<int>> T3(m + 1, vector<int>(n + 1));
    
    // Initialize the tables
    for (size_t i = 0; i <= m; i++) {
        T1[i][0] = -1;
        T2[i][0] = -1;
        T3[i][0] = -1;
    }
    for (size_t j = 0; j <= n; j++) {
        T1[0][j] = -1;
        T2[0][j] = -1;
        T3[0][j] = -1;
    }
    T1[0][0] = 0;

    fill_tables_parallel(A, B, m, n, T1, T2, T3, g, h, p);
    partition = find_partition_parallel(T1, T2, T3, m, n, p);
}

int main() {
    char A[] = "AGTACGCA";
    char B[] = "TATGC";
    size_t m = 8;
    size_t n = 5;
    size_t p = 4;
    double g = 1.0;
    double h = 2.0;
    
    vector<align> partition;
    find_partial_balanced_partition_parallel(A, B, m, n, p, g, h, partition);
    
    for (const auto &align : partition) {
        cout << "(" << align.i << ", " << align.j << ", " << align.t << ")\n";
    }
    
    return 0;
}

