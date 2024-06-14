#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <tuple>

using namespace std;

typedef struct alignment_point {
    size_t i;
    size_t j;
    int t;
    struct alignment_point* next = nullptr;
} align;

const double MATCH_SCORE = 2;
const double MISMATCH_PENALTY = 0;
const double GAP_OPEN_PENALTY = 2;
const double GAP_EXTENSION_PENALTY = 1;

double match_score(char a, char b) {
    return (a == b) ? MATCH_SCORE : MISMATCH_PENALTY;
}

void initialize_tables(vector<vector<double>>& T1, vector<vector<double>>& T2, vector<vector<double>>& T3, int m, int n) {
    T1.assign(m + 1, vector<double>(n + 1, -INFINITY));
    T2.assign(m + 1, vector<double>(n + 1, -INFINITY));
    T3.assign(m + 1, vector<double>(n + 1, -INFINITY));
}

void fill_tables(const string& A, const string& B, vector<vector<double>>& T1, vector<vector<double>>& T2, vector<vector<double>>& T3) {
    int m = A.length();
    int n = B.length();
    
    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j <= n; ++j) {
            if (i == 0 || j == 0) {
                T1[i][j] = T2[i][j] = T3[i][j] = 0;
            } else {
                T1[i][j] = match_score(A[i-1], B[j-1]) + max({T1[i-1][j-1], T2[i-1][j-1], T3[i-1][j-1]});
                T2[i][j] = max({T1[i][j-1] - GAP_OPEN_PENALTY, T2[i][j-1] - GAP_EXTENSION_PENALTY, T3[i][j-1] - GAP_OPEN_PENALTY});
                T3[i][j] = max({T1[i-1][j] - GAP_OPEN_PENALTY, T2[i-1][j] - GAP_OPEN_PENALTY, T3[i-1][j] - GAP_EXTENSION_PENALTY});
            }
        }
    }
}

void traceback_path(const vector<vector<double>>& T1, const vector<vector<double>>& T2, const vector<vector<double>>& T3, vector<align>& partial_bp, int start_i, int start_j) {
    int m = T1.size() - 1;
    int n = T1[0].size() - 1;
    int i = m, j = n;

    while (i > 0 && j > 0) {
        if (T1[i][j] >= T2[i][j] && T1[i][j] >= T3[i][j]) {
            size_t i_ = start_i + i;
            size_t j_ = start_j + j;
            partial_bp.push_back({i_, j_, 1});
            --i;
            --j;
        } else if (T2[i][j] >= T1[i][j] && T2[i][j] >= T3[i][j]) {
            size_t i_ = start_i + i;
            size_t j_ = start_j + j;
            partial_bp.push_back({i_, j_, 2});
            --j;
        } else {
            size_t i_ = start_i + i;
            size_t j_ = start_j + j;
            partial_bp.push_back({i_, j_, 3});
            --i;
        }
    }
    partial_bp.push_back({start_i, start_j, 0}); // Add start point with a special value for t
    reverse(partial_bp.begin(), partial_bp.end());
}

// This function computes the partition points for the sequences based on the number of processors.
vector<align> compute_partitions(const string& A, const string& B, int num_partitions) {
    int m = A.length();
    int n = B.length();
    vector<align> partitions;
    int step_i = m / num_partitions;
    int step_j = n / num_partitions;

    for (int t = 0; t < num_partitions; ++t) {
        int start_i = t * step_i;
        int end_i = (t == num_partitions - 1) ? m : (t + 1) * step_i - 1;
        int start_j = t * step_j;
        int end_j = (t == num_partitions - 1) ? n : (t + 1) * step_j - 1;
        
        partitions.push_back({start_i, start_j, 0});
        partitions.push_back({end_i, end_j, 0});
    }

    return partitions;
}

int main() {
    string A = "ATGTCGA";
    string B = "AGAATCTA";
    
    int num_partitions = 4; // Example number of partitions
    vector<align> partitions = compute_partitions(A, B, num_partitions);

    for (const auto& partition : partitions) {
        cout << "Partition: (" << partition.i << ", " << partition.j << ", " << partition.t << ")\n";
    }

    return 0;
}

