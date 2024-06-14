#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

struct align {
    size_t i;
    size_t j;
    int t;
    align* next = nullptr;
};

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

void traceback_path(const vector<vector<double>>& T1, const vector<vector<double>>& T2, const vector<vector<double>>& T3, vector<align>& partial_bp) {
    int m = T1.size() - 1;
    int n = T1[0].size() - 1;
    int i = m, j = n;

    while (i > 0 && j > 0) {
        if (T1[i][j] >= T2[i][j] && T1[i][j] >= T3[i][j]) {
            size_t i_ = i;
            size_t j_ = j;
            partial_bp.push_back({i_, j_, 1});
            --i;
            --j;
        } else if (T2[i][j] >= T1[i][j] && T2[i][j] >= T3[i][j]) {
            size_t i_ = i;
            size_t j_ = j;
            partial_bp.push_back({i_, j_, 2});
            --j;
        } else {
            size_t i_ = i;
            size_t j_ = j;
            partial_bp.push_back({i_, j_, 3});
            --i;
        }
    }
    partial_bp.push_back({0, 0, -1}); // Add start point
    reverse(partial_bp.begin(), partial_bp.end());
}

void optimal_partition(const string& A, const string& B, vector<align>& partial_bp) {
    int m = A.length();
    int n = B.length();
    
    vector<vector<double>> T1, T2, T3;
    initialize_tables(T1, T2, T3, m, n);
    fill_tables(A, B, T1, T2, T3);
    traceback_path(T1, T2, T3, partial_bp);
}

int main() {
    string A = "ATGTCGA";
    string B = "AGAATCTA";

    vector<align> partial_bp;
    optimal_partition(A, B, partial_bp);

    for (const auto& point : partial_bp) {
        cout << "(" << point.i << ", " << point.j << ", " << point.t << ")\n";
    }

    return 0;
}
