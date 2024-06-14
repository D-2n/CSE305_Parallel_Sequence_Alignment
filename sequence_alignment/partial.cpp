#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <algorithm>
#include <cmath>
#include <climits>

using namespace std;

mutex mtx;

struct Cell {
    int i, j, type;
};

struct Subproblem_p {
    int i0, j0, i1, j1;
    int startType, endType;
};
// Example scoring and gap penalty functions
int score(char a, char b) {
    return (a == b) ? 2 : -1;
}

int gapPenalty(int k) {
    return 2 + k;
}

// Function to initialize tables
void initializeTables(vector<vector<int>>& T1, vector<vector<int>>& T2, vector<vector<int>>& T3, int m, int n, int startType) {
    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j <= n; ++j) {
            T1[i][j] = T2[i][j] = T3[i][j] = INT_MIN;
        }
    }
    T1[0][0] = T2[0][0] = T3[0][0] = 0;
}

// Function to fill the tables
void fillTables(vector<vector<int>>& T1, vector<vector<int>>& T2, vector<vector<int>>& T3, const string& A, const string& B) {
    int m = A.size();
    int n = B.size();
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            T1[i][j] = score(A[i-1], B[j-1]) + max({T1[i-1][j-1], T2[i-1][j-1], T3[i-1][j-1]});
            T2[i][j] = max({T1[i][j-1] - gapPenalty(1), T2[i][j-1] - 1, T3[i][j-1] - gapPenalty(1)});
            T3[i][j] = max({T1[i-1][j] - gapPenalty(1), T2[i-1][j] - gapPenalty(1), T3[i-1][j] - 1});
        }
    }
}

// Function to find the partial balanced partition
void findPartition(vector<Subproblem_p>& partitions, const string& A, const string& B, int p) {
    int m = A.size();
    int n = B.size();
    vector<vector<int>> T1(m + 1, vector<int>(n + 1));
    vector<vector<int>> T2(m + 1, vector<int>(n + 1));
    vector<vector<int>> T3(m + 1, vector<int>(n + 1));

    initializeTables(T1, T2, T3, m, n, -1);
    fillTables(T1, T2, T3, A, B);

    int rowPartSize = m / p;
    int colPartSize = n / p;

    for (int i = 0; i < p; ++i) {
        int i0 = i * rowPartSize;
        int i1 = (i + 1) * rowPartSize - 1;
        int j0 = i * colPartSize;
        int j1 = (i + 1) * colPartSize - 1;

        int startType = (i == 0) ? -1 : (i % 3) + 1;
        int endType = ((i + 1) % 3) + 1;

        partitions.push_back({i0, j0, i1, j1, startType, endType});
    }
}

// Thread function to process Subproblem_ps
void processSubproblems(vector<Subproblem_p>& partitions, int threadId, int p) {
    int subproblemSize = partitions.size() / p;
    int start = threadId * subproblemSize;
    int end = (threadId + 1) * subproblemSize;

    for (int i = start; i < end; ++i) {
        // Processing the subproblem (for demonstration, just print the details)
        mtx.lock();
        cout << "Thread " << threadId << " processing subproblem: (" << partitions[i].i0 << "," << partitions[i].j0 << ") to ("
             << partitions[i].i1 << "," << partitions[i].j1 << "), startType: " << partitions[i].startType << ", endType: " << partitions[i].endType << endl;
        mtx.unlock();
    }
}
