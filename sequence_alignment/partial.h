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

struct Cell {
    int i, j, type;
};

struct Subproblem_p {
    int i0, j0, i1, j1;
    int startType, endType;
};


// Example scoring and gap penalty functions
int score(char a, char b);

int gapPenalty(int k);

// Function to initialize tables
void initializeTables(vector<vector<int>>& T1, vector<vector<int>>& T2, vector<vector<int>>& T3, int m, int n, int startType);

// Function to fill the tables
void fillTables(vector<vector<int>>& T1, vector<vector<int>>& T2, vector<vector<int>>& T3, const string& A, const string& B);


// Function to find the partial balanced partition
void findPartition(vector<Subproblem_p>& partitions, const string& A, const string& B, int p);

// Thread function to process subproblems
void processSubproblems(vector<Subproblem_p>& partitions, int threadId, int p);

#endif