#ifndef SUBPROBLEM_ALIGNMENT_H
#define SUBPROBLEM_ALIGNMENT_H

#include <stdio.h>
#include <thread>
#include <vector>

class Subproblem {
public: 
    char *A;
    char *B;
    size_t m; // lenA
    size_t n; // lenB
    size_t id_A;
    size_t id_B; // offsets from the original sequences
    bool invert; // if m > n in the constructor
    size_t p; // num_processors
    int start_type;
    int end_type;
    int g;
    int h;
    std::vector<std::vector<long int>> T1;
    std::vector<std::vector<long int>> T2;
    std::vector<std::vector<long int>> T3;

    Subproblem (char* _A, char* _B, size_t _m, size_t _n, size_t _id_A, size_t _id_B, size_t _p, int start, int end, int _g, int _h) {
        if (_m <= n) {
            A = _A;
            B = _B;
            m = _m;
            n = _n;
            id_A = _id_A;
            id_B = _id_B;
            invert = false;
        }
        else {
            A = _B;
            B = _A;
            m = _n;
            n = _m;
            id_A = _id_B;
            id_B = _id_A;
            invert = true;
        }
        p = _p;
        start_type = start;
        end_type = end;
        g = _g;
        h = _h;

        
        /* initialize the tables */

        //first row

        for (size_t i = 0; i <= m; i++) {
            std::vector<long int> t1(n+1);
            std::vector<long int> t2(n+1);
            std::vector<long int> t3(n+1);
            T1.push_back(t1);
            T2.push_back(t2);
            T3.push_back(t3);
        }
    }

    static void ComputeRowMapThread13(Subproblem *subp, size_t i, size_t start, size_t end);
    static void ComputeOmegaMapThread(Subproblem *subp, size_t i, size_t start, size_t end, std::vector<long int> &omega);
    static void ComputeRowMapThread2(Subproblem *subp, size_t i, size_t start, size_t end, std::vector<long int> &partial);
    static void ComputeFirstRowMapThread(Subproblem *subp, size_t start, size_t end);
    void compute_row(size_t i);
    void compute_tables();
    long int f(size_t i, size_t j) {
        if (A[id_A+i] == B[id_B + j]) {
            return 1;
        }
        return 0;
    }
};

#endif //SUBPROBLEM_ALIGNMENT_H