#ifndef UTILS_H_
#define UTILS_H_

#include <string>

#include "camd.h"
#include "cholmod.h"

// returns the number of non-zero elements after amd
int nnz_amd(cholmod_sparse* A, cholmod_common* cp, int* P, int* Parent,
    int* Post, int* ColCount, int* First, int* Level);

// returns the number of increased fill-ins after amd
int fill_in_amd(cholmod_sparse* A, cholmod_common* cp, int* P, int* Parent,
    int* Post, int* ColCount, int* First, int* Level, int A_nnz);

void add_tri_entryd(cholmod_triplet* T, int r, int c, double x);

// P: the ordering of columns
void verification(cholmod_sparse* A, cholmod_common* cp, int* P);

// checks if A through P has the same # of reduced nnz as X
bool check(cholmod_sparse* A, cholmod_common* cp, int* P, int x);

int* test_parl(cholmod_sparse* A, cholmod_common* cp);

void output_file(int* P, int n);

#endif