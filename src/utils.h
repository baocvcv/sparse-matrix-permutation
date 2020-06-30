#ifndef UTILS_H_
#define UTILS_H_

#include <string>
#include <mpi.h>

#include "camd.h"
#include "cholmod.h"

/**
 * Returns the number of non-zero elements after amd is performed on A.
 * 
 * @param P Out. The elmination sequence.
 */
int nnz_amd(cholmod_sparse* A, cholmod_common* cp, int* P, int* Parent,
    int* Post, int* ColCount, int* First, int* Level);

/**
 * Returns the number of increased fill-ins after amd is performed on A.
 * 
 * @param P Out. The elmination sequence.
 */
int fill_in_amd(cholmod_sparse* A, cholmod_common* cp, int* P, int* Parent,
    int* Post, int* ColCount, int* First, int* Level, int A_nnz);

/** Add entry x to T[r][c]. */
void add_tri_entryd(cholmod_triplet* T, int row, int col, double val);

/**
 * Compare AMD and custom sequence P.
 * 
 * @param P In. The sequence of elmination.
 */
void verification(cholmod_sparse* A, cholmod_common* cp, int* P);

/**
 * A simpler checker.
 * 
 * Checks if A permutated using sequence P has #(reduced nnz) as x.
 *
 * @param P In. The squence of elimination.
 * @param x In. The corrrect #(reduced nnz).
 */
bool check(cholmod_sparse* A, cholmod_common* cp, int* P, int x, int log_level=0);

/** Test. */
int* test_parl(cholmod_sparse* A, cholmod_common* cp);

/** Output elimination sequence to file. */
void output_file(int* P, int len);

/**
 * Timer for MPI program.
 * 
 * When first called, record the current time and return 0;
 * when called again, return the max difference in time among all processes.
*/
double time_toggle(MPI_Comm comm);

#endif