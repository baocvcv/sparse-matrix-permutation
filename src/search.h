#ifndef SEARCH_H_
#define SEARCH_H_

#include <memory>

#include "camd.h"
#include "cholmod.h"

// perform A* using amd as heuristic
//TODO: add "step"
int* A_star_amd(cholmod_sparse* A, cholmod_common* cp, int w, int Nnum, int log_level=0);

// perform A_start_amd with mpi
void mpi_A_star_amd(int pid, int nproc, cholmod_sparse* A, cholmod_common* cp,
                    int w, int Nnum, std::unique_ptr<int[]>& P_result, int log_level=0);

// release memory
void A_star_free_all(cholmod_sparse* A, cholmod_common* cp, int* P);

/*
    Eliminate a node from the current matrix
    @param S: in, the original matrix
    @param cp: in, cholmod environment
    @param col: in, which column of S to eliminate
    @param S_nnz: in, #nnz in S
    @param adjacants: out, adjacant nodes to node[col]
    @param adj_size: out, length of adjacants[]
    @param diag: out, whether S[col][col] is non-zero
    @param log_level: in, 0-no log, 1-info
*/
cholmod_sparse* eliminate_node(cholmod_sparse* S, cholmod_common* cp, int col, int S_nnz,
                    int* adjacants, int& adj_size, bool* diag, int log_level);

void swap_P(int idx1, int idx2, int* P);

#endif