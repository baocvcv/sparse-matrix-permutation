#ifndef SEARCH_H_
#define SEARCH_H_

#include <memory>

#include "camd.h"
#include "cholmod.h"

// 
//TODO: add "step"
/**
 * Perform A* search using amd as heuristic.
 * 
 * @param w f(x) = g(x) + w * h(x)
 * @param Nnum The no. of children evaluated at each node
 * @param log_level 0-no log, 1-info, 2-debug
 */
int* A_star_amd(cholmod_sparse* A, cholmod_common* cp,
                int w, int Nnum, int log_level=0);

/**
 * Perform A_star_amd with mpi
 * 
 * @param pid Process rank
 * @param nproc #processes
 * @param P_result Out : Stores the elmination sequence
 * @param w f(x) = g(x) + w * h(x)
 * @param Nnum The no. of children evaluated at each node
 * @param step After selecting the best children, the algorithm
 * performs the elmination indicated by the child, along with 
 * #(step-1) elminations suggested by the AMD algorithm.
 */
void mpi_A_star_amd(int pid, int nproc, cholmod_sparse* A, cholmod_common* cp,
                    std::unique_ptr<int[]>& P_result,
                    int w, int Nnum, int step=1, int log_level=0);

/** Frees memory used by A_start_amd */
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

/** Swaps P[idx1] and P[idx2] */
void swap_P(int idx1, int idx2, int* P);

#endif