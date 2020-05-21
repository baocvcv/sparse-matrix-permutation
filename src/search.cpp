
#include "utils.h"
#include "search.h"

int* A_star_amd(cholmod_sparse* A, cholmod_common* cp, int w, int Nnum)
{
    int* P = new int[A->nrow];
    int* f_cost_array = new int[A->nrow]; // g(x)
    int* g_cost_array = new int[A->nrow]; // h(x)
    // int Nnum = 32 ;
    int* part_j = new int[Nnum];
    bool* non_min = new bool[Nnum];

    for (int i = 0; i < (int)A->nrow; ++i) {
        P[i] = i;
        f_cost_array[i] = 0;
        g_cost_array[i] = 0;
    }

    // dimension of S changes every loop
    auto S = cholmod_copy_sparse(A, cp);
    int min_res = INT_MAX;

    int* P_opt = new int[S->nrow];

    // TODO: apply multiple amd recommendations
    for (int i = 0; i < (int)A->nrow - 1; ++i) { // run for nrow-1 times!
        printf("i=%d\n", i);
        /*FILE* fp = fopen (std::string("out/Matrix_" + std::to_string(i) +
        ".txt").c_str(), "w") ; if (!fp) { puts ("fail to write!") ;
        }
        for (int u = 0; u < (int)S->ncol; ++ u) {
            for (int v = ((int*)S->p)[u]; v < ((int*)S->p)[u + 1]; ++ v) {
                fprintf (fp, "%d %d 1\n", ((int*)S->i)[v], u) ;
            }
        }
        fclose (fp) ;*/
        int S_nnz = cholmod_nnz(S, cp);
        printf("S_nnz: %d\nf_score: %d\n", S_nnz, f_cost_array[i]);
        int full = (S->nrow - i) * (S->nrow - i + 1) / 2;

        if (S_nnz == full) {
            break;
        }

        #pragma omp parallel for
        for (int j = 0; j < Nnum; ++j) { // generate Nnum random column nums
            // part_j is mapped to column idx in original mat
            part_j[j] = rand() % ((int)A->nrow - i) + i;
            non_min[j] = false;
        }

        //#pragma omp parallel num_threads(num_thread)
        //#pragma omp for
        #pragma omp parallel for
        // for (int j = i; j < (int)A->nrow; ++ j) {
        for (int jj = 0; jj < Nnum; ++jj) { // consider each candidate
            bool diag = false;
            int j = part_j[jj]; // the actual column

            for (int k = 0; k < Nnum; ++k) { // check for duplicates
                if (j == part_j[k] && jj > k) {
                    non_min[jj] = true;
                    break;
                }
            }
            if (non_min[jj]) { // not necessary ?
                continue;
            }

            int* Parent = new int[S->nrow];
            int* Post = new int[S->nrow];
            int* ColCount = new int[S->nrow];
            int* First = new int[S->nrow];
            int* Level = new int[S->nrow];
            int* P_tmp = new int[S->nrow];
            int* adjacants = new int[S->nrow];

            cholmod_common c;
            cholmod_start(&c);

            int adj_size = 0;
            int lid = j - i; // row to use for elimination in current S
            auto F = eliminate_node(S, &c, lid, S_nnz,  adjacants, adj_size, &diag, 0);
            // cholmod_print_sparse (F, "F0", cp) ;
            // cholmod_print_sparse (F, "F1", cp) ;
            int F_nnz = cholmod_nnz(F, &c);

            printf("%d: %d\n", j - i, f_cost_array[j]);
            // F does not contain column[lid]
            // so we need to add the #nnz of column[lid] to f[j]
            // for a fair comparison
            if (diag) {
                f_cost_array[j] += (F_nnz - S_nnz + adj_size + 1);
            } else {
                f_cost_array[j] += (F_nnz - S_nnz + adj_size);
            }
            printf("%d: %d\n", j - i, f_cost_array[j]);
            g_cost_array[j] = fill_in_amd(F, &c, P_tmp, Parent, Post, ColCount, First,
                Level, F_nnz);
            // printf ("g_cost_array: %d\n", g_cost_array[j]) ;
            // TODO: can cache F and use outside the loop
            cholmod_free_sparse(&F, &c);

            cholmod_free_work(&c);

            delete[] P_tmp;
            delete[] adjacants;
            delete[] Parent;
            delete[] Post;
            delete[] ColCount;
            delete[] First;
            delete[] Level;
        }

        // find the min h
        int min_cost = w * f_cost_array[part_j[0]] + g_cost_array[part_j[0]];
        int min_idx = part_j[0];
        // for (int j = i + 1; j < (int)A->nrow; ++ j) {
        for (int jj = 1; jj < Nnum; ++jj) {
            if (non_min[jj]) {
                continue;
            }
            int j = part_j[jj];
            printf("min_cost: %d, j-cost: %d\n", min_cost,
                (w * f_cost_array[j] + g_cost_array[j]));
            if (min_cost > (w * f_cost_array[j] + g_cost_array[j])) {
                min_idx = j;
                min_cost = w * f_cost_array[j] + g_cost_array[j];
            }
        }
        printf("f_score: %d, g_score: %d\n", (int)f_cost_array[min_idx],
            (int)g_cost_array[min_idx]);

        for (int k = min_idx; k > i; --k) {
            swap_P(k - 1, k, P);
        }

        for (int k = i; k < (int)A->nrow; ++k) { //TODO: copy all the way up???
            f_cost_array[k] = f_cost_array[min_idx];
            // f_cost_array[k] = 0 ;
        }

        // Apply 1 permutation
        int* adjacants = new int[S->nrow];
        int* P_tmp = new int[S->nrow];
        int adj_size = 0;
        int lid = min_idx - i;


        auto T = eliminate_node(S, cp, lid, S_nnz, adjacants, adj_size, nullptr, 0);
        cholmod_free_sparse(&S, cp);
        S = T;

        // printf ("min_cost[%d]: %d, min_idx[%d]: %d\n", i, min_cost, i, min_idx -
        // i) ;
        for (int u = 0; u < i + 1; ++u) {
            P_opt[u] = P[u];
        }
        // perform amd on new S
        (void)camd_order(S->nrow, (int*)S->p, (int*)S->i, P_tmp, nullptr, nullptr,
            nullptr);
        for (int u = i + 1; u < (int)A->ncol; ++u) {
            P_opt[u] = P[P_tmp[u - i - 1] + i + 1];
        }
        // maintain a global min_cost (= min_res)
        if (min_res > min_cost) {
            min_res = min_cost;
        }
        if (!check(A, cp, P_opt, min_cost)) {
            for (int k = 0; k < adj_size; ++k) {
                printf("%d ", adjacants[k]);
            }
            puts("");
            printf("lid: %d\n", lid);
            printf("F_nnz: %d\n", cholmod_nnz(S, cp));
            puts("dump last matrix");
            FILE* fp = fopen(
                std::string("out/Matrix_" + std::to_string(i + 1) + ".txt").c_str(),
                "w");
            if (!fp) {
                puts("fail to write!");
            }
            for (int u = 0; u < (int)S->ncol; ++u) {
                for (int v = ((int*)S->p)[u]; v < ((int*)S->p)[u + 1]; ++v) {
                    fprintf(fp, "%d %d 1\n", ((int*)S->i)[v], u);
                }
            }
            fclose(fp);
            exit(1);
        }

        delete[] adjacants;
        delete[] P_tmp;
    }

    delete[] P;
    P = P_opt;

    delete[] f_cost_array;
    delete[] g_cost_array;
    delete[] part_j;
    delete[] non_min;

    return P;
}

void A_star_free_all(cholmod_sparse* A, cholmod_common* cp, int* P)
{
    delete[] P;
    cholmod_free_sparse(&A, cp);
    cholmod_finish(cp); /* finish CHOLMOD */

    return;
}

cholmod_sparse* eliminate_node(cholmod_sparse* S, cholmod_common* cp, int lid, int S_nnz,
                    int* adjacants, int& adj_size, bool* diag, int log_level)
{
    //int lid = j - i; // row to use for elimination in current S
    // discover adjacants in S
    //TODO: can count row first and then column
    for (int k = 0; k < (int)S->ncol; ++k) { // for each column of S
        if (k == lid) {
            // Apply change indicated by lid
            // for each element in column lid
            for (int u = ((int*)S->p)[k]; u < ((int*)S->p)[k + 1]; ++u) {
                // upper matrix's row index is smaller than column index
                if (((int*)S->i)[u] == lid) {
                    if (diag != nullptr)
                        *diag = true;
                    continue;
                }
                if (log_level)
                    printf("delete edge (%d, %d)\n", ((int*)S->i)[u], lid);
                // record the adjacent columns
                adjacants[adj_size] = ((int*)S->i)[u];
                ++adj_size;
            }
            continue;
        }

        for (int u = ((int*)S->p)[k]; u < ((int*)S->p)[k + 1]; ++u) {
            if (((int*)S->i)[u] == lid) {
                if (log_level)
                    printf("delete edge (%d, %d)\n", ((int*)S->i)[u], k);
                // upper matrix's row index is smaller than column index
                adjacants[adj_size] = k - 1;
                // adjacants[adj_size] = k ;
                ++adj_size;
            }
        }
    }
    // allocates a triplet matrix
    // // TODO: adj_size * adj_size is too large a boundary
    auto T = cholmod_allocate_triplet(
        S->ncol - 1, S->nrow - 1,
        S_nnz + adj_size * adj_size, // max # of nonzeros
        S->stype, S->xtype, cp);
    // cholmod_print_sparse (S, "S0", cp) ;
    // cholmod_print_triplet (T, "T0", cp) ;
    // add S - column[lid] to T
    for (int k = 0; k < (int)S->ncol; ++k) {
        // Gaussian elimination
        if (k == lid) {
            continue;
        }

        int c = k > lid ? (k - 1) : k;
        for (int u = ((int*)S->p)[k]; u < ((int*)S->p)[k + 1]; ++u) {
            if (((int*)S->i)[u] != lid) {
                int r = ((int*)S->i)[u] > lid ? (((int*)S->i)[u] - 1)
                                                : ((int*)S->i)[u];
                add_tri_entryd(T, r, c, ((double*)S->x)[u]);
            }
        }
    }
    if (log_level == 2)
        cholmod_print_triplet (T, "T3", cp) ;

    // add fill_in entries
    // TODO: qsort adjacants first??
    for (int k = 0; k < adj_size - 1; ++k) {
        for (int u = k + 1; u < adj_size; ++u) {
            if (adjacants[u] < adjacants[k]) {
                int tmp = adjacants[u];
                adjacants[u] = adjacants[k];
                adjacants[k] = tmp;
            }
            add_tri_entryd(T, adjacants[k], adjacants[u], 1);
            if (log_level) {
                int r = adjacants[k] >= lid ? adjacants[k] + 1 : adjacants[k];
                int c = adjacants[u] >= lid ? adjacants[u] + 1 : adjacants[u];
                printf("add edge (%d, %d)\n", r, c);
            }
        }
    }

    // F still includes t
    // TODO: check if has a parallel version
    auto F = cholmod_triplet_to_sparse(T, T->nnz, cp);
    cholmod_free_triplet(&T, cp);
    return F;
}

void swap_P(int idx1, int idx2, int* P)
{
    int tmp = P[idx1];
    P[idx1] = P[idx2];
    P[idx2] = tmp;

    return;
}