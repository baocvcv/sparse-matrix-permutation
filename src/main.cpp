/* ========================================================================== */
/* === Demo/cholmod_simple ================================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Demo Module.  Copyright (C) 2005-2006, Timothy A. Davis
 * -------------------------------------------------------------------------- */

/* Read in a real symmetric or complex Hermitian matrix from stdin in
 * MatrixMarket format, solve Ax=b where b=[1 1 ... 1]', and print the residual.
 * Usage: cholmod_simple < matrixfile
 */
#include <limits.h>
#include <omp.h>

#include <string>

#include "camd.h"
#include "cholmod.h"

// return the number of fillings
int fill_in_amd(cholmod_sparse* A, cholmod_common* cp, int* P, int* Parent,
    int* Post, int* ColCount, int* First, int* Level, int A_nnz);
// no of non-zeros elements
int nnz_amd(cholmod_sparse* A, cholmod_common* cp, int* P, int* Parent,
    int* Post, int* ColCount, int* First, int* Level);
// perform A* using amd as heuristic
int* A_star_amd(cholmod_sparse* A, cholmod_common* cp, int w, int Nnum);
// release memory
void A_star_free_all(cholmod_sparse* A, cholmod_common* cp, int* P);
// COMPLICATED
void add_tri_entryd(cholmod_triplet* T, int r, int c, double x);
void swap_P(int idx1, int idx2, int* P);
void verification(cholmod_sparse* A, cholmod_common* cp, int* P);
bool check(cholmod_sparse* A, cholmod_common* cp, int* P, int x);
int* test_parl(cholmod_sparse* A, cholmod_common* cp);
void output_file(int* P, int n);

int main(int argc, char** argv)
{
    if (argc < 3) {
        puts("Please provide matrix market format file with -f!");
        return (0);
    }

    // printf ("Case: %s\n", argv[2]) ;   // print 1
    FILE* fp;
    fp = fopen(argv[2], "r");
    if (!fp) {
        puts("can't find trip matrix!");
        return (0);
    }

    int weight = 1; // default value5
    int candidates = 32;
    if (argc == 5) {
        // weight = atoi (argv [4]) ;
        // printf ("user-set weight: %d", weight) ;
        candidates = atoi(argv[4]);
        printf("user-set candidate number: %d", candidates);
    }

    cholmod_sparse* A;
    cholmod_common c;
    cholmod_start(&c); /* start CHOLMOD */
    // c.print = 5 ;                       /* set print level to the highest */
    c.supernodal = 2;

    // A = cholmod_read_triplet (fp, &c) ;
    A = cholmod_read_sparse(fp, &c); /* read in a matrix */
    if (A == NULL || A->stype == 0) /* A must be symmetric */
    {
        puts("not symmetric!");
        cholmod_free_sparse(&A, &c);
        cholmod_finish(&c);
        return (0);
    }

    double beg = omp_get_wtime();
    auto P = A_star_amd(A, &c, weight, candidates);
    // auto P = test_parl (A, &c) ;
    double end = omp_get_wtime();
    // printf ("Time elapsed is %lf seconds.\n", end - beg) ;   // print 2
    printf("%lf ", end - beg);

    // cholmod_print_perm (P, A->nrow, A->nrow, "P_res", &c) ;

    // verification
    verification(A, &c, P);
    output_file(P, A->ncol);

    A_star_free_all(A, &c, P);
    fclose(fp);
    return (0);
}

void output_file(int* P, int n)
{
    FILE* fp = fopen("label.out", "a");

    for (int i = 0; i < n; ++i) {
        fprintf(fp, "%d ", P[i]);
    }
    fprintf(fp, "\n");

    fclose(fp);
    return;
}

// P: the ordering of columns
void verification(cholmod_sparse* A, cholmod_common* cp, int* P)
{
    int A_nnz = cholmod_nnz(A, cp);

    int* P_tmp = new int[A->nrow];
    int* Parent = new int[A->nrow];
    int* Post = new int[A->nrow];
    int* ColCount = new int[A->nrow];
    int* First = new int[A->nrow];
    int* Level = new int[A->nrow];

    cholmod_analyze_ordering(A, 1, P, nullptr, 0, Parent, Post, ColCount, First,
        Level, cp);
    int nnz = 0;
    for (int i = 0; i < (int)A->nrow; ++i) {
        nnz += ColCount[i];
    }
    // auto A_nnz = cholmod_nnz (A, cp) ;

    // printf ("The final fill-in verification: %d\n", (int)(nnz - A_nnz) ) ;

    // printf ("The amd fill-in: %d\n", fill_in_amd (A, cp, P_tmp, Parent, Post,
    // ColCount, First, Level, A_nnz) ) ;

    // printf ("The final nz verification: %d\n", (int)(nnz) ) ;

    // printf ("The amd nz: %d\n", nnz_amd (A, cp, P_tmp, Parent, Post, ColCount,
    // First, Level) ) ;

    cp->nmethods = 1;
    cp->method[0].ordering = CHOLMOD_METIS;
    cholmod_analyze(A, cp);
    // printf ("The metis nz: %lf\n", cp->method[0].lnz) ;
    printf("%d %d %lf %d %d\n", (int)A->ncol, (int)A_nnz, cp->method[0].lnz,
        nnz_amd(A, cp, P_tmp, Parent, Post, ColCount, First, Level),
        (int)(nnz));

    delete[] P_tmp;
    delete[] Parent;
    delete[] Post;
    delete[] ColCount;
    delete[] First;
    delete[] Level;

    return;
}

// checks if A through P has the same # of reduced nnz as X
bool check(cholmod_sparse* A, cholmod_common* cp, int* P, int x)
{
    int A_nnz = cholmod_nnz(A, cp);

    int* P_tmp = new int[A->nrow];
    int* Parent = new int[A->nrow];
    int* Post = new int[A->nrow];
    int* ColCount = new int[A->nrow];
    int* First = new int[A->nrow];
    int* Level = new int[A->nrow];

    cholmod_analyze_ordering(A, 1, P, nullptr, 0, Parent, Post, ColCount, First,
        Level, cp);
    int nnz = 0;
    for (int i = 0; i < (int)A->nrow; ++i) {
        nnz += ColCount[i];
    }
    printf("The ver fill-in: %d\n", nnz - A_nnz);
    printf("The xxx fill-in: %d\n", x);

    delete[] P_tmp;
    delete[] Parent;
    delete[] Post;
    delete[] ColCount;
    delete[] First;
    delete[] Level;

    return x == (nnz - A_nnz);
}

// returns the number of increased fill-ins after amd
int fill_in_amd(cholmod_sparse* A, cholmod_common* cp, int* P, int* Parent,
    int* Post, int* ColCount, int* First, int* Level, int A_nnz)
{
    return (nnz_amd(A, cp, P, Parent, Post, ColCount, First, Level) - A_nnz);
}

// returns the number of non-zero elements after amd
int nnz_amd(cholmod_sparse* A, cholmod_common* cp, int* P, int* Parent,
    int* Post, int* ColCount, int* First, int* Level)
{
    // cholmod_print_sparse (A, "A in fiil_in_amd", cp) ;		    /* print the
    // matrix */

    (void)camd_order(A->nrow, (int*)A->p, (int*)A->i, P, nullptr, nullptr,
        nullptr);
    // cholmod_print_perm (P, A->nrow, A->nrow, "P_amd", cp) ;

    cholmod_analyze_ordering(A, 1, P, nullptr, 0, Parent, Post, ColCount, First,
        Level, cp);

    int nnz = 0;
    for (int i = 0; i < (int)A->nrow; ++i) {
        nnz += ColCount[i];
    }

    return (nnz);
}

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
            int j = part_j[jj];

            for (int k = 0; k < Nnum; ++k) { // check for duplicates
                if (j == part_j[k] && jj > k) {
                    non_min[jj] = true;
                    break;
                }
            }
            if (non_min[jj]) { // not necessary ?
                continue;
            }

            cholmod_common c;
            cholmod_start(&c);

            int* P_tmp = new int[S->nrow];
            int* Parent = new int[S->nrow];
            int* Post = new int[S->nrow];
            int* ColCount = new int[S->nrow];
            int* First = new int[S->nrow];
            int* Level = new int[S->nrow];
            int* adjacants = new int[S->nrow];

            int adj_size = 0;
            int lid = j - i; // row to use for elimination in current S
            // discover adjacants in S
            //TODO: can count row first and then column
            for (int k = 0; k < (int)S->ncol; ++k) { // for each column of S
                if (k == lid) {
                    // for each element in column lid
                    for (int u = ((int*)S->p)[k]; u < ((int*)S->p)[k + 1]; ++u) {
                        // upper matrix's row index is smaller than column index
                        if (((int*)S->i)[u] == lid) {
                            diag = true;
                            continue;
                        }
                        // record the adjacent columns
                        adjacants[adj_size] = ((int*)S->i)[u];
                        ++adj_size;
                    }
                    continue;
                }

                for (int u = ((int*)S->p)[k]; u < ((int*)S->p)[k + 1]; ++u) {
                    if (((int*)S->i)[u] == lid) {
                        // upper matrix's row index is smaller than column index
                        adjacants[adj_size] = k - 1;
                        // adjacants[adj_size] = k ;
                        ++adj_size;
                    }
                }
            }
            // allocates a triplet matrix
            auto T = cholmod_allocate_triplet(
                S->ncol - 1, S->nrow - 1,
                S_nnz + adj_size * adj_size, // max # of nonzeros
                S->stype, S->xtype, &c);
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
                }
            }

            // cholmod_print_triplet (T, "T1", cp) ;

            // convert back to sparse mat
            // TODO: F still includes t
            // TODO: implement a parallel version for this
            auto F = cholmod_triplet_to_sparse(T, T->nnz, &c);
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
            cholmod_free_sparse(&F, &c);

            cholmod_free_triplet(&T, &c);

            cholmod_free_work(&c);

            delete[] P_tmp;
            delete[] Parent;
            delete[] Post;
            delete[] ColCount;
            delete[] First;
            delete[] Level;
            delete[] adjacants;
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

        // !!! you should pick up from here and modify the following code!!!
        int* adjacants = new int[S->nrow];
        int* P_tmp = new int[S->nrow];
        int lid = min_idx - i;
        int adj_size = 0;
        for (int k = 0; k < (int)S->ncol; ++k) {
            // Apply change indicated by min_idx
            if (k == lid) {
                for (int u = ((int*)S->p)[k]; u < ((int*)S->p)[k + 1]; ++u) {
                    // upper matrix's row index is smaller than column index
                    if (((int*)S->i)[u] == lid) {
                        continue;
                    }
                    printf("delete edge (%d, %d)\n", ((int*)S->i)[u], lid);
                    adjacants[adj_size] = ((int*)S->i)[u];
                    ++adj_size;
                }
                continue;
            }

            for (int u = ((int*)S->p)[k]; u < ((int*)S->p)[k + 1]; ++u) {
                if (((int*)S->i)[u] == lid) {
                    // upper matrix's row index is smaller than column index
                    printf("delete edge (%d, %d)\n", ((int*)S->i)[u], k);
                    adjacants[adj_size] = k - 1;
                    ++adj_size;
                }
            }
        }

        // TODO: adj_size * adj_size is too large a boundary
        auto T = cholmod_allocate_triplet(S->ncol - 1, S->nrow - 1,
            S_nnz + adj_size * adj_size, S->stype,
            S->xtype, cp);
        // add S - column[lid] to T
        for (int k = 0; k < (int)S->ncol; ++k) {
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

        // add fill_in entries
        for (int k = 0; k < adj_size - 1; ++k) {
            for (int u = k + 1; u < adj_size; ++u) {
                if (adjacants[u] < adjacants[k]) {
                    int tmp = adjacants[u];
                    adjacants[u] = adjacants[k];
                    adjacants[k] = tmp;
                }
                add_tri_entryd(T, adjacants[k], adjacants[u], 1);
                int r = adjacants[k] >= lid ? adjacants[k] + 1 : adjacants[k];
                int c = adjacants[u] >= lid ? adjacants[u] + 1 : adjacants[u];
                printf("add edge (%d, %d)\n", r, c);
            }
        }

        // cholmod_print_triplet (T, "T3", cp) ;

        cholmod_free_sparse(&S, cp);
        // TODO: check if has a parallel version
        S = cholmod_triplet_to_sparse(T, T->nnz, cp);
        cholmod_free_triplet(&T, cp);

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

void add_tri_entryd(cholmod_triplet* T, int r, int c, double x)
{
    if (T->nzmax <= T->nnz) {
        puts("No memory of cholmod_tripple for more entry");
        return;
    }

    ((int*)T->i)[T->nnz] = r;
    ((int*)T->j)[T->nnz] = c;
    ((double*)T->x)[T->nnz] = x;

    ++T->nnz;

    return;
}

void swap_P(int idx1, int idx2, int* P)
{
    int tmp = P[idx1];
    P[idx1] = P[idx2];
    P[idx2] = tmp;

    return;
}

int* test_parl(cholmod_sparse* A, cholmod_common* cp)
{
    int* Perm = new int[A->nrow];
    for (int i = 0; i < (int)A->nrow; ++i) {
        Perm[i] = i;
    }

    auto S = cholmod_copy_sparse(A, cp);
    int N = A->nrow * A->ncol;

    cholmod_sparse* F;

    #pragma omp parallel for private(F)
    for (int i = 0; i < N; ++i) {
        int* P = new int[A->ncol];
        int* Parent = new int[A->ncol];
        int* Post = new int[A->ncol];
        int* ColCount = new int[A->ncol];
        int* First = new int[A->ncol];
        int* Level = new int[A->ncol];

        cholmod_common c;
        cholmod_start(&c);

        auto T = cholmod_allocate_triplet(S->ncol, S->nrow, cholmod_nnz(S, &c),
            S->stype, S->xtype, &c);

        for (int k = 0; k < (int)S->ncol; ++k) {
            for (int u = ((int*)S->p)[k]; u < ((int*)S->p)[k + 1]; ++u) {
                add_tri_entryd(T, ((int*)S->i)[u], k, ((double*)S->x)[u]);
            }
        }

        F = cholmod_triplet_to_sparse(T, T->nnz, &c);

        //(void) camd_order (F->nrow, (int*)F->p, (int*)F->i, P, nullptr, nullptr,
        // nullptr) ; cholmod_analyze_ordering (F, 1, P, nullptr, 0, Parent, Post,
        // ColCount, First, Level, &c) ;

        cholmod_free_sparse(&F, &c);
        cholmod_free_triplet(&T, &c);
        cholmod_free_work(&c);

        delete[] P;
        delete[] Parent;
        delete[] Post;
        delete[] ColCount;
        delete[] First;
        delete[] Level;
    }

    cholmod_free_sparse(&S, cp);

    return Perm;
}
