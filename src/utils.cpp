#include <string>

#include "utils.h"

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

int fill_in_amd(cholmod_sparse* A, cholmod_common* cp, int* P, int* Parent,
    int* Post, int* ColCount, int* First, int* Level, int A_nnz)
{
    return (nnz_amd(A, cp, P, Parent, Post, ColCount, First, Level) - A_nnz);
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