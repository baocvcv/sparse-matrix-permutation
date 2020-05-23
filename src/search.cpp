#include <omp.h>
#include <mpi.h>
#include <time.h>
#include <unordered_set>
#include <memory>
#include <cstring>
#include <climits>
#include <iostream>
#include <cstdio>

#include "utils.h"
#include "search.h"

int* A_star_amd(cholmod_sparse* A, cholmod_common* cp, int w, int Nnum, int log_level)
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
        if (log_level)
            printf("i: %5d S_nnz: %7d f_score: %7d\n", i, S_nnz, f_cost_array[i]);
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
            auto F = eliminate_node(S, &c, lid, S_nnz, adjacants, adj_size, &diag, 0);
            // cholmod_print_sparse (F, "F0", cp) ;
            // cholmod_print_sparse (F, "F1", cp) ;
            int F_nnz = cholmod_nnz(F, &c);

            // if (log_level == 2)
            //     printf("%d: %d\n", j - i, f_cost_array[j]);
            // F does not contain column[lid]
            // so we need to add the #nnz of column[lid] to f[j]
            // for a fair comparison
            if (diag) {
                f_cost_array[j] += (F_nnz - S_nnz + adj_size + 1);
            } else {
                f_cost_array[j] += (F_nnz - S_nnz + adj_size);
            }
            if (log_level == 2)
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
            if (log_level == 2)
                printf("min_cost: %d, j-cost: %d\n", min_cost,
                (w * f_cost_array[j] + g_cost_array[j]));
            if (min_cost > (w * f_cost_array[j] + g_cost_array[j])) {
                min_idx = j;
                min_cost = w * f_cost_array[j] + g_cost_array[j];
            }
        }
        if (log_level)
            printf("f_score: %d, g_score: %d\n", (int)f_cost_array[min_idx],
                (int)g_cost_array[min_idx]);

        // Apply 1 permutation
        for (int k = min_idx; k > i; --k) {
            swap_P(k - 1, k, P);
        }

        for (int k = i; k < (int)A->nrow; ++k) { //TODO: copy all the way up???
            f_cost_array[k] = f_cost_array[min_idx];
            // f_cost_array[k] = 0 ;
        }

        int* adjacants = new int[S->nrow];
        int* P_tmp = new int[S->nrow];
        int adj_size = 0;
        int lid = min_idx - i;


        auto T = eliminate_node(S, cp, lid, S_nnz, adjacants, adj_size, nullptr, 0);
        cholmod_free_sparse(&S, cp);
        S = T;

        // printf ("min_cost[%d]: %d, min_idx[%d]: %d\n", i, min_cost, i, min_idx -
        // i) ;
        for (int u = 0; u < i + 1; ++u) { // make a copy of global minimum ordering
            P_opt[u] = P[u];
        }
        // perform amd on new S and get order
        (void)camd_order(S->nrow, (int*)S->p, (int*)S->i, P_tmp, nullptr, nullptr,
            nullptr);
        for (int u = i + 1; u < (int)A->ncol; ++u) {
            P_opt[u] = P[P_tmp[u - i - 1] + i + 1];
        }

        // maintain a global min_cost (= min_res)
        if (min_res > min_cost) {
            min_res = min_cost;
        }
        if (!check(A, cp, P_opt, min_cost, log_level)) {
            if (log_level == 1) {
                for (int k = 0; k < adj_size; ++k) {
                    printf("%d ", adjacants[k]);
                }
                puts("");
                printf("lid: %d\n", lid);
                printf("F_nnz: %d\n", cholmod_nnz(S, cp));
                puts("dump last matrix");
            }
            FILE* fp = fopen(
                std::string("out/Matrix_" + std::to_string(i + 1) + ".txt").c_str(),
                "w");
            if (!fp) { puts("fail to write!"); }
            for (int u = 0; u < (int)S->ncol; ++u) {
                for (int v = ((int*)S->p)[u]; v < ((int*)S->p)[u + 1]; ++v) {
                    fprintf(fp, "%d %d 1\n", ((int*)S->i)[v], u);
                }
            }
            fclose(fp);
            return nullptr;
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

bool mpi_all_agree (std::unique_ptr<int[]>& t, int len, int pid)
{
    if (pid == 0) {
        auto t_min = std::make_unique<int[]>(len);
        auto t_max = std::make_unique<int[]>(len);
        MPI_Reduce(t.get(), t_min.get(), len, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(t.get(), t_max.get(), len, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        for (int i = 0; i < len; i++) {
            if (t_min[i] != t_max[i])
                return false;
        }
        return true;
    } else {
        MPI_Reduce(t.get(), nullptr, len, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(t.get(), nullptr, len, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        return true;
    }
}

struct Info {
    int min_cost;
    int col_no; // the column to eliminate
    int f_cost;
};
// custom min function for reduction
void cs_min (void *in, void *inout, int *len, MPI_Datatype *d_type)
{
    //TODO: check if works
    Info* in_ = static_cast<Info*>(in);
    Info* inout_ = static_cast<Info*>(inout);
    for (int i = 0; i < *len; i++) {
        if (in_->min_cost < inout_->min_cost) {
            *inout_ = *in_;
        } else if (in_->min_cost == inout_->min_cost && in_->col_no < inout_->col_no) {
            *inout_ = *in_;
        }
        in_++, inout_++;
    }
}

void mpi_A_star_amd(int pid, int nproc, cholmod_sparse* A, cholmod_common* cp,
                    int w, int Nnum, std::unique_ptr<int[]>& P_result, int log_level)
{
    // custom op and datatype used later in reduction
    MPI_Op myOp;
    MPI_Datatype myType;
    MPI_Type_contiguous(3, MPI_INT, &myType);
    MPI_Type_commit(&myType);
    MPI_Op_create(cs_min, true, &myOp);

    // each process maintains its own cost_array
    // int* P = new int[A->nrow];
    // int* f_cost_array = new int[A->nrow]; // g(x)
    // int* g_cost_array = new int[A->nrow]; // h(x)
    // int* part_j = new int[Nnum];
    // bool* non_min = new bool[Nnum];

    auto P = std::make_unique<int[]>(A->nrow);
    auto f_cost_array = std::make_unique<int[]>(A->nrow);
    auto g_cost_array = std::make_unique<int[]>(A->nrow);
    auto part_j = std::make_unique<int[]>(Nnum);
    auto non_min = std::make_unique<bool[]>(Nnum);

    for (int i = 0; i < (int)A->nrow; ++i) {
        P[i] = i;
        f_cost_array[i] = 0;
        g_cost_array[i] = 0;
    }

    // dimension of S changes every loop
    auto S = cholmod_copy_sparse(A, cp);
    int min_res = INT_MAX;

    auto P_opt = std::make_unique<int[]>(S->nrow);

    for (int i = 0; i < (int)A->nrow - 1; i++) {
        // check if full
        int S_nnz = cholmod_nnz(S, cp);
        if (log_level)
            printf("i: %5d S_nnz: %7d f_score: %7d\n", i, S_nnz, f_cost_array[i]);
        int full = (S->nrow - i) * (S->nrow - i + 1) / 2;
        if (S_nnz == full) {
            break;
        }

        printf("i=%d pid=%d\n", i, pid);
        
        // at each step, process 0 generates nproc * Nnum random numbers in [0, i]
        // scatters the job to #nproc processes
        int thread_limit = nproc * Nnum;
        int max_child_no = (int)A->nrow - i;
        int jobs_len = std::min(thread_limit, max_child_no);
        auto send_cts = std::make_unique<int[]>(nproc); // the number of jobs each process takes
        int stride = jobs_len / nproc;
        int remainder = jobs_len % nproc;
        for (int j = 0; j < nproc; j++) {
            send_cts[j] = stride + (j < remainder ? 1 : 0);
        }
        if (pid == 0) {
#ifndef DEBUG
            srand(time(nullptr));
#endif
            auto jobs = std::make_unique<int[]>(jobs_len);
            if (thread_limit < max_child_no) {
                // generate random numbers
                std::unordered_set<int> nums;
                while (nums.size() < thread_limit) {
                    nums.insert(rand() % (max_child_no) + i);
                }
                int j = 0;
                printf("Jobs:");
                for (auto num: nums) {
                    jobs[j++] = num;
                }
            } else {
                // use [i, A->nrow) directly
                printf("Jobs_low:");
                for (int j = 0; j < max_child_no; j++) {
                    jobs[j] = i + j;
                }
            }
            for (int k = 0; k < jobs_len; k++) {
                printf("%d ", jobs[k]);
            }
            puts("");

            // scatter jobs into part_j[]
            auto displs = std::make_unique<int[]>(nproc);
            displs[0] = 0;
            for (int j = 1; j < nproc; j++) {
                displs[j] = displs[j-1] + send_cts[j-1];
            }
            printf("send_cts(displs): ");
            for (int k = 0; k < nproc; k++) {
                printf("%d(%d) ", send_cts[k], displs[k]);
            }
            printf("\n");

            MPI_Scatterv(jobs.get(), send_cts.get(), displs.get(), MPI_INT,
                         part_j.get(), send_cts[pid], MPI_INT, 0, MPI_COMM_WORLD);

        } else {
            MPI_Scatterv(nullptr, nullptr, nullptr, MPI_INT,
                         part_j.get(), send_cts[pid], MPI_INT, 0, MPI_COMM_WORLD);
        }
        // TODO: test that the jobs are received correctly
        auto print = [&](int now_pid) {
            if (pid == now_pid) {
                std::cout << "pid=" << pid << ": ";
                for (int i = 0; i < send_cts[pid]; i++) {
                    std::cout << part_j[i] << " ";
                }
                std::cout << std::endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        };
        for (int i = 0; i < nproc; i++) {
            print(i);
        }

        // consider each candidate using omp parallel for
        #pragma omp parallel for
        // consider each candidate
        // send_cts[pid] is the #jobs assigned to this process
        for (int jj = 0; jj < send_cts[pid]; ++jj) {
            bool diag = false;
            int j = part_j[jj]; // the actual column

            auto Parent = std::make_unique<int[]>(S->nrow);
            auto Post = std::make_unique<int[]>(S->nrow);
            auto ColCount = std::make_unique<int[]>(S->nrow);
            auto First = std::make_unique<int[]>(S->nrow);
            auto Level = std::make_unique<int[]>(S->nrow);
            auto P_tmp = std::make_unique<int[]>(S->nrow);
            auto adjacants = std::make_unique<int[]>(S->nrow);

            cholmod_common c;
            cholmod_start(&c);

            int adj_size = 0;
            int lid = j - i; // row to use for elimination in current S
            auto F = eliminate_node(S, &c, lid, S_nnz, adjacants.get(), adj_size, &diag, 0);
            // cholmod_print_sparse (F, "F0", cp) ;
            // cholmod_print_sparse (F, "F1", cp) ;
            int F_nnz = cholmod_nnz(F, &c);

            // if (log_level == 2)
            //     printf("%d: %d\n", j - i, f_cost_array[j]);
            // F does not contain column[lid]
            // so we need to add the #nnz of column[lid] to f[j]
            // for a fair comparison
            if (diag) {
                f_cost_array[j] += (F_nnz - S_nnz + adj_size + 1);
            } else {
                f_cost_array[j] += (F_nnz - S_nnz + adj_size);
            }
            if (log_level == 2)
                printf("%d: %d\n", j - i, f_cost_array[j]);
            g_cost_array[j] = fill_in_amd(F, &c, P_tmp.get(), Parent.get(), Post.get(),
                                          ColCount.get(), First.get(), Level.get(), F_nnz);
            // printf ("g_cost_array: %d\n", g_cost_array[j]) ;
            // TODO: can cache F and use outside the loop
            cholmod_free_sparse(&F, &c);
            cholmod_free_work(&c);
        }

        // get the best child
        int min_cost = w * f_cost_array[part_j[0]] + g_cost_array[part_j[0]];
        int min_idx = part_j[0];
        for (int jj = 1; jj < send_cts[pid]; ++jj) {
            int j = part_j[jj];
            if (log_level == 2) {
                printf("min_cost: %d, j-cost: %d\n", min_cost,
                (w * f_cost_array[j] + g_cost_array[j]));
            }
            auto cost = w * f_cost_array[j] + g_cost_array[j];
            if (min_cost > cost) {
                min_idx = j;
                min_cost = cost;
            }
        }
        if (log_level) {
            printf("f_score: %d, g_score: %d\n", (int)f_cost_array[min_idx],
                (int)g_cost_array[min_idx]);
        }
        printf("Finish evaluating pid=%d min_cost=%d min_idx=%d\n", pid, min_cost, min_idx);

        // gather the result and select the best overall child
        // scatter this information and each node should update its matrix
        // TODO: try only update the matrix at that node and scatter it to other nodes
        Info myInfo, bestInfo;
        myInfo.min_cost = min_cost;
        myInfo.col_no = min_idx;
        myInfo.f_cost = f_cost_array[min_idx];
        MPI_Allreduce(&myInfo, &bestInfo, 1, myType, myOp, MPI_COMM_WORLD);
        min_idx = bestInfo.col_no;
        min_cost = bestInfo.min_cost;
        f_cost_array[min_idx] = bestInfo.f_cost;
        printf("Reduction finished pid=%d min_cost=%d min_idx=%d\n", pid, min_cost, min_idx);

        // update matrix using bestInfo
        for (int k = min_idx, tmp; k > i; --k) {
            tmp = P[k], P[k] = P[k-1], P[k-1] = tmp;
        }

        for (int k = i; k < (int)A->nrow; ++k) { //TODO: copy all the way up???
            f_cost_array[k] = f_cost_array[min_idx];
            // f_cost_array[k] = 0 ;
        }
        // TODO: elminate node 
        auto P_tmp = std::make_unique<int[]>(S->nrow);
        auto adjacants = std::make_unique<int[]>(S->nrow);
        int adj_size = 0;
        int lid = min_idx - i;

        auto T = eliminate_node(S, cp, lid, S_nnz, adjacants.get(), adj_size, nullptr, 0);
        cholmod_free_sparse(&S, cp);
        S = T;

        memcpy(P_opt.get(), P.get(), (i+1) * sizeof(int));
        camd_order(S->nrow, (int*)S->p, (int*)S->i, P_tmp.get(), nullptr, nullptr,
            nullptr);
        // record ordering proposed by amd
        for (int u = i + 1; u < (int)A->ncol; ++u) {
            P_opt[u] = P[P_tmp[u - i - 1] + i + 1];
        }

        if (mpi_all_agree(f_cost_array, A->nrow, pid) && pid == 0) printf("f agrees\n");
        else if (pid == 0) printf("f does not agree\n");
        if (mpi_all_agree(P, A->nrow, pid) && pid == 0) printf("P agrees\n");
        else if (pid == 0) printf("P does not agree\n");
        if (mpi_all_agree(P_opt, S->nrow, pid) && pid == 0) printf("P_opt agrees\n");
        else if (pid == 0) printf("P_opt does not agree\n");
        printf("Finish eliminating node pid=%d\n", pid);

        // check for errors
        // TODO: not needed?
        if (!check(A, cp, P_opt.get(), min_cost, log_level)) {
            if (pid == 0) {
                if (log_level == 1) {
                    for (int k = 0; k < adj_size; ++k) {
                        printf("%d ", adjacants[k]);
                    }
                    puts("");
                    printf("lid: %d\n", lid);
                    printf("F_nnz: %d\n", cholmod_nnz(S, cp));
                    puts("dump last matrix");
                }
                FILE* fp = fopen(
                    std::string("out/Matrix_" + std::to_string(i + 1) + ".txt").c_str(),
                    "w");
                if (!fp) { puts("fail to write!"); }
                for (int u = 0; u < (int)S->ncol; ++u) {
                    for (int v = ((int*)S->p)[u]; v < ((int*)S->p)[u + 1]; ++v) {
                        fprintf(fp, "%d %d 1\n", ((int*)S->i)[v], u);
                    }
                }
                fclose(fp);
            }
            P_result.reset();
            return;
        }
        if (min_res > min_cost) { // maintain a global min cost
            min_res = min_cost;
        }
        printf("Finish checking pid=%d\n", pid); MPI_Barrier(MPI_COMM_WORLD);
    }
    // process 0 should set P and all process should clean up and return
    P_result = std::move(P_opt);

    MPI_Type_free(&myType);
    MPI_Op_free(&myOp);
}


void A_star_free_all(cholmod_sparse* A, cholmod_common* cp, int* P)
{
    if (P != nullptr)
        delete[] P;
    cholmod_free_sparse(&A, cp);
    cholmod_finish(cp); /* finish CHOLMOD */

    return;
}

cholmod_sparse* eliminate_node(cholmod_sparse* S, cholmod_common* cp, int lid, int S_nnz,
                    int* adjacants, int& adj_size, bool* diag, int log_level=0)
{
    //int lid = j - i; // row to use for elimination in current S
    // discover adjacants in S
    //TODO: can count row first and then column
    for (int k = 0; k < (int)S->ncol; ++k) { // for each column of S
        int start = ((int*)S->p)[k];
        int finish = ((int*)S->p)[k+1];
        if (k == lid) {
            // Apply change indicated by lid
            // for each element in column lid
            for (int u = start; u < finish; ++u) {
                // upper matrix's row index is smaller than column index
                if (((int*)S->i)[u] == lid) {
                    if (diag != nullptr)
                        *diag = true;
                    continue;
                }
                if (log_level)
                    printf("delete edge (%d, %d)\n", ((int*)S->i)[u], lid);
                // record the adjacent columns
                adjacants[adj_size++] = ((int*)S->i)[u];
            }
            continue;
        }

        for (int u = start; u < finish; ++u) {
            if (((int*)S->i)[u] == lid) {
                if (log_level)
                    printf("delete edge (%d, %d)\n", ((int*)S->i)[u], k);
                // upper matrix's row index is smaller than column index
                adjacants[adj_size++] = k - 1;
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
        int start = ((int*)S->p)[k];
        int finish = ((int*)S->p)[k+1];

        int c = k > lid ? (k - 1) : k;
        for (int u = start; u < finish; ++u) {
            int col = ((int*)S->i)[u];
            if (col != lid) {
                int r = col > lid ? (col - 1) : col;
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

inline void swap_P(int idx1, int idx2, int* P)
{
    int tmp = P[idx1];
    P[idx1] = P[idx2];
    P[idx2] = tmp;

    return;
}