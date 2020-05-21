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
#include <iostream>

#include "camd.h"
#include "cholmod.h"

#include "args.hxx"
#include "utils.h"
#include "search.h"

int main(int argc, char** argv)
{
    args::ArgumentParser parser("Sparse matrix permutation parameters", "");
    args::ValueFlag<std::string> arg_input(parser, "filename", "Path to input file", {'f', "file"});
    args::ValueFlag<int> arg_child_num(parser, "num_child", "Num of child of each A* node", {'n', "num_child"});
    args::ValueFlag<int> arg_step(parser, "step", "Steps to take at each node", {'s', "step"});
    args::HelpFlag help(parser, "help", "Display help menu", {'h', "help"});

    try {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help) {
        std::cout << parser;
        return 0;
    }
    catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    catch (args::ValidationError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    std::string input_file = "";
    int child_num = 32;
    int step = 1;

    if (arg_input) input_file = args::get(arg_input);
    if (arg_child_num) child_num = args::get(arg_child_num);
    if (arg_step) step = args::get(arg_step);

    std::cout << "Input file: " << input_file << std::endl;
    std::cout << "Child num: " << child_num << std::endl;
    std::cout << "Step: " << step << std::endl;

    FILE* fp;
    fp = fopen(input_file.c_str(), "r");
    if (!fp) {
        puts("Can't find trip matrix!");
        return 1;
    }

    int weight = 1; // default value5
    //TODO: should weight be used??
    // if (argc == 5) {
    //     // weight = atoi (argv [4]) ;
    //     // printf ("user-set weight: %d", weight) ;
    //     candidates = atoi(argv[4]);
    //     printf("user-set candidate number: %d", candidates);
    // }

    cholmod_sparse* A;
    cholmod_common c;
    cholmod_start(&c); /* start CHOLMOD */
    // c.print = 5 ;                       /* set print level to the highest */
    c.supernodal = 2;

    // A = cholmod_read_triplet (fp, &c) ;
    A = cholmod_read_sparse(fp, &c); /* read in a matrix */
    fclose(fp);
    if (A == NULL || A->stype == 0) /* A must be symmetric */
    {
        puts("not symmetric!");
        cholmod_free_sparse(&A, &c);
        cholmod_finish(&c);
        return 0;
    }

    double beg = omp_get_wtime();
    auto P = A_star_amd(A, &c, weight, child_num);
    // auto P = test_parl (A, &c) ;
    double end = omp_get_wtime();
    printf ("Time elapsed is %lf seconds.\n", end - beg) ;
    // printf("%lf ", end - beg);

    // cholmod_print_perm (P, A->nrow, A->nrow, "P_res", &c) ;

    // verification
    verification(A, &c, P);
    output_file(P, A->ncol);

    A_star_free_all(A, &c, P);
    return 0;
}
