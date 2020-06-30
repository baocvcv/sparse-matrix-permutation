#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char** argv) {
    srand (time (NULL)) ;
    FILE* fp = fopen (argv[1], "w") ;
    if (!fp) {
        puts ("fail to open the assigned file!") ;
    }

    int row = atoi (argv[2]) ;
    int col = atoi (argv[3]) ;
    int nnz = atoi (argv[4]) ;

    fprintf (fp, "%%MatrixMarket matrix coord real symmetric\n") ;
    fprintf (fp, "%d %d %d\n", row, col, nnz) ;

    for (int i = 0; i < row ; ++ i) {
        fprintf (fp, "%d %d %d\n", i, i, 1000) ;
    }

    for (int i = 0; i < nnz - row; ++ i) {
        int c = (rand () % row ) ;
        int r = (rand () % row ) ;

        if (r > c) {
            int tmp = r ;
            r = c ;
            c = tmp ;
        } else if (r == c) {
            c ++ ;
        }

        fprintf (fp, "%d %d %d\n", r, c, 1) ;
    }


    fclose (fp) ;

    return 0 ;
}
