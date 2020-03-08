#include <stdio.h>
#include <stdlib.h>

#include"input.h"
#include "factorization.h"

int main(int argc, char* argv[])
{
    input_values* init;
    double** L, **R; 
 
    if (argc != 2){
        printf("ERROR: inserted more than 1 input file.\n");
        exit(0);
    };

    init = read_input(argv[1]);

    L = MatrixInit(init->nU, init->nF);
    R = MatrixInit(init->nF, init->nI);

    random_fill_LR(&L, &R, init->nU, init->nI, init->nF);

    // só para verificar se está correto 
    // TODO: apagar ...
    for (int row=0; row<init->nU; row++)
    {
        for(int columns=0; columns<init->nF; columns++)
            {
                printf("%lf     ", L[row][columns]);
            }
        printf("\n");
    }

    printf("\n");

    for (int row=0; row<init->nF; row++)
    {
        for(int columns=0; columns<init->nI; columns++)
            {
                printf("%lf     ", R[row][columns]);
            }
        printf("\n");
    }

    return 0;
}