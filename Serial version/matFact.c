#include <stdio.h>
#include <stdlib.h>

#include"input.h"
#include "factorization.h"

int main(int argc, char* argv[])
{
    input_values* init;
    double** L, **R;
    double **pre_L,**pre_R; 
    double **A;
    double **B,**pre_B;
 
    if (argc != 2){
        printf("ERROR: inserted more than 1 input file.\n");
        exit(0);
    };

    init = read_input(argv[1]);

    L = MatrixInit(init->nU, init->nF);
    R = MatrixInit(init->nF, init->nI);
    pre_L = MatrixInit(init->nU, init->nF);
    pre_R = MatrixInit(init->nF, init->nI);
    B=MatrixInit(init->nU, init->nI);
    pre_B=MatrixInit(init->nU, init->nI);

    A = init->matrix;

    random_fill_LR(&pre_L, &pre_R, init->nU, init->nI, init->nF);
    matrix_mul(&pre_L, &pre_R,&pre_B,init->nU, init->nI, init->nF);

    for (int row=0; row<init->nU; row++)
    {
        for(int columns=0; columns<init->nF; columns++)
            {
                printf("%lf     ", pre_L[row][columns]);
            }
        printf("\n");
    }

    printf("\n");

    for (int row=0; row<init->nF; row++)
    {
        for(int columns=0; columns<init->nI; columns++)
            {
                printf("%lf     ", pre_R[row][columns]);
            }
        printf("\n");
    }
    printf("\n");    
    recalculate_Matrix(&L,&R,&pre_L,&pre_R,&A,&B,&pre_B,init->nU, init->nI, init->nF,init->iter,init->alpha);
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
    printf("\n");

    /*for (int row=0; row<init->nU; row++)
    {
        for(int columns=0; columns<init->nI; columns++)
            {
                printf("%lf     ", pre_B[row][columns]);
            }
        printf("\n");
    }*/

    return 0;
}