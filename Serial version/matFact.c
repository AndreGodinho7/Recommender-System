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

    recalculate_Matrix(&L,&R,&pre_L,&pre_R,&A,&B,&pre_B,init->nU, init->nI, init->nF,init->iter,init->alpha);
    
    return 0;
}