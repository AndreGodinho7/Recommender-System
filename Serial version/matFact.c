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
    B = MatrixInit(init->nU, init->nI);
    pre_L = MatrixInit(init->nU, init->nF); // Matrix that stores the previous iteration of L
    pre_R = MatrixInit(init->nF, init->nI); // Matrix that stores the previous iteration of R
    pre_B = MatrixInit(init->nU, init->nI); // Matrix that stores the previous iteration of B

    A = init->matrix;
    
    random_fill_LR(pre_L, pre_R, init->nU, init->nI, init->nF);

    /* Para ajudar a ver os resultados */
    printf("Matrix iniciais=== \n");
    printf("\n====  pre_L  =====\n");
    printMatrix(pre_L,init->nU,init->nF);
    printf("\n====  pre_R  =====\n");
    printMatrix(pre_R,init->nF,init->nI);

    matrix_mul(pre_L, pre_R, pre_B, init->nU, init->nI, init->nF);
    printf("\n====  pre_B  =====\n");
    printMatrix(pre_B,init->nU,init->nI);

   /*Do all iterations */
    for(int i = 0 ; i < init->iter ; i++){
        /*update the matrix*/
        if(i > 0 ){
            pre_L = L;
            pre_B = B;
        }
        
        recalculate_Matrix(L, R, pre_L, pre_R, A, B, pre_B,init->nU, init->nI, init->nF, init->iter,init->alpha,init->v ,init->non_zeros);
        matrix_mul(L, R, B, init->nU, init->nI, init->nF);

        if(i<5){
            /*Para ajudar a ver os resultados */ 
            printf("Matrix iter = %d",i);
            printf("\n====  L  =====\n");
            printMatrix(L,init->nU,init->nF);
            printf("\n====  R  =====\n");
            printMatrix(R,init->nF,init->nI);
            printf("\n====  B  =====\n");
            printMatrix(B,init->nU,init->nI);
        }

    }

    printf("\n====  B FINAL  =====\n");
    printMatrix(B,init->nU,init->nI);

    return 0;
}