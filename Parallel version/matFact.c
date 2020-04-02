#include <stdio.h>
#include <stdlib.h>

#include"input.h"
#include "factorization.h"
#include <time.h>
#include <limits.h>
#include<omp.h>

int main(int argc, char* argv[])
{   

    double begin = omp_get_wtime(); 
    input_values* init;
    double sum;
    
    // for allocating matrices
    double** L, **R;
    double **L_hold,**R_hold; 

    // auxiliar pointers to avoid copying matrices
    double** L1, **L2;
    double** R1, **R2; 
    double** tmp;

    // double **A;
    // double **B;
    int i,j,k;

    if (argc != 2){
        printf("ERROR: inserted more than 1 input file.\n");
        exit(0);
    };

    init = read_input(argv[1]);
    // A = init->matrix;
    // B = MatrixInit(init->nU, init->nI);
    L = MatrixInit(init->nU, init->nF);
    R = MatrixInit(init->nF, init->nI);
    L_hold = MatrixInit(init->nU, init->nF); // Matrix that stores the previous iteration of L
    R_hold = MatrixInit(init->nF, init->nI); // Matrix that stores the previous iteration of R
    
    #pragma omp parallel
    {
        #pragma omp single
        {
            random_fill_LR(L, R, init->nU, init->nI, init->nF);
        }
        
        #pragma omp single
        {
            L1 = L;
            L2 = L_hold;
        }

        R = transpose(R, init->nF, init->nI); 
        R_hold = transpose(R_hold, init->nF, init->nI); 

        #pragma omp single
        {
            R1 = R;
            R2 = R_hold;
        }

        matrix_mul(L1, R1, init->v, init->num_zeros, init->nF);

        for(int i = 0 ; i < init->iter ; i++){
            /*update the matrix*/
            #pragma omp single
            {
                tmp = L1;
                L1 = L2;
                L2 = tmp;

                tmp = R1;
                R1 = R2;
                R2 = tmp; 
            }
            recalculate_Matrix(L1, R1, L2, R2, init->nU, init->nI, init->nF, init->alpha,init->v ,init->num_zeros);
            matrix_mul(L1, R1, init->v, init->num_zeros, init->nF);
        }
        
    }
    create_output(init->v, init->nU, init->nI, init->nF, argv[1], L1, R1, init->num_zeros);
    
    for (int i = 0; i < init->nU; i++)
    {
        free(L1[i]);
        free(L2[i]);
    }
    free(L1);
    free(L2);

    for (int i = 0; i < init->nI; i++)
    {
        free(R1[i]);
        free(R2[i]);
    }

    free(R1);
    free(R2);

    free(init);

    double end = omp_get_wtime();
    double time_spent = end - begin;
    printf("Execution time: %lf seconds\n", time_spent);
    return 0;
}