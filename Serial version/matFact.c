#include <stdio.h>
#include <stdlib.h>

#include"input.h"
#include "factorization.h"
#include <time.h>

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
    for(int j=0;j<init->non_zeros;j++){
        printf("row :%d ",init->v[j].row);
        printf("column :%d ",init->v[j].column);
        printf("value :%f \n",init->v[j].value);
    }
    L = MatrixInit(init->nU, init->nF);
    //R = MatrixInit(init->nF, init->nI);
    R = MatrixInit(init->nI, init->nF);
    B = MatrixInit(init->nU, init->nI);
    pre_L = MatrixInit(init->nU, init->nF); // Matrix that stores the previous iteration of L
    pre_R = MatrixInit(init->nF, init->nI); // Matrix that stores the previous iteration of R
    pre_B = MatrixInit(init->nU, init->nI); // Matrix that stores the previous iteration of B

    A = init->matrix;
    
    random_fill_LR(pre_L, pre_R, init->nU, init->nI, init->nF);
    copy_matrix(pre_L,L,init->nU,init->nF);
    copy_matrix(pre_R,R,init->nF,init->nI); 
    pre_R = transpose(pre_R, init->nF, init->nI); //add 

    matrix_mul(pre_L, pre_R, pre_B, init->nU, init->nI, init->nF);
    

   /*Do all iterations */
   clock_t begin = clock();
    for(int i = 0 ; i < init->iter ; i++){
        /*update the matrix*/
        if(i > 0 ){
            copy_matrix(L,pre_L,init->nU,init->nF);
            //copy_matrix(R,pre_R,init->nF,init->nI);
            copy_matrix(R, pre_R, init->nI, init->nF); //add
            copy_matrix(B,pre_B,init->nU,init->nI);
            
        }

        recalculate_Matrix(L, R, pre_L, pre_R, A, B, pre_B,init->nU, init->nI, init->nF, init->iter,init->alpha,init->v ,init->non_zeros);

        matrix_mul(L, R, B, init->nU, init->nI, init->nF);

        

    }
  
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Execution time: %lf seconds\n", time_spent);
    //matrix_mul(L, R, B, init->nU, init->nI, init->nF);
    //printf("LInha 11 \n\n");
    //printMatrix(B,init->nU,init->nI);
    //printf("LInha L \n\n");
    //printMatrix(L,init->nU,init->nF);
    //printf("LInha R \n\n");
    //printMatrix(R,init->nF,init->nI);    
    create_output(B,init->nU, init->nI,argv[1],A);



    return 0;
}