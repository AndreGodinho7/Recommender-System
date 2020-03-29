#include <stdio.h>
#include <stdlib.h>
#include<omp.h>
#include "input.h"
#define RAND01 ((double) random() / (double) RAND_MAX)

/******************************************************************************
* printMatrix()
*
* Arguments: matrix - double pointer (matrix)
*            rows - number of rows of matrix
*            columns - number of columns of matrix
*
* Returns: void
*										
* Side-Effects: 
*
* Description: prints matrix
*
*****************************************************************************/
void printMatrix(double** matrix, int rows, int columns)
{
    for (int i = 0; i < rows; i++)
    {
        for(int j = 0; j < columns; j++)
                
            printf("%f     ", matrix[i][j]);
                

        printf("\n");
    }
    printf("\n");
}

/******************************************************************************
* MatrixInit()
*
* Arguments: rows - number of rows of matrix
*            columns - number of columns of matrix
*            
* Returns: pointer to pointer that points to the matrix
*										
* Side-Effects: if any allocation fails exits with exit(1)
*
* Description: alocates memory for a matrix with size [rows x columns]
*
*****************************************************************************/

double** MatrixInit(int rows, int columns)
{   
    double** matrix; 

    matrix = (double **)calloc(rows, sizeof(double *)); 
    if (matrix == NULL){
      printf("ERROR: Out of memory.\n");
		  exit(1);
    }
    for (int i = 0 ; i < rows ; i++) {
         matrix[i] = (double *)calloc(columns, sizeof(double));
    }

    return matrix;
}

/******************************************************************************
* drand()
*
* Arguments: low - lower boundary
*            high - upper boundary         
*
* Returns: double
*										
* Side-Effects: 
*
* Description: generates a random number between the two boundaries
*
*****************************************************************************/

double drand ( double low, double high )
{
    return ( (double)rand() * ( high - low ) ) / (double)RAND_MAX + low;
}

/******************************************************************************
* random_fill_LR()
*
* Arguments: L - passing by reference pointer to pointer of matrix L
*            R - passing by reference pointer to pointer of matrix R
*            nU - number of users          
*            nI - number of items
*            nF - number of features
*
* Returns: void
*										
* Side-Effects: 
*
* Description: random initialization of matrices L and R
*
*****************************************************************************/
void random_fill_LR(double** L_, double** R_, int nU, int nI, int nF)
{   
    srandom(0);
    for(int i = 0; i < nU; i++)
        for(int j = 0; j < nF; j++)
            L_[i][j] = RAND01 / (double) nF;

    for(int i = 0; i < nF; i++)
        for(int j = 0; j < nI; j++)
            R_[i][j] = RAND01 / (double) nF;
}

/******************************************************************************
* transpose()
*
* Arguments: matrix - pointer to pointer of matrix to be tranposed
*            rows - number of rows of matrix
*            columns - number of columns of matrix
*            
* Returns: pointer to pointer of matrix
*										
* Side-Effects: if any allocation fails exits with exit(1)
*               memory of input matrix is freed
*
* Description: allocates memory for a matrix [columns x rows] and copies matrix[i][j]
*              values to the allocated matrix result[j]
*****************************************************************************/

double** transpose(double** matrix, int rows, int columns)
{   
    double** result = MatrixInit(columns, rows);
    int i,j;
    #pragma omp parallel for private(j)
    for (i = 0; i < rows; ++i)
        for (j = 0; j < columns; ++j)
            result[j][i] = matrix[i][j];
    return result;
}


/******************************************************************************
* matrix_mul()
*
* Arguments: firstMatrix - passing by reference pointer to pointer of matrix 
*            secondMatrix - passing by reference pointer to pointer of matrix 
*            matrix3 - passing by reference pointer to pointer of matrix 
*            nU - number of users          
*            nI - number of items
*            nF - number of features
*
* Returns: void
*										
* Side-Effects: 
*
* Description: multiple two matrix and give the result to the pointer matrix3
*
*****************************************************************************/

void matrix_mul(double **firstMatrix, double **secondMatrix, non_zero* v, int num_zeros ,int nF){
    // result[n1][n3] = firstMatrix[n1][n2] x secondMatrix[n3][n2] (transposed)
    //int n1 = nU; // sizeof(matrix_1)/sizeof(matrix_1[0]);
    //int n2 = nF; // sizeof(matrix_1[0])/sizeof(matrix_1[0][0]);
    //int n3 = nI; //sizeof(matrix_2[0])/sizeof(matrix_2[0][0]);
    //double tmp;
    
/*    for(int i = 0; i < n1; i++)
		for(int j = 0; j < n3; j++)
            {
                tmp = 0;
			    for(int k = 0; k < n2; k++)    
                    tmp += firstMatrix[i][k] * secondMatrix[k][j];         
				    result[i][j] = tmp;
            } */
    /*
    for(int i = 0; i < n1; i++)
		for(int j = 0; j < n3; j++)
            {
                tmp = 0;
			    for(int k = 0; k < n2; k++)
                    tmp += firstMatrix[i][k] * secondMatrix[j][k];       
				    result[i][j] = tmp;
            }*/
    
    for (int z = 0; z < num_zeros; z++){
        int i = v[z].row;
        int j = v[z].column;
        v[z].B = 0;
        for (int k = 0; k < nF; k++)
            v[z].B += firstMatrix[i][k]*secondMatrix[j][k];
    }

}

void zero_LR(double** L, double** R, int nU, int nI, int nF){
    int u,i,f;
    //#pragma omp parallel for private(f)
    for(u = 0; u < nU; u++)
        for(f = 0; f < nF; f++)
            L[u][f] = 0;


    //#pragma omp parallel for private(f)
    for(i = 0; i < nI; i++)
        for(f = 0; f < nF; f++)
            R[i][f] = 0;

}


/******************************************************************************
* recalculate_Matrix()
*
* Arguments: L - passing by reference pointer to pointer of matrix L
*            R - passing by reference pointer to pointer of matrix R
*            pre_L - passing by reference pointer to pointer of matrix pre_L
*            pre_R - passing by reference pointer to pointer of matrix pre_R
*            A - passing by reference pointer to pointer of matrix A
*            B - passing by reference pointer to pointer of matrix B
*            pre_B - passing by reference pointer to pointer of matrix pre_B
*            nU - number of users          
*            nI - number of items
*            nF - number of features
*            iter - maximum number of iterations
*            alpha - converge rate
*            v - array of type _non_zero with all information about non zero values in matrix A
*            non_zer0 - number of non zero values in matrix A
*
* Returns: void
*										
* Side-Effects: 
*
* Description: computes the algorithm( minimizing the difference between A and B)
*
*****************************************************************************/

void recalculate_Matrix(double** L, double** R, double** pre_L, double** pre_R, int nU, int nI, int nF, double alpha, non_zero *v, int num_zeros){
    int i, j, z;
    double a, b;
/*
    for(int num_zero = 0 ; num_zero < non_zero ;num_zero++){   
        //printf(" linha: %d  juna : %d\n",v[k].row,v[k].column);
        for(int k = 0; k < nF; k ++){
            double sum_L = 0;
            double sum_R = 0;
            for(int j = 0; j < nI; j++){    
                if(A[v[num_zero].row][j] != 0)
                    //sum_L += (2*((A[v[k].row][j]-pre_B[v[k].row][j])*(-pre_R[feature][j])));
                    //sum_L += (20*((A[v[k].row][j]-internal_product(pre_L[v[k].row],pre_R[j],nF))*(-pre_R[j][k])));
                    // alterar o *2 para o final
                    sum_L += (2*((A[v[num_zero].row][j]-pre_B[v[num_zero].row][j])*(-pre_R[j][k]))); 
                    
            }                                  
            for(int i=0; i < nU; i++){
                if(A[i][v[num_zero].column] != 0)
                    //sum_R += (2*((A[i][v[k].column]-internal_product(pre_L[i],pre_R[v[k].column],nF))*(-pre_L[i][k])));
                    sum_R += (2*((A[i][v[num_zero].column]-pre_B[i][v[num_zero].column])*(-pre_L[i][k])));

            }
            L[v[num_zero].row][k] = pre_L[v[num_zero].row][k] - alpha*sum_L;
            //R[feature][v[k].column] = pre_R[feature][v[k].column] - alpha*sum_R;
            R[v[num_zero].column][k] = pre_R[v[num_zero].column][k] - alpha*sum_R; //add
        }
    } */

    zero_LR(L, R, nU, nI, nF);

    
    for (z = 0; z < num_zeros; z++){
        i = v[z].row;
        j = v[z].column;
        a = v[z].A;
        b = v[z].B;
        //REDUCTION HERE MAYBE
        for(int k = 0; k < nF; k++){
            L[i][k] += (a-b)*(pre_R[j][k]);
            R[j][k] += (a-b)*(pre_L[i][k]);
        }
    }
    
    for(int u = 0; u < nU; u++)
        for(int f = 0; f < nF; f++)
            L[u][f] = pre_L[u][f] + alpha*2*L[u][f];

    
    for(int i = 0; i < nI; i++)
        for(int f = 0; f < nF; f++)
            R[i][f] = pre_R[i][f] + alpha*2*R[i][f]; 

}

void copy_matrix(double** original, double** copied,int rows, int columns){
    for(int i= 0; i<rows;i++){
        for(int j=0;j<columns;j++){
            copied[i][j]=original[i][j];
        }
    }
}

void create_output(double** B,int rows, int columns,char* filename,double** A){
    FILE* fp = fopen("recsystem.out", "w");
    //int size = strlen(filename);
    //printf("%d",size); 
    int i;
    #pragma omp parallel for
    for(i = 0 ; i < rows ;i++){
        double max=0;
        int item;
        for(int j=0;j<columns;j++){
            if(A[i][j]==0.00){
                if(B[i][j]>max){
                    max=B[i][j];
                    item=j;
                }
            }
        }

        fprintf(fp,"%d\n",item);
    }
    fclose(fp);
}