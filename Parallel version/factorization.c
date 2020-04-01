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
    int i, j;
        
    for(i = 0; i < nU; i++)
        for(j = 0; j < nF; j++)
            L_[i][j] = RAND01 / (double) nF;
    
    for(i = 0; i < nF; i++)
        for(j = 0; j < nI; j++)
            R_[i][j] = RAND01 / (double) nF;
}

/******************************************************************************
* transpose()
*
* Arguments: matrix - pointer to pointer of matrix to be tranposed
*            rows - number of rows of matrix
*            columns - number of columns of matrix
*            
* Returns: double pointer to transposed matrix
*										
* Side-Effects: 
*              
*
* Description: transposes matrix
*              
*****************************************************************************/

double** transpose(double** matrix, int rows, int columns)
{
    double** result = MatrixInit(columns, rows);
    int i,j;
    #pragma omp parallel for private(i,j)
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
*            v - vector which has elements of B that are non-zero
*            nF - number of features
*
* Returns: void
*										
* Side-Effects: 
*
* Description: calculates one element of matrix B by multiplying a row of firstmatrix and 
*              a column of secondmatrix
*
*****************************************************************************/

void matrix_mul(double **firstMatrix, double **secondMatrix, non_zero* v, int num_zeros ,int nF){
    int z, k, i, j;
    double sum = 0;

    // dúvida: inner for k não é dividido entre threads
    #pragma omp for private(z,k,i,j) 
    for (z = 0; z < num_zeros; z++){
        i = v[z].row;
        j = v[z].column;
        sum = 0;
        for (k = 0; k < nF; k++)
            sum += firstMatrix[i][k]*secondMatrix[j][k];
        v[z].B = sum;
    }
}

/******************************************************************************
* zero_LR()
*
* Arguments: L - matrix L 
*            R - matrix R
*            nU - number of users          
*            nI - number of items
*            nF - number of features
*
* Returns: void
*										
* Side-Effects: 
*
* Description: puts elements of matrices L and R as 0
*
*****************************************************************************/
void zero_LR(double** L, double** R, int nU, int nI, int nF){
    int u,i,f;
    
    #pragma omp for private(f) nowait
    for(u = 0; u < nU; u++)
        for(f = 0; f < nF; f++)
            L[u][f] = 0;

    #pragma omp for private(f) 
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
*            nU - number of users          
*            nI - number of items
*            nF - number of features
*            iter - maximum number of iterations
*            alpha - converge rate
*            v - array of type _non_zero with all information about non zero values in matrix A and B
*            num_zeros - number of non zero values in matrix A
*
* Returns: void
*										
* Side-Effects: 
*
* Description: computes the algorithm (minimizing the difference between A and B)
*
*****************************************************************************/

void recalculate_Matrix(double** L, double** R, double** pre_L, double** pre_R, int nU, int nI, int nF, double alpha, non_zero *v, int num_zeros){
    int i, j, z, f, k;
    double a, b;

    zero_LR(L, R, nU, nI, nF);
    
    // for (z = 0; z < num_zeros; z++)
    //     printf("%lf ", v[z].B);

    #pragma omp for private(z,a,b,i,j,k)
    for (z = 0; z < num_zeros; z++){
            i = v[z].row;
            j = v[z].column;
            a = v[z].A;
            b = v[z].B;
        for(k = 0; k < nF; k++){
           #pragma omp atomic
            L[i][k] += (a-b)*(pre_R[j][k]);
           #pragma omp atomic
            R[j][k] += (a-b)*(pre_L[i][k]);
        }
    }
    
    #pragma omp for private(f) nowait
    for(int u = 0; u < nU; u++)
        for(f = 0; f < nF; f++)
            L[u][f] = pre_L[u][f] + alpha*2*L[u][f];

    #pragma omp for private(f) 
    for(int i = 0; i < nI; i++)
        for(f = 0; f < nF; f++)
            R[i][f] = pre_R[i][f] + alpha*2*R[i][f]; 
}

/******************************************************************************
* create_output()
*
* Arguments: B - final matrix B
*            A - matrix A
*            rows - number of rows of matrix B
*            columns - number of columns of matrix B
*
* Returns: void
*										
* Side-Effects: 
*
* Description: creates output file with recommendations
*
*****************************************************************************/

void create_output(double** B,int rows, int columns,char* filename,double** A){
    FILE* fp = fopen("recsystem.out", "w");
    //int size = strlen(filename);
    //printf("%d",size); 
    int i;
    //#pragma omp parallel for
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