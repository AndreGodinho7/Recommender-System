#include <stdio.h>
#include <stdlib.h>
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
void random_fill_LR(double** L, double** R, int nU, int nI, int nF)
{   
    srandom(0);
    for(int i = 0; i < nU; i++)
        for(int j = 0; j < nF; j++)
            L[i][j] = RAND01 / (double) nF;

    for(int i = 0; i < nF; i++)
        for(int j = 0; j < nI; j++)
            R[i][j] = RAND01 / (double) nF;
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
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < columns; ++j)
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

void matrix_mul(double **firstMatrix, double **secondMatrix, double **result,int nU, int nI, int nF ){
    // result[n1][n3] = firstMatrix[n1][n2] x secondMatrix[n3][n2] (transposed)
    int n1 = nU; // sizeof(matrix_1)/sizeof(matrix_1[0]);
    int n2 = nF; // sizeof(matrix_1[0])/sizeof(matrix_1[0][0]);
    int n3 = nI; //sizeof(matrix_2[0])/sizeof(matrix_2[0][0]);
    double tmp;
    
/*    for(int i = 0; i < n1; i++)
		for(int j = 0; j < n3; j++)
            {
                tmp = 0;
			    for(int k = 0; k < n2; k++)    
                    tmp += firstMatrix[i][k] * secondMatrix[k][j];         
				    result[i][j] = tmp;
            } */
    
    for(int i = 0; i < n1; i++)
		for(int j = 0; j < n3; j++)
            {
                tmp = 0;
			    for(int k = 0; k < n2; k++)
                    tmp += firstMatrix[i][k] * secondMatrix[j][k];       
				    result[i][j] = tmp;
            }
}

double internal_product( double *row ,double *column, int nF){
    double sum=0.0;
    for(int i=0;i<nF;i++){
        sum+=row[i]*column[i];
    }
    return sum;
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
/* TODO 
    - Fazer para todas as iteracoes
    - actualizar pre_L, pre_R etc ...
    - mais algumas verificacoes maybe
*/

void recalculate_Matrix(double** L, double** R,double** pre_L, double** pre_R,double **A,double** B, double** pre_B,int nU, int nI, int nF,int iter, double alpha, _non_zero *v, int non_zero){
    double sum_L = 0;
    double sum_R = 0;
    int k = 0;
    int feature = 0;
    int col = 0;
    
    for(k = 0 ; k < non_zero ;k++){   
        //printf(" linha: %d  coluna : %d\n",v[k].row,v[k].column);
        for(feature = 0; feature < nF; feature ++){
            sum_L = 0;
            sum_R = 0;
                                       
            for(col = 0; col < nI; col++){
                
                if(A[v[k].row][col]!=0){
                    //sum_L += (2*((A[v[k].row][col]-pre_B[v[k].row][col])*(-pre_R[feature][col])));
                    
                    sum_L += (2*((A[v[k].row][col]-pre_B[v[k].row][col])*(-pre_R[col][feature])));
                    //sum_L += (2*((A[v[k].row][col]-internal_product(pre_L[v[k].row],pre_R[col],nF))*(-pre_R[col][feature])));
                }  
            }                                  
            for(int line=0; line < nU ; line++){
                if(A[line][v[k].column]!=0){
                    sum_R += (2*((A[line][v[k].column]-pre_B[line][v[k].column])*(-pre_L[line][feature])));
                    //sum_R += (2*((A[line][v[k].column]-internal_product(pre_L[line],pre_R[v[k].column],nF))*(-pre_L[line][feature])));
                }
            }
            L[v[k].row][feature] = pre_L[v[k].row][feature] - alpha*sum_L;
            //R[feature][v[k].column] = pre_R[feature][v[k].column] - alpha*sum_R;
            R[v[k].column][feature] = pre_R[v[k].column][feature] - alpha*sum_R; //add
        }
    } 

}



void copy_matrix(double** original, double** copied,int rows, int columns){
    for(int i= 0; i<rows;i++){
        for(int j=0;j<columns;j++){
            copied[i][j]=original[i][j];
        }
    }
}

void create_output(double** B,int rows, int columns,char* filename,double** A){
    
    //int size = strlen(filename);
    //printf("%d",size); 

    for(int i = 0 ; i < rows ;i++){
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
    
        printf("%d\n",item);
    }
}