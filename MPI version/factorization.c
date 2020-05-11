#include <stdio.h>
#include <stdlib.h>
#include "input.h"
#include <float.h>
#include <mpi.h>

#define RAND01 ((double) random() / (double) RAND_MAX)
#define INDEX(row,column,num_column) ((row*num_column)+column)






int find_upper_bound(int lower_row,int upper_row,int lower_bound,non_zero *v,int num_zeros){
    
    for(int i =lower_bound; i<num_zeros;i++){
        
        if(v[i].row==upper_row){
            
            return i;
        }
    }
    return num_zeros ;
}





/***************************************
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
***************************************/
void printMatrix(double* matrix, int rows, int columns)
{
    for (int i = 0; i < rows; i++)
    {
        for(int j = 0; j < columns; j++)
                
            printf("%f     ", matrix[INDEX(i,j,columns)]);
                

        printf("\n");
    }
    printf("\n");
}

/***************************************
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
***************************************/

double* MatrixInit(int rows, int columns)
{   
    double* matrix; 

    matrix = (double *)calloc(rows*columns, sizeof(double *)); 
    if (matrix == NULL){
      printf("ERROR: Out of memory.\n");
		  exit(1);
    }


    return matrix;
}


/***************************************
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
***************************************/
void random_fill_LR(double* L_, double* R_, int nU, int nI, int nF)
{   
    srandom(0);
    for(int i = 0; i < nU; i++)
        for(int j = 0; j < nF; j++)
            L_[INDEX(i,j,nF)] = RAND01 / (double) nF;

    for(int i = 0; i < nF; i++)
        for(int j = 0; j < nI; j++)
            R_[INDEX(i,j,nI)] = RAND01 / (double) nF;
}


/***************************************
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
***************************************/

double* transpose(double* matrix, int rows, int columns)
{   
    double* result = MatrixInit(columns, rows);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < columns; ++j)
            result[INDEX(j,i,rows)] = matrix[INDEX(i,j,columns)];


    free(matrix);

    return result;
}


/***************************************
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
***************************************/

void matrix_mul(double *firstMatrix, double *secondMatrix, non_zero* v, int num_zeros ,int nF){
    int start= v[0].row;
    for (int z = 0; z < num_zeros; z++){
        int i = v[z].row;
        int j = v[z].column;
        int pos_row=0;
        v[z].B = 0;
        for (int k = 0; k < nF; k++){
            pos_row=i-start;
            
            if(firstMatrix[INDEX(pos_row,k,nF)]*secondMatrix[INDEX(j,k,nF)]<0){
                printf("ZEROO - com L= %f and R = %f\n",firstMatrix[INDEX(pos_row,k,nF)],secondMatrix[INDEX(j,k,nF)]);
                
                
            }
        
            v[z].B += firstMatrix[INDEX(pos_row,k,nF)]*secondMatrix[INDEX(j,k,nF)];
        }
    }
}

/***************************************
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
***************************************/
void zero_LR(double* L, double* R, int nU, int nI, int nF){
    for(int u = 0; u < nU; u++)
        for(int f = 0; f < nF; f++)
            L[INDEX(u,f,nF)] = 0;

    for(int i = 0; i < nI; i++)
        for(int f = 0; f < nF; f++)
            R[INDEX(i,f,nF)] = 0;

}

/***************************************
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
***************************************/


void recalculate_Matrix(double* L, double* R, double* pre_L, double* pre_R, int nU, int nI, int nF, double alpha, non_zero *v, int num_zeros,int id,int p, division_mpi *slaves){
    int i, j, z;
    double a, b;
 
    zero_LR(L, R, nU, nI, nF);
    //printf("PUS A ZERO no processo %d\n",id);
    int start= v[0].row;
    //printf("O processo %d vai começar em == %d\n",id,start);
    int pos_row;
    for (z = 0; z < num_zeros; z++){
        i = v[z].row;
        j = v[z].column;
        a = v[z].A;
        b = v[z].B;

        for(int k = 0; k < nF; k++){
            pos_row=i-start;
            if(pre_R[INDEX(j,k,nF)]<0){
                printf("\n\nSomething Wrong R %f\n\n",pre_R[INDEX(j,k,nF)]);
            }
            if(pre_L[INDEX(pos_row,k,nF)]<0){
                printf("\n\nSomething Wrong L %f\n\n",pre_L[INDEX(pos_row,k,nF)]);
            }
            L[INDEX(pos_row,k,nF)] += (a-b)*(pre_R[INDEX(j,k,nF)]);
            R[INDEX(j,k,nF)] += (a-b)*(pre_L[INDEX(pos_row,k,nF)]);
        }
    }

    /*for(int u = 0; u < nU; u++)
        for(int f = 0; f < nF; f++)
            L[INDEX(u,f,nF)] = pre_L[INDEX(u,f,nF)] + alpha*2*L[INDEX(u,f,nF)];


    for(int i = 0; i < nI; i++)
        for(int f = 0; f < nF; f++)
            R[INDEX(i,f,nF)] = pre_R[INDEX(i,f,nF)] + alpha*2*R[INDEX(i,f,nF)]; */

}


/***************************************
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
***************************************/

void create_output(non_zero *v, int nU, int nI, int nF, double* L, double* R, int num_zeros){
    
    int i,j,k;
    int z = 0;
    double element;
    int position;

    for(i = 0 ; i < nU ;i++){
        double max = -DBL_MAX;
        element = 0;
        position = -1;

        for(j = 0 ; j < nI ;j++){
            if(v[z].row == i && v[z].column == j && z < num_zeros){
                z++;
                continue;
            }

            element = 0;
            for (k = 0; k < nF; k++)
                element += L[INDEX(i,k,nF)] * R[INDEX(j,k,nF)];

            if (element > max){
                max = element;
                position = j;
            }
        }
        printf("%d\n",position);
    }
    
}