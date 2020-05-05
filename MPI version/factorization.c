#include <stdio.h>
#include <stdlib.h>
#include "input.h"
#include <float.h>
#include <mpi.h>

#define RAND01 ((double) random() / (double) RAND_MAX)

#define MASTER_TO_SLAVE_TAG 1 //tag for messages sent from master to slaves
#define SLAVE_TO_MASTER_TAG 4 //tag for messages sent from slaves to master

MPI_Status status; // store status of a MPI_Recv
MPI_Request request; //capture request of a MPI_Isend
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
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < columns; ++j)
            result[j][i] = matrix[i][j];

    for (int k = 0; k < rows; k++)
    {
        free(matrix[k]);
        
    }
    free(matrix);

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
void matrix_mul_mpi(double **firstMatrix, double **secondMatrix, non_zero* v, int num_zeros,int nF,int id, int p){
    
    
    const int nitems=4;
    int          blocklengths[4] = {1,1,1,1};
    MPI_Datatype types[4] = {MPI_INT, MPI_INT,MPI_DOUBLE,MPI_DOUBLE};
    MPI_Datatype mpi_non_zero;
    MPI_Aint     offsets[4];

    offsets[0] = offsetof(non_zero, row);
    offsets[1] = offsetof(non_zero, column);
    offsets[2] = offsetof(non_zero, A);
    offsets[3] = offsetof(non_zero, B);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_non_zero);
    MPI_Type_commit(&mpi_non_zero);

    int length_vector=num_zeros;
    //int left_overs=length_vector%p;
    int portion = length_vector/(p-1);
    //int display[p];    //For MPI_GATHERV() if used
    //int recv_counts[p];
    int lower_bound;
    int upper_bound;

    if(id==0){ //master process
        //display[0]=0; //For MPI_GATHERV() if used
        //recv_counts[0]=0;
        for(int i = 1 ; i<p ;i++){
            lower_bound=(i-1)*portion;
            if((i+1)==p && (length_vector%(p-1) !=0)){
                upper_bound=length_vector;
            }
            else{
                upper_bound=lower_bound+portion;
            }
            //display[i]=lower_bound;
            //recv_counts[i]=upper_bound-lower_bound;
            printf("process %d is sending to %d\n", id,i);
            MPI_Isend(&lower_bound, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &request);
            MPI_Isend(&upper_bound, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG + 1, MPI_COMM_WORLD, &request);
            MPI_Isend(&v[lower_bound], (upper_bound - lower_bound), mpi_non_zero, i, MASTER_TO_SLAVE_TAG + 2, MPI_COMM_WORLD, &request);
        }
    }

    if(id>0){ //slave processors

        MPI_Recv(&lower_bound, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&upper_bound, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG + 1, MPI_COMM_WORLD, &status);
        non_zero* aux = malloc((upper_bound-lower_bound) * sizeof(non_zero));
        MPI_Recv(&aux[0], (upper_bound - lower_bound) , mpi_non_zero, 0, MASTER_TO_SLAVE_TAG + 2, MPI_COMM_WORLD, &status);
        
        printf("process %d has %d elements to compute\n",id,upper_bound-lower_bound);
        //printf("vector of non zero v starts at row %d, column %d\n",init->v[lower_bound].row,init->v[lower_bound].column);
    
    for (int z = 0; z < upper_bound - lower_bound; z++){
        int i = aux[z].row;
        int j = aux[z].column;
        aux[z].B = 0;
        for (int k = 0; k < nF; k++)
        aux[z].B += firstMatrix[i][k]*secondMatrix[j][k];
    }
        MPI_Isend(&lower_bound, 1, MPI_INT, 0, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD, &request);
        MPI_Isend(&upper_bound, 1, MPI_INT, 0, SLAVE_TO_MASTER_TAG + 1, MPI_COMM_WORLD, &request);
        MPI_Isend(&aux[0], (upper_bound - lower_bound), mpi_non_zero, 0, SLAVE_TO_MASTER_TAG + 2, MPI_COMM_WORLD, &request);
        
        free(aux);
    }
    if(id==0){  // so esta aqui porque so a parte de cima esta em paralelo, e para nao dar 4 vezes o resultado
        for(int i = 1 ; i<p ;i++){ // master process receives all results
            MPI_Recv(&lower_bound, 1, MPI_INT, i, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD, &status);
            MPI_Recv(&upper_bound, 1, MPI_INT, i, SLAVE_TO_MASTER_TAG + 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&v[lower_bound], (upper_bound - lower_bound) , mpi_non_zero, i, SLAVE_TO_MASTER_TAG + 2, MPI_COMM_WORLD, &status);
        
        }
    
    }
    MPI_Barrier(MPI_COMM_WORLD);


}
void matrix_mul(double **firstMatrix, double **secondMatrix, non_zero* v, int num_zeros ,int nF){
    for (int z = 0; z < num_zeros; z++){
        int i = v[z].row;
        int j = v[z].column;
        v[z].B = 0;
        for (int k = 0; k < nF; k++){
            v[z].B += firstMatrix[i][k]*secondMatrix[j][k];
        }
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
    for(int u = 0; u < nU; u++)
        for(int f = 0; f < nF; f++)
            L[u][f] = 0;

    for(int i = 0; i < nI; i++)
        for(int f = 0; f < nF; f++)
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


void recalculate_Matrix(double** L, double** R, double** pre_L, double** pre_R, int nU, int nI, int nF, double alpha, non_zero *v, int num_zeros,int id){
    int i, j, z;
    double a, b;

    zero_LR(L, R, nU, nI, nF);

    //printf("o processador no recalculate que esta aqui e %d\n",id);
    for (z = 0; z < num_zeros; z++){
        i = v[z].row;
        j = v[z].column;
        a = v[z].A;
        b = v[z].B;

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

void create_output(non_zero *v, int nU, int nI, int nF, double** L, double** R, int num_zeros){
    
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
                element += L[i][k] * R[j][k];

            if (element > max){
                max = element;
                position = j;
            }
        }
        printf("%d\n",position);
    }
    
}