#include <stdio.h>
#include <stdlib.h>

#include"input.h"
#include "factorization.h"
#include <limits.h>
#include <mpi.h>

#define MASTER_TO_SLAVE_TAG 1 //tag for messages sent from master to slaves
#define SLAVE_TO_MASTER_TAG 4 //tag for messages sent from slaves to master

MPI_Status status; // store status of a MPI_Recv
MPI_Request request; //capture request of a MPI_Isend

int main(int argc, char* argv[])
{   
    MPI_Init(&argc, &argv);


    int id, p;

	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

    /* create a type for struct _non_zero */
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

/*==========================================*/

    input_values* init;
    
    


    // for allocating matrices
    double** L, **R;
    double **L_hold,**R_hold; 

    // auxiliar pointers to avoid copying matrices
    double** L1, **L2;
    double** R1, **R2; 
    double** tmp;

    if (argc != 2){
        printf("ERROR: inserted more than 1 input file.\n");
        exit(0);
    };

    init = read_input(argv[1]);
    L = MatrixInit(init->nU, init->nF);
    R = MatrixInit(init->nF, init->nI);
    L_hold = MatrixInit(init->nU, init->nF); // Matrix that stores the previous iteration of L
    R_hold = MatrixInit(init->nF, init->nI); // Matrix that stores the previous iteration of R
    random_fill_LR(L, R, init->nU, init->nI, init->nF);

    L1 = L;
    L2 = L_hold;

    R = transpose(R, init->nF, init->nI); 
    R_hold = transpose(R_hold, init->nF, init->nI); 

    R1 = R;
    R2 = R_hold;

    /*MPI division for all slaves */
    int length_vector=init->num_zeros;
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
            MPI_Isend(&init->v[lower_bound], (upper_bound - lower_bound), mpi_non_zero, i, MASTER_TO_SLAVE_TAG + 2, MPI_COMM_WORLD, &request);
        }
    }

    if(id>0){ //slave processors

        MPI_Recv(&lower_bound, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&upper_bound, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG + 1, MPI_COMM_WORLD, &status);
        non_zero* aux = malloc((upper_bound-lower_bound) * sizeof(non_zero));
        MPI_Recv(&aux[0], (upper_bound - lower_bound) , mpi_non_zero, 0, MASTER_TO_SLAVE_TAG + 2, MPI_COMM_WORLD, &status);
        
        printf("process %d has %d elements to compute\n",id,upper_bound-lower_bound);
        printf("vector of non zero v starts at row %d, column %d\n",init->v[lower_bound].row,init->v[lower_bound].column);
        
        matrix_mul(L1, R1, aux,(upper_bound - lower_bound) , init->nF);
        
        MPI_Isend(&lower_bound, 1, MPI_INT, 0, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD, &request);
        MPI_Isend(&upper_bound, 1, MPI_INT, 0, SLAVE_TO_MASTER_TAG + 1, MPI_COMM_WORLD, &request);
        MPI_Isend(&aux[0], (upper_bound - lower_bound), mpi_non_zero, 0, SLAVE_TO_MASTER_TAG + 2, MPI_COMM_WORLD, &request);
        
        free(aux);
    }

    

    
    //MPI_Gatherv(init->v, rank, MPI_INT, init->v, recv_counts, display, MPI_INT, 0, MPI_COMM_WORLD);
    
    if(id==0){  // so esta aqui porque so a parte de cima esta em paralelo, e para nao dar 4 vezes o resultado
        for(int i = 1 ; i<p ;i++){ // master process receives all results
            MPI_Recv(&lower_bound, 1, MPI_INT, i, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD, &status);
            MPI_Recv(&upper_bound, 1, MPI_INT, i, SLAVE_TO_MASTER_TAG + 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&init->v[lower_bound], (upper_bound - lower_bound) , mpi_non_zero, i, SLAVE_TO_MASTER_TAG + 2, MPI_COMM_WORLD, &status);
        
        }
        /*printf("matrix L\n");
        printMatrix(L1,init->nU, init->nF);
        printf("matrix R\n");
        printMatrix(R1,init->nI, init->nF);

        printf("Depois da mul\n\n");
        for (int z = 0; z < init->num_zeros; z++){
            int i = init->v[z].row;
            int j = init->v[z].column;
            printf("row: %d  column:%d  value:%f \n",i,j,init->v[z].B);
            
            
        }*/
    
    
    for(int i = 0 ; i < init->iter ; i++){
        /*update the matrix*/
            tmp = L1;
            L1 = L2;
            L2 = tmp;

            tmp = R1;
            R1 = R2;
            R2 = tmp; 
        recalculate_Matrix(L1, R1, L2, R2, init->nU, init->nI, init->nF, init->alpha,init->v ,init->num_zeros);
        matrix_mul(L1, R1, init->v, init->num_zeros, init->nF);
    }
    
/*     for(int i = 0 ; i < init->iter ; i++){
        
        if(i==10)
            break;
        if(id==0){
            tmp = L1;
            L1 = L2;
            L2 = tmp;

            tmp = R1;
            R1 = R2;
            R2 = tmp; 
            recalculate_Matrix(L1, R1, L2, R2, init->nU, init->nI, init->nF, init->alpha,init->v ,init->num_zeros,id);
            printf("\nRECALCULEI A MATRIX\n");
            fflush(stdout);
            for(int j = 1 ; j<p ;j++){
                lower_bound=(j-1)*portion;
                if((j+1)==p && (length_vector%(p-1) !=0)){
                    upper_bound=length_vector;
                }
                else{
                    upper_bound=lower_bound+portion;
                }
                //display[i]=lower_bound;
                //recv_counts[i]=upper_bound-lower_bound;
                //printf("process %d is sending to %d\n", id,i);
                MPI_Send(&lower_bound, 1, MPI_INT, j, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD);
                MPI_Send(&upper_bound, 1, MPI_INT, j, MASTER_TO_SLAVE_TAG + 1, MPI_COMM_WORLD);
                MPI_Send(&init->v[lower_bound], (upper_bound - lower_bound), mpi_non_zero, j, MASTER_TO_SLAVE_TAG + 2, MPI_COMM_WORLD);
                printf("ENVIEI INICIAL para %d non zero v row %d, column %d value=%f \n",j,init->v[lower_bound].row,init->v[lower_bound].column,init->v[lower_bound].B);
                fflush(stdout);
            }

        }
        if(id>0){ //slave processors

            MPI_Recv(&lower_bound, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);
            MPI_Recv(&upper_bound, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG + 1, MPI_COMM_WORLD, &status);
            non_zero* aux = malloc((upper_bound-lower_bound) * sizeof(non_zero));
            MPI_Recv(&aux[0], (upper_bound - lower_bound) , mpi_non_zero, 0, MASTER_TO_SLAVE_TAG + 2, MPI_COMM_WORLD, &status);
            printf("RECEBI INICIAL de 0 e sou o %d para 0 non zero v row %d, column %d value=%f \n",id,aux[0].row,aux[0].column,aux[0].B);
            //printf("process %d has %d elements to compute\n",id,upper_bound-lower_bound);
            //printf("non zero v row %d, column %d  and aux row %d, column %d \n",init->v[lower_bound].row,init->v[lower_bound].column,aux[0].row,aux[0].column);
            
            matrix_mul(L1, R1, aux,(upper_bound - lower_bound) , init->nF);
            
            MPI_Send(&lower_bound, 1, MPI_INT, 0, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD);
            MPI_Send(&upper_bound, 1, MPI_INT, 0, SLAVE_TO_MASTER_TAG + 1, MPI_COMM_WORLD);
            MPI_Send(&aux[0], (upper_bound - lower_bound), mpi_non_zero, 0, SLAVE_TO_MASTER_TAG + 2, MPI_COMM_WORLD);
            printf("ENVIEI de %d para 0 non zero v row %d, column %d value=%f \n",id,aux[0].row,aux[0].column,aux[0].B);
            fflush(stdout);
            free(aux);
        }
        if(id==0){  // so esta aqui porque so a parte de cima esta em paralelo, e para nao dar 4 vezes o resultado
            for(int k = 1 ; k<p ;k++){ // master process receives all results
                MPI_Recv(&lower_bound, 1, MPI_INT, k, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD, &status);
                MPI_Recv(&upper_bound, 1, MPI_INT, k, SLAVE_TO_MASTER_TAG + 1, MPI_COMM_WORLD, &status);
                MPI_Recv(&init->v[lower_bound], (upper_bound - lower_bound) , mpi_non_zero, k, SLAVE_TO_MASTER_TAG + 2, MPI_COMM_WORLD, &status);
                printf("RECEBI de %d non zero v row %d, column %d value=%f \n",k,init->v[lower_bound].row,init->v[lower_bound].column,init->v[lower_bound].B);
                fflush(stdout);
            }
        }
        //matrix_mul(L1, R1, init->v,init->num_zeros , init->nF);
        //matrix_mul_mpi(L1, R1, init->v,init->num_zeros , init->nF,id,p);
    } */
    //if(id==0){
        /**for (int z = 0; z < init->num_zeros; z++){
            int i = init->v[z].row;
            int j = init->v[z].column;
            printf("row: %d  column:%d  value:%f \n",i,j,init->v[z].B);
            
            
        }*/
        create_output(init->v, init->nU, init->nI, init->nF, L1, R1, init->num_zeros);
    }
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
    free(init->v);
    free(init);
    
    MPI_Finalize ();
    return 0;
}