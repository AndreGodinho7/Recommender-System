#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"input.h"
#include "factorization.h"
#include <limits.h>
#include <mpi.h>

#define MASTER_TO_SLAVE_TAG 1 //tag for messages sent from master to slaves
#define SLAVE_TO_MASTER_TAG 20 //tag for messages sent from slaves to master
#define INDEX(row,column,num_column) ((row*num_column)+column)

MPI_Status status; // store status of a MPI_Recv
MPI_Request request; //capture request of a MPI_Isend
MPI_Request request_send_0;
MPI_Request request_send_1;
MPI_Request request_rec_0;
MPI_Request request_rec_1;

int main(int argc, char* argv[])
{   

    MPI_Init(&argc, &argv);


    int id, p;

	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);


    /* create a type for struct typedef _division_mpi */
    const int nitemsV2=5;
    int          blocklengthsV2[5] = {1,1,1,1,1};
    MPI_Datatype typesV2[5] = {MPI_INT, MPI_INT,MPI_INT,MPI_INT,MPI_INT};
    MPI_Datatype mpi_division_slaves;
    MPI_Aint     offsetsV2[5];

    offsetsV2[0] = offsetof(division_mpi, lower_bound);
    offsetsV2[1] = offsetof(division_mpi, upper_bound);
    offsetsV2[2] = offsetof(division_mpi, initial_L);
    offsetsV2[3] = offsetof(division_mpi, final_L);
    offsetsV2[4] = offsetof(division_mpi, num_elements_L);

    MPI_Type_create_struct(nitemsV2, blocklengthsV2, offsetsV2, typesV2, &mpi_division_slaves);
    MPI_Type_commit(&mpi_division_slaves);

    /* create a type for struct _non_zero */
    const int nitems=5;
    int          blocklengths[5] = {1,1,1,1,1};
    MPI_Datatype types[5] = {MPI_INT, MPI_INT,MPI_DOUBLE,MPI_DOUBLE,MPI_INT};
    MPI_Datatype mpi_non_zero;
    MPI_Aint     offsets[5];

    offsets[0] = offsetof(non_zero, row);
    offsets[1] = offsetof(non_zero, column);
    offsets[2] = offsetof(non_zero, A);
    offsets[3] = offsetof(non_zero, B);
    offsets[4] = offsetof(non_zero, process);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_non_zero);
    MPI_Type_commit(&mpi_non_zero);
//====================================================
    int lower_bound;
    int upper_bound; 
    double* L, *R;
    double *L_hold,*R_hold; 

        // auxiliar pointers to avoid copying matrices
    double* L1, *L2;
    double* R1, *R2; 
    double* tmp; 

    double* aux_sum_L;
    double* aux_sum_R;

    // Allocacao de todas as variaveis para a master
    input_values* init ;
    //non_zero* aux;
    division_mpi* slaves; //struct that saves position of L and elements that each processor gets;
    slaves = malloc(p * sizeof(division_mpi));
    int num_Features;
    int num_Items;
    int num_Users;
    double alpha_value;
    int iterations;


    if(id==0){
        
        if (argc != 2){
                printf("ERROR: inserted more than 1 input file.\n");
                exit(0);
        };        
        // for allocating matrice
        init = read_input(argv[1]);
        
        mark_process_in_nonzero(init->num_zeros,init->v,p);



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
        //==============================================================
        //divisao das matrizes para os slaves
        
       
       
        int my_up;
        int my_lower;

     
        my_lower=0;
        my_up=getProcessUpBoundary(init->v,init->num_zeros,id);
        
        printf("Existem %d processos\n",p);
        printf("MASTER : processo %d gets from %d to %d\n",id,my_lower,my_up);
        lower_bound=my_up;
        //int message_tag
        for(int i=1;i<p;i++){
 
            
            
            upper_bound=getProcessUpBoundary(init->v,init->num_zeros,i);

            slaves[i].lower_bound=lower_bound;
            slaves[i].upper_bound=upper_bound;

            printf("MASTER : processo %d gets from %d to %d\n",i,slaves[i].lower_bound,slaves[i].upper_bound);
            fflush(stdout);

            MPI_Isend(&slaves[i].lower_bound, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG +i, MPI_COMM_WORLD, &request);
            MPI_Isend(&slaves[i].upper_bound, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG + 1 +i, MPI_COMM_WORLD, &request);
            MPI_Isend(&init->nF, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG + 2 +i, MPI_COMM_WORLD, &request);
            MPI_Isend(&init->v[slaves[i].lower_bound], (slaves[i].upper_bound - slaves[i].lower_bound), mpi_non_zero, i, MASTER_TO_SLAVE_TAG + 3 +i, MPI_COMM_WORLD, &request);
            
            
            int initial_row=init->v[slaves[i].lower_bound].row;

            int final_row=init->v[slaves[i].upper_bound-1].row + 1;

            int portion_L=(final_row-initial_row)*init->nF; 


            slaves[i].initial_L=INDEX(initial_row,0,init->nF);
            slaves[i].final_L=INDEX(final_row,0,init->nF);
            slaves[i].num_elements_L=portion_L;

            //printf("vou enviar desde %d com %d para %d\n",slaves[i].initial_L,slaves[i].num_elements_L,i);
            //fflush(stdout);

            MPI_Isend(&L1[INDEX(initial_row,0,init->nF)], portion_L, MPI_DOUBLE, i, MASTER_TO_SLAVE_TAG , MPI_COMM_WORLD, &request);
            MPI_Isend(&L2[INDEX(initial_row,0,init->nF)], portion_L, MPI_DOUBLE, i, MASTER_TO_SLAVE_TAG + 1 , MPI_COMM_WORLD, &request);
            
            

            MPI_Isend(&init->iter, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG  , MPI_COMM_WORLD, &request);
            MPI_Isend(&init->alpha, 1, MPI_DOUBLE, i, MASTER_TO_SLAVE_TAG + 2 , MPI_COMM_WORLD, &request); 
            MPI_Isend(&init->nU, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG + 3 , MPI_COMM_WORLD, &request);
            MPI_Isend(&init->nI, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG , MPI_COMM_WORLD, &request);
            MPI_Isend(&R[0], init->nF*init->nI , MPI_DOUBLE, i, MASTER_TO_SLAVE_TAG +2 , MPI_COMM_WORLD, &request);
            MPI_Isend(&R_hold[0], init->nF*init->nI, MPI_DOUBLE, i, MASTER_TO_SLAVE_TAG + 3 , MPI_COMM_WORLD, &request);
            
                    
            lower_bound=upper_bound;
            
            
        }

        lower_bound=my_lower;
        upper_bound=my_up;  
        int initial_row=init->v[lower_bound].row;
        int final_row=init->v[upper_bound-1].row + 1;
        int portion_L=(final_row-initial_row)*init->nF;

        slaves[0].lower_bound=lower_bound;
        slaves[0].upper_bound=upper_bound;        
        slaves[0].num_elements_L=portion_L;
        slaves[0].initial_L=INDEX(initial_row,0,init->nF);
        slaves[0].final_L=INDEX(final_row,0,init->nF);

        for(int i = 1 ; i<p ;i++){ // master process receives all results
            MPI_Isend(&slaves[0], p, mpi_division_slaves, i, MASTER_TO_SLAVE_TAG , MPI_COMM_WORLD, &request);
            
        }  
        num_Users=init->nU;
        num_Items=init->nI;
        num_Features=init->nF;
        alpha_value=init->alpha;
        iterations=init->iter;


         
    }
    if(id>0){
    
        MPI_Recv(&lower_bound, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG +id, MPI_COMM_WORLD, &status);
        MPI_Recv(&upper_bound, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG + 1+id, MPI_COMM_WORLD, &status);
        MPI_Recv(&num_Features, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG + 2 +id, MPI_COMM_WORLD, &status);
        init = (input_values*)malloc(sizeof(input_values));
        init->v = malloc((upper_bound-lower_bound) * sizeof(non_zero));
        
        
        MPI_Recv(&init->v[0], (upper_bound - lower_bound) , mpi_non_zero, 0, MASTER_TO_SLAVE_TAG + 3+id, MPI_COMM_WORLD, &status);
        /*printf("V array received in process %d\n",id);
        fflush(stdout);*/
        printf("lower bound = %d upper_bound=%d on process %d\n",lower_bound,upper_bound,id);
        fflush(stdout);

        int initial_row=init->v[0].row;

        int final_row=init->v[upper_bound - lower_bound-1].row + 1 ;

        int portion_L=(final_row-initial_row)*num_Features;

        printf("Estou a ver da linha %d ate %d no processo %d\n",initial_row,final_row,id);
        fflush(stdout);
        int init_row=final_row-initial_row;
        L = MatrixInit(init_row, num_Features);
        L_hold = MatrixInit(init_row, num_Features);

        MPI_Recv(&L[0], portion_L , MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG , MPI_COMM_WORLD, &status);
        MPI_Recv(&L_hold[0], portion_L , MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG + 1, MPI_COMM_WORLD, &status);
        printf("L matrix received in process %d\n",id);
        fflush(stdout);

        MPI_Recv(&iterations, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG  , MPI_COMM_WORLD, &status);
        MPI_Recv(&alpha_value, 1, MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG + 2 , MPI_COMM_WORLD, &status);
        MPI_Recv(&num_Users, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG + 3 , MPI_COMM_WORLD, &status);
        MPI_Recv(&num_Items, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG , MPI_COMM_WORLD, &status);


        R = MatrixInit( num_Items,num_Features);
        R_hold = MatrixInit(  num_Items,num_Features);
        
        
        MPI_Recv(&R[0], num_Items*num_Features , MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG + 2, MPI_COMM_WORLD, &status);
        MPI_Recv(&R_hold[0], num_Items*num_Features , MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG + 3, MPI_COMM_WORLD, &status);
        MPI_Recv(&slaves[0], p , mpi_division_slaves, 0, MASTER_TO_SLAVE_TAG , MPI_COMM_WORLD, &status);
        printf("slave lower = %d slaves upper =%d on process %d\n",slaves[id].lower_bound,slaves[id].upper_bound,id);
        //printf("R matrix received in process %d\n",id);
        fflush(stdout);       
        L1 = L;
        L2 = L_hold; 
        R1 = R;
        R2 = R_hold;  
       
        matrix_mul(L1,R1,init->v,upper_bound-lower_bound,num_Features);
     
        MPI_Isend(&lower_bound, 1, MPI_INT, 0, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD, &request);
        MPI_Isend(&upper_bound, 1, MPI_INT, 0, SLAVE_TO_MASTER_TAG + 1, MPI_COMM_WORLD, &request);
        MPI_Isend(&init->v[0], (upper_bound - lower_bound), mpi_non_zero, 0, SLAVE_TO_MASTER_TAG + 2, MPI_COMM_WORLD, &request);
        

           
             
               
    }
    
    if(id==0){
        matrix_mul(L1, R1, init->v, upper_bound-lower_bound , init->nF); 
        for(int i = 1 ; i<p ;i++){ // master process receives all results
            MPI_Recv(&lower_bound, 1, MPI_INT, i, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD, &status);
            MPI_Recv(&upper_bound, 1, MPI_INT, i, SLAVE_TO_MASTER_TAG + 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&init->v[lower_bound], (upper_bound - lower_bound) , mpi_non_zero, i, SLAVE_TO_MASTER_TAG + 2, MPI_COMM_WORLD, &status);
        
        }   

        
        aux_sum_L= MatrixInit(init->nU, init->nF); // auxiliar matriz to sum matrix L
        aux_sum_R= MatrixInit(init->nI, init->nF); // auxiliar matriz to sum matrix R
        
        //printMatrix(L1,init->nU, init->nF);
    }

    int *displs = (int *)malloc(p*sizeof(int)); 
    int *scounts = (int *)malloc(p*sizeof(int));
    int *rcounts = (int *)malloc(p*sizeof(int));  
    for (int i=0; i<p; i++) { 
        displs[i] = slaves[i].initial_L; 
        scounts[i] = slaves[i].num_elements_L; 
        rcounts[i] = slaves[i].num_elements_L;
    } 
    rcounts[0]=num_Users*num_Features;

    int zeros_count;
    
       
    /*printf("process %d\n",id);
    fflush(stdout);
    printf("lower = %d\n",slaves[id].lower_bound);
    fflush(stdout);
    printf("upper = %d\n",slaves[id].upper_bound);
    fflush(stdout);
    printf("initial_L = %d\n",slaves[id].initial_L);
    fflush(stdout);
    printf("final_L = %d\n",slaves[id].final_L);
    fflush(stdout);
    printf("portion = %d\n",slaves[id].num_elements_L);
    fflush(stdout);*/
        
    
    int message_tag=0;
    int user_portion;
    for(int i = 0 ; i < iterations ; i++){
        //update the matrix
        /*if(i==250){
            MPI_Barrier(MPI_COMM_WORLD);
            break;
        }*/
        tmp = L1;
        L1 = L2;
        L2 = tmp;

        tmp = R1;
        R1 = R2;
        R2 = tmp; 
       
        message_tag=(i*2);
        zeros_count= slaves[id].upper_bound - slaves[id].lower_bound;


        if(id==0){
            user_portion=num_Users;
        }
        if(id>0){
            int initial_row=init->v[0].row;
            int final_row=init->v[zeros_count-1].row + 1 ;
            user_portion= final_row-initial_row;  

        }
        recalculate_Matrix(L1, R1, L2, R2, user_portion, num_Items,num_Features,alpha_value, init->v ,zeros_count ,id,p,slaves);

        MPI_Gatherv(&L1[0], slaves[id].num_elements_L, MPI_DOUBLE, &aux_sum_L[0], scounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if(id==0){
            //printMatrix(aux_sum_L,num_Users,num_Features);
            for(int u = 0; u < num_Users; u++){               
                for(int f = 0; f < num_Features; f++){           
                    aux_sum_L[INDEX(u,f,num_Features)] = L2[INDEX(u,f,num_Features)] + alpha_value*2*aux_sum_L[INDEX(u,f,num_Features)];              

                }
            }
        }

        
        MPI_Reduce(&R1[0], &aux_sum_R[0], num_Items*num_Features, MPI_DOUBLE, MPI_SUM, 0,MPI_COMM_WORLD);

        if(id==0){
            //printMatrix(aux_sum_R,num_Items,num_Features);
            for(int k = 0; k < num_Items; k++){
                for(int f = 0; f < num_Features; f++){
                    R1[INDEX(k,f,num_Features)] = R2[INDEX(k,f,num_Features)] + alpha_value*2*aux_sum_R[INDEX(k,f,num_Features)]; 
                } 
            }
            /*if(iterations-i < 2){
                printf("Matrix L na iteracao %d\n",i);
                printMatrix(aux_sum_L,num_Users,num_Features);
                printf("Matrix R na iteracao %d\n",i);     
                printMatrix(R1,num_Items,num_Features);
            }*/

                
        }        
        

        MPI_Scatterv(&aux_sum_L[0], rcounts, displs, MPI_DOUBLE, &L1[0] , rcounts[id], MPI_DOUBLE, 0, MPI_COMM_WORLD);
 
        MPI_Bcast(&R1[0],num_Items*num_Features , MPI_DOUBLE, 0, MPI_COMM_WORLD );

        matrix_mul(L1,R1,init->v,slaves[id].upper_bound-slaves[id].lower_bound,num_Features);

        if(id>0){
           MPI_Isend(&init->v[0], (slaves[id].upper_bound-slaves[id].lower_bound), mpi_non_zero, 0, SLAVE_TO_MASTER_TAG + message_tag, MPI_COMM_WORLD, &request); 
        }
        if(id==0){
            for(int k = 1 ; k<p ;k++){ // master process receives all results
                MPI_Recv(&init->v[slaves[k].lower_bound], (slaves[k].upper_bound-slaves[k].lower_bound) , mpi_non_zero, k, SLAVE_TO_MASTER_TAG +message_tag, MPI_COMM_WORLD, &status);
        
            }             
        }


        
        /*printf("Matrix L na iteracao %d\n",i);
        printMatrix(L1,init->nU, init->nF); 
        printf("Matrix R na iteracao %d\n",i);       
        printMatrix(R1,init->nI, init->nF); */
        

    }
    if(id==0){
        create_output(init->v, init->nU, init->nI, init->nF, L1, R1, init->num_zeros);
        free(init->v);
        free(init);
        
    }
    free(L1);
    free(L2);
    free(slaves);
    free(R1);
    free(R2);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize ();
    return 0;
}