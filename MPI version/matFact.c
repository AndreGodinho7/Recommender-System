#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"input.h"
#include "factorization.h"
#include <limits.h>
#include <mpi.h>

#define MASTER_TO_SLAVE_TAG 1 //tag for messages sent from master to slaves
#define SLAVE_TO_MASTER_TAG 5 //tag for messages sent from slaves to master
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
    non_zero* aux;
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

        mark_process_in_nonzero(init->num_zeros, init->v, p);
        

        L = MatrixInit(init->nU, init->nF);
        //printf("FIZ primerio init\n");
        R = MatrixInit(init->nF, init->nI);
        L_hold = MatrixInit(init->nU, init->nF); // Matrix that stores the previous iteration of L
        R_hold = MatrixInit(init->nF, init->nI); // Matrix that stores the previous iteration of R
        random_fill_LR(L, R, init->nU, init->nI, init->nF);
        //printf("FIZ random fill\n");
        L1 = L;
        L2 = L_hold;

        R = transpose(R, init->nF, init->nI); 
        R_hold = transpose(R_hold, init->nF, init->nI); 
        //printf("FIZ transpose\n");
        R1 = R;
        R2 = R_hold;
        //==============================================================
        //divisao das matrizes para os slaves
        int portion ;
        double num_ceil;
        int my_up;
        int my_lower;
        int upper_row;
        
        mark_process_in_nonzero(init->num_zeros, init->v, p);
        my_up = getProcessUpBoundary(init->v, init->num_zeros, 0);
        my_lower = 0;

        printf("MASTER : processo %d gets from %d to %d\n",id,my_lower,my_up);
        lower_bound=my_up;

        int curr_p=3;
        for(int i=1;i<p;i++){
            upper_bound = getProcessUpBoundary(init->v, init->num_zeros, i);

            slaves[i].lower_bound=lower_bound;
            slaves[i].upper_bound=upper_bound;

            printf("MASTER : processo %d gets from %d to %d\n",i,lower_bound,upper_bound);
            fflush(stdout);
            MPI_Isend(&lower_bound, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &request);
            MPI_Isend(&upper_bound, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG + 1, MPI_COMM_WORLD, &request);
            MPI_Isend(&init->nF, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG + 2, MPI_COMM_WORLD, &request);
            MPI_Isend(&init->v[lower_bound], (upper_bound - lower_bound), mpi_non_zero, i, MASTER_TO_SLAVE_TAG + 3, MPI_COMM_WORLD, &request);
            
            int initial_row=init->v[lower_bound].row;

            int final_row=init->v[upper_bound-1].row + 1;

            int portion_L=(final_row-initial_row)*init->nF; 


            slaves[i].initial_L=INDEX(initial_row,0,init->nF);
            slaves[i].final_L=INDEX(final_row,0,init->nF);
            slaves[i].num_elements_L=portion_L;

            MPI_Isend(&L1[INDEX(initial_row,0,init->nF)], portion_L, MPI_DOUBLE, i, MASTER_TO_SLAVE_TAG , MPI_COMM_WORLD, &request);
            MPI_Isend(&L2[INDEX(initial_row,0,init->nF)], portion_L, MPI_DOUBLE, i, MASTER_TO_SLAVE_TAG + 1 , MPI_COMM_WORLD, &request);


            MPI_Isend(&init->iter, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG  , MPI_COMM_WORLD, &request);
            MPI_Isend(&init->alpha, 1, MPI_DOUBLE, i, MASTER_TO_SLAVE_TAG + 2 , MPI_COMM_WORLD, &request); 
            MPI_Isend(&init->nU, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG + 3 , MPI_COMM_WORLD, &request);
            MPI_Isend(&init->nI, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG , MPI_COMM_WORLD, &request);
            MPI_Isend(&R[0], init->nF*init->nI , MPI_DOUBLE, i, MASTER_TO_SLAVE_TAG +2 , MPI_COMM_WORLD, &request);
            MPI_Isend(&R_hold[0], init->nF*init->nI, MPI_DOUBLE, i, MASTER_TO_SLAVE_TAG + 3 , MPI_COMM_WORLD, &request);

            lower_bound=upper_bound;
            curr_p--;
            
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

        printf("Iniciar slaves\n");  
        fflush(stdout);     
        MPI_Recv(&lower_bound, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&upper_bound, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG + 1, MPI_COMM_WORLD, &status);
        MPI_Recv(&num_Features, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG + 2, MPI_COMM_WORLD, &status);

        aux = malloc((upper_bound-lower_bound) * sizeof(non_zero));
        MPI_Recv(&aux[0], (upper_bound - lower_bound) , mpi_non_zero, 0, MASTER_TO_SLAVE_TAG + 3, MPI_COMM_WORLD, &status);
        int initial_row=aux[0].row;

        int final_row=aux[upper_bound - lower_bound-1].row + 1 ;

        int portion_L=(final_row-initial_row)*num_Features;


        L = MatrixInit((final_row-initial_row), num_Features);
        L_hold = MatrixInit((final_row-initial_row), num_Features);

        MPI_Recv(&L[0], portion_L , MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG , MPI_COMM_WORLD, &status);
        MPI_Recv(&L_hold[0], portion_L , MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG + 1, MPI_COMM_WORLD, &status);

        MPI_Recv(&iterations, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG  , MPI_COMM_WORLD, &status);
        MPI_Recv(&alpha_value, 1, MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG + 2 , MPI_COMM_WORLD, &status);
        MPI_Recv(&num_Users, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG + 3 , MPI_COMM_WORLD, &status);
        MPI_Recv(&num_Items, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG , MPI_COMM_WORLD, &status);


        R = MatrixInit(num_Features, num_Items);
        R_hold = MatrixInit( num_Features, num_Items);

        MPI_Recv(&R[0], num_Items*num_Features , MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG + 2, MPI_COMM_WORLD, &status);
        MPI_Recv(&R_hold[0], num_Items*num_Features , MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG + 3, MPI_COMM_WORLD, &status);
        MPI_Recv(&slaves[0], p , mpi_division_slaves, 0, MASTER_TO_SLAVE_TAG , MPI_COMM_WORLD, &status);

        L1 = L;
        L2 = L_hold; 
        R1 = R;
        R2 = R_hold;  
       
        matrix_mul(L1,R1,aux,upper_bound-lower_bound,num_Features);
     
        MPI_Isend(&lower_bound, 1, MPI_INT, 0, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD, &request);
        MPI_Isend(&upper_bound, 1, MPI_INT, 0, SLAVE_TO_MASTER_TAG + 1, MPI_COMM_WORLD, &request);
        MPI_Isend(&aux[0], (upper_bound - lower_bound), mpi_non_zero, 0, SLAVE_TO_MASTER_TAG + 2, MPI_COMM_WORLD, &request);
        fflush(stdout);   
             
               
    }
    
    if(id==0){
        matrix_mul(L1, R1, init->v, upper_bound-lower_bound , init->nF); 
        for(int i = 1 ; i<p ;i++){ // master process receives all results
            MPI_Recv(&lower_bound, 1, MPI_INT, i, SLAVE_TO_MASTER_TAG, MPI_COMM_WORLD, &status);
            MPI_Recv(&upper_bound, 1, MPI_INT, i, SLAVE_TO_MASTER_TAG + 1, MPI_COMM_WORLD, &status);
            MPI_Recv(&init->v[lower_bound], (upper_bound - lower_bound) , mpi_non_zero, i, SLAVE_TO_MASTER_TAG + 2, MPI_COMM_WORLD, &status);
        
        }   

        aux_sum_L= MatrixInit(init->nU, init->nF); // auxiliar matriz to sum matrix L
        aux_sum_R= MatrixInit(init->nF, init->nI); // auxiliar matriz to sum matrix L
        printf("Matrix L inteira\n");
        //printMatrix(L1,init->nU, init->nF);
    }

    int zeros_count;
    
        
    printf("process %d\n",id);
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
    fflush(stdout);
        
    
    int message_tag=0;

    for(int i = 0 ; i < iterations ; i++){
        //update the matrix
        if(i==0){
            MPI_Barrier(MPI_COMM_WORLD);
            break;
        }
        tmp = L1;
        L1 = L2;
        L2 = tmp;

        tmp = R1;
        R1 = R2;
        R2 = tmp; 
       
        message_tag=(i*2);
        zeros_count= slaves[id].upper_bound - slaves[id].lower_bound;

        if( id > 0){

            int initial_row=aux[0].row;
            int final_row=aux[zeros_count-1].row + 1 ;
            int user_portion= final_row-initial_row;  
            printf("user portion = %d com inicial row =%d final=%d id=%d\n",user_portion,initial_row,final_row,id);
            printMatrix(L1,user_portion,num_Features);
            recalculate_Matrix(L1, R1, L2, R2, user_portion, num_Items,num_Features,alpha_value, aux ,zeros_count ,id,p,slaves);
 
            
            MPI_Isend(&L1[0],slaves[id].num_elements_L, MPI_DOUBLE, 0, SLAVE_TO_MASTER_TAG + message_tag, MPI_COMM_WORLD, &request_send_0);
            MPI_Isend(&R1[0],num_Items*num_Features, MPI_DOUBLE, 0, SLAVE_TO_MASTER_TAG + message_tag +1 , MPI_COMM_WORLD, &request_send_1);

            MPI_Recv(&L1[0], slaves[id].num_elements_L , MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG +message_tag , MPI_COMM_WORLD, &status);
            MPI_Recv(&R1[0], num_Items*num_Features , MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG + message_tag+1, MPI_COMM_WORLD, &status);
            
     
      
            matrix_mul(L1,R1,aux,slaves[id].upper_bound-slaves[id].lower_bound,num_Features);


            MPI_Isend(&aux[0], (slaves[id].upper_bound-slaves[id].lower_bound), mpi_non_zero, 0, SLAVE_TO_MASTER_TAG + message_tag, MPI_COMM_WORLD, &request);
        }

        if(id==0){   
            recalculate_Matrix(L1, R1, L2, R2, num_Users, num_Items,num_Features,alpha_value, init->v ,zeros_count ,id,p,slaves);
            

            
            for(int k=1;k<p;k++){
                
                

                MPI_Recv(&aux_sum_L[slaves[k].initial_L],slaves[k].num_elements_L ,MPI_DOUBLE, k, SLAVE_TO_MASTER_TAG + message_tag, MPI_COMM_WORLD,&status);
                MPI_Recv(&aux_sum_R[0],num_Items*num_Features,MPI_DOUBLE, k, SLAVE_TO_MASTER_TAG + message_tag +1, MPI_COMM_WORLD,&status);
                

                for(int u = 0; u < num_Users; u++)
                    for(int f = 0; f < num_Features; f++){
                        L1[INDEX(u,f,num_Features)]+=aux_sum_L[INDEX(u,f,num_Features)];
                        /*if(L1[INDEX(u,f,num_Features)]<0)
                            printf("\nWRONG L\n");*/
                    }
                
                for(int c = 0; c < num_Items; c++)
                    for(int f = 0; f < num_Features; f++){
                        R1[INDEX(c,f,num_Features)] += aux_sum_R[INDEX(c,f,num_Features)];  
                        /*if(R1[INDEX(c,f,num_Features)]<0)
                            printf("\nWRONG R\n");*/
                    }    
                
                zero_LR(aux_sum_L, aux_sum_R, init->nU, init->nI, init->nF);
               


            }

            for(int u = 0; u < num_Users; u++)
                for(int f = 0; f < num_Features; f++){
                    L1[INDEX(u,f,num_Features)] = L2[INDEX(u,f,num_Features)] + alpha_value*2*L1[INDEX(u,f,num_Features)];
                    //printf("L1 value= %f row: %d column: %d\n",L1[INDEX(u,f,num_Features)],u,f);
                    /*if(L1[INDEX(u,f,num_Features)]<0)
                        printf("\nWRONG no recalc L\n");*/
                }

            for(int k = 0; k < num_Items; k++)
                for(int f = 0; f < num_Features; f++){
                    R1[INDEX(k,f,num_Features)] = R2[INDEX(k,f,num_Features)] + alpha_value*2*R1[INDEX(k,f,num_Features)]; 
                    /*if(R1[INDEX(k,f,num_Features)]<0)
                        printf("\nWRONG no recal R\n");*/
                }  

            for(int k=1;k<p;k++){
                MPI_Isend(&L1[slaves[k].initial_L], slaves[k].num_elements_L, MPI_DOUBLE, k, MASTER_TO_SLAVE_TAG + message_tag , MPI_COMM_WORLD, &request);
                MPI_Isend(&R1[0], init->nF*init->nI , MPI_DOUBLE, k, MASTER_TO_SLAVE_TAG +message_tag+1 , MPI_COMM_WORLD, &request);
               
            }
            matrix_mul(L1, R1, init->v, slaves[id].upper_bound-slaves[id].lower_bound, init->nF);
        
            for(int k = 1 ; k<p ;k++){ // master process receives all results
                MPI_Recv(&init->v[slaves[k].lower_bound], (slaves[k].upper_bound-slaves[k].lower_bound) , mpi_non_zero, k, SLAVE_TO_MASTER_TAG +message_tag, MPI_COMM_WORLD, &status);
        
            }   
            printf("Matrix L na iteracao %d\n",i);
            printMatrix(L1,init->nU, init->nF); 
            printf("Matrix R na iteracao %d\n",i);       
            printMatrix(R1,init->nI, init->nF); 
        
        
        }




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