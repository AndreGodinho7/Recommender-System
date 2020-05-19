#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"input.h"
#include "factorization.h"
#include <limits.h>
#include <mpi.h>
#include <omp.h>

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
    
    //MPI_Init(&argc, &argv);

    int provided;
    MPI_Init_thread(&argc,&argv, MPI_THREAD_FUNNELED, &provided);
    
    int id, p;

	MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

    omp_set_num_threads(4);
    /* create a type for struct typedef _division_mpi */
    const int nitemsV2=6;
    int          blocklengthsV2[6] = {1,1,1,1,1,1};
    MPI_Datatype typesV2[6] = {MPI_INT, MPI_INT,MPI_INT,MPI_INT,MPI_INT,MPI_INT};
    MPI_Datatype mpi_division_slaves;
    MPI_Aint     offsetsV2[6];

    offsetsV2[0] = offsetof(division_mpi, lower_bound);
    offsetsV2[1] = offsetof(division_mpi, upper_bound);
    offsetsV2[2] = offsetof(division_mpi, initial_L);
    offsetsV2[3] = offsetof(division_mpi, num_elements_L);
    offsetsV2[4] = offsetof(division_mpi, initial_row);
    offsetsV2[5] = offsetof(division_mpi, final_row);

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

    double* aux_sum_R;


    int* result_per_process;
    int* all_results;

    // Allocacao de todas as variaveis para a master
    input_values* init ;
    division_mpi* slaves; //struct that saves position of L and elements that each processor gets;
    
    slaves = malloc(p * sizeof(division_mpi));
    
    int num_Features;
    int num_Items;
    int num_Users;
    double alpha_value;
    int iterations;
    int initial_row;
    int final_row;


    if(id==0){
        
        if (argc != 2){
                printf("ERROR: inserted more than 1 input file.\n");
                exit(0);
        };        
        // for allocating matrice
        init = read_input(argv[1]);
        mark_process_in_nonzero(init->num_zeros, init->v, p);
        
        non_zero* aux;


        //L_sent = MatrixInit(init->nU, init->nF);
        R = MatrixInit(init->nF, init->nI);
        
        //aux_sum_L = L_sent; 

        R_hold = MatrixInit(init->nF, init->nI); // Matrix that stores the previous iteration of R

        //random_fill_LR(aux_sum_L, R, init->nU, init->nI, init->nF);



        //==============================================================
        //divisao das matrizes para os slaves
        
       
       
        all_results=(int*)malloc(init->nU*sizeof(int));

        int my_up;
        int my_lower;

     
        my_lower=0;
        my_up=getProcessUpBoundary(init->v,init->num_zeros,id);

        initial_row=init->v[my_lower].row;
        
        final_row=init->v[my_up-1].row + 1 ;//+1
        int portion_L=(final_row-initial_row)*init->nF;
        int user_rows=final_row-initial_row;
        
        result_per_process=(int*)malloc(user_rows*sizeof(int));
        //printf("num users=%d no process =%d\n",user_rows,id);
        //fflush(stdout);

        L_hold = MatrixInit(user_rows, init->nF);
        L=MatrixInit(user_rows,init->nF);

        random_fill_LR(L,R,init->nU,init->nI,init->nF,initial_row,final_row);
        
        

        R = transpose(R, init->nF, init->nI);  
        R_hold = transpose(R_hold, init->nF, init->nI);         
       
        R1 = R;
        R2 = R_hold;
        L2 = L_hold;
        L1=L;


        /*printf("Estou a ver da linha %d ate %d no processo %d\n",initial_row,final_row,id);
        fflush(stdout);*/
        
        slaves[0].lower_bound=my_lower;
        slaves[0].upper_bound=my_up;        
        slaves[0].num_elements_L=portion_L;
        slaves[0].initial_L=INDEX(initial_row,0,init->nF);
        slaves[0].initial_row=initial_row;
        slaves[0].final_row=final_row;   
        //printf("Existem %d processos\n",p);


        /*printf("MASTER : processo %d gets from %d to %d\n",id,my_lower,my_up);
        fflush(stdout);*/
        
        lower_bound=my_up;
        //int message_tag
        for(int i=1;i<p;i++){
 
            
            
            upper_bound=getProcessUpBoundary(init->v,init->num_zeros,i);

            slaves[i].lower_bound=lower_bound;
            slaves[i].upper_bound=upper_bound;

            /*printf("MASTER : processo %d gets from %d to %d\n",i,slaves[i].lower_bound,slaves[i].upper_bound);
            fflush(stdout);*/

            MPI_Isend(&slaves[i].lower_bound, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG +i, MPI_COMM_WORLD, &request);
            MPI_Isend(&slaves[i].upper_bound, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG + 1 +i, MPI_COMM_WORLD, &request);
            MPI_Isend(&init->nF, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG + 2 +i, MPI_COMM_WORLD, &request);
            MPI_Isend(&init->v[slaves[i].lower_bound], (slaves[i].upper_bound - slaves[i].lower_bound), mpi_non_zero, i, MASTER_TO_SLAVE_TAG + 3 +i, MPI_COMM_WORLD, &request);
            


            if(init->v[slaves[i].lower_bound].row != init->v[slaves[i-1].upper_bound-1].row+1){
                initial_row=init->v[slaves[i-1].upper_bound-1].row+1;
            }
            else{
                initial_row=init->v[slaves[i].lower_bound].row;
            }
            final_row=init->v[slaves[i].upper_bound-1].row + 1;

            slaves[i].initial_row=initial_row;
            slaves[i].final_row=final_row;
            int portion_L=(slaves[i].final_row-slaves[i].initial_row)*init->nF; 
            /*printf("VOU enviar da linha %d ate %d no processo %d\n",slaves[i].initial_row,slaves[i].final_row,i);
            fflush(stdout);*/

            MPI_Isend(&slaves[i].initial_row, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG +3+i, MPI_COMM_WORLD, &request);
            MPI_Isend(&slaves[i].final_row, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG + 4 +i, MPI_COMM_WORLD, &request);
            

            slaves[i].initial_L=INDEX(initial_row,0,init->nF);
            slaves[i].num_elements_L=portion_L;


            //MPI_Isend(&aux_sum_L[INDEX(initial_row,0,init->nF)], portion_L, MPI_DOUBLE, i, MASTER_TO_SLAVE_TAG , MPI_COMM_WORLD, &request);
            //MPI_Isend(&aux_sum_L[INDEX(initial_row,0,init->nF)], portion_L, MPI_DOUBLE, i, MASTER_TO_SLAVE_TAG + 1 , MPI_COMM_WORLD, &request);
            
            

            MPI_Isend(&init->iter, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG  , MPI_COMM_WORLD, &request);
            MPI_Isend(&init->alpha, 1, MPI_DOUBLE, i, MASTER_TO_SLAVE_TAG + 2 , MPI_COMM_WORLD, &request); 
            MPI_Isend(&init->nU, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG + 3 , MPI_COMM_WORLD, &request);
            MPI_Isend(&init->nI, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG , MPI_COMM_WORLD, &request);

                    
            lower_bound=slaves[i].upper_bound;
            
            
        }

        int OK;
        for(int i = 1 ; i<p ;i++){ // master process receives all results
            MPI_Recv(&OK, 1 , MPI_INT, i, SLAVE_TO_MASTER_TAG+5, MPI_COMM_WORLD, &status);
        
        }  

        lower_bound=my_lower;
        upper_bound=my_up;  
        int portion_zero=upper_bound-lower_bound;



        aux =(non_zero*) malloc( portion_zero* sizeof(non_zero));
        for(int y=0;y<portion_zero;y++){
            aux[y]=init->v[y];
        }

        free(init->v);
        init->v =(non_zero*) malloc( portion_zero* sizeof(non_zero));
        for(int y=0;y<portion_zero;y++){
            init->v[y]=aux[y];
        } 
        free(aux);       

        for(int i = 1 ; i<p ;i++){ 
            MPI_Isend(&slaves[0], p, mpi_division_slaves, i, MASTER_TO_SLAVE_TAG , MPI_COMM_WORLD, &request);
            
        }  
        num_Users=init->nU;
        num_Items=init->nI;
        num_Features=init->nF;
        alpha_value=init->alpha;
        iterations=init->iter;

        /*printf("ENVIEI TUDO DA MASTER \n");
        fflush(stdout);*/
         
    }
    if(id>0){

        MPI_Recv(&lower_bound, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG +id, MPI_COMM_WORLD, &status);
        MPI_Recv(&upper_bound, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG + 1+id, MPI_COMM_WORLD, &status);
        MPI_Recv(&num_Features, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG + 2 +id, MPI_COMM_WORLD, &status);
      
        init = (input_values*)malloc(sizeof(input_values));
        init->v = (non_zero*)malloc((upper_bound-lower_bound) * sizeof(non_zero));

        
        MPI_Recv(&init->v[0], (upper_bound - lower_bound) , mpi_non_zero, 0, MASTER_TO_SLAVE_TAG + 3+id, MPI_COMM_WORLD, &status);
 
        MPI_Recv(&initial_row, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG +3+id, MPI_COMM_WORLD, &status);
        MPI_Recv(&final_row, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG + 4+id, MPI_COMM_WORLD, &status);
        
        /*printf("Recebi da linha %d ate %d no processo %d\n",initial_row,final_row,id);
        fflush(stdout);*/

        


        int init_row=final_row-initial_row;

        result_per_process=(int*)malloc(init_row*sizeof(int));

        L = MatrixInit(init_row, num_Features);
      
        L_hold = MatrixInit(init_row, num_Features);




        MPI_Recv(&iterations, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG  , MPI_COMM_WORLD, &status);
        MPI_Recv(&alpha_value, 1, MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG + 2 , MPI_COMM_WORLD, &status);
        MPI_Recv(&num_Users, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG + 3 , MPI_COMM_WORLD, &status);
        MPI_Recv(&num_Items, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG , MPI_COMM_WORLD, &status);

        int OK=0;
        MPI_Isend(&OK, 1, MPI_INT, 0, SLAVE_TO_MASTER_TAG+5  , MPI_COMM_WORLD, &request);

        R = MatrixInit( num_Features,num_Items);
        R_hold = MatrixInit(  num_Features,num_Items);

        

        
       
        MPI_Recv(&slaves[0], p , mpi_division_slaves, 0, MASTER_TO_SLAVE_TAG , MPI_COMM_WORLD, &status);


        random_fill_LR(L,R,num_Users,num_Items,num_Features,slaves[id].initial_row,slaves[id].final_row);

        R = transpose(R,num_Features, num_Items);  
        R_hold = transpose(R_hold, num_Features, num_Items); 
         

        L1 = L;
        L2 = L_hold; 
        R1 = R;
        R2 = R_hold;  

        /*printf(" antes do mul Matrix L (process %d) \n",id);
        printMatrix(L1,init_row, num_Features); 
        printf("antes do mul Matrix R no process %d\n",id);       
        printMatrix(R1,num_Items, num_Features); */
        //printf("antes do mul ;start = %d para processo =%d\n",)

        matrix_mul(L1,R1,init->v,upper_bound-lower_bound,num_Features,slaves[id].initial_row,init_row);
    

        

           
             
               
    }
    
    if(id==0){
        matrix_mul(L1, R1, init->v, upper_bound-lower_bound , init->nF,slaves[id].initial_row,slaves[id].final_row-slaves[id].initial_row); 

  
    }

    aux_sum_R= MatrixInit(num_Items, num_Features);

    
    int *displs = (int *)malloc(p*sizeof(int)); 
    int *scounts = (int *)malloc(p*sizeof(int));

    /*printf("Chegeui á parallel section\n");
    fflush(stdout);*/

    #pragma omp parallel 
    {
        int l;
        #pragma omp for private(l) 
        for ( l=0; l<p; l++) { 
            displs[l] = slaves[l].initial_L; 
            scounts[l] = slaves[l].num_elements_L; 
        
        } 

        /*printf("fiz o displys and counts\n");
        fflush(stdout);*/


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
            
        
        
        int user_portion;
        user_portion= slaves[id].final_row-slaves[id].initial_row;  
        /*printf("iteraçoes = %d\n",iterations);
        fflush(stdout);
        printf("Matrix L (process %d) inicial\n",id);
        printMatrix(L1,user_portion, num_Features); 
        printf("Matrix R \n");       
        printMatrix(R1,num_Items, num_Features);*/

        for(int i = 0 ; i < iterations ; i++){ 
            //printf("For iteracao = %d\n",i);
            fflush(stdout);
            //update the matrix
            /*if(i==5){
                MPI_Barrier(MPI_COMM_WORLD);
                break;
            }*/
            #pragma omp single
            {        
                tmp = L1;
                L1 = L2;
                L2 = tmp;

                tmp = R1;
                R1 = R2;
                R2 = tmp; 
            }
            //printf("mudei vectores\n");
            //fflush(stdout);
            zeros_count= slaves[id].upper_bound - slaves[id].lower_bound;
            //printf("zeros =%d\n",zeros_count);
            //fflush(stdout);


            recalculate_Matrix(L1, R1, L2, R2, user_portion, num_Items,num_Features,alpha_value, init->v ,zeros_count ,id,p,slaves);
            //printf("fiz o recalculate\n");
            //fflush(stdout);
            #pragma omp barrier
            #pragma omp master
            {
                MPI_Allreduce(R1, aux_sum_R, num_Items*num_Features, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            }
            #pragma omp barrier
            
            int k,f;
            #pragma omp for private(f,k) 
            for( k = 0; k < num_Items; k++){
                for( f = 0; f < num_Features; f++){
                    R1[INDEX(k,f,num_Features)] = R2[INDEX(k,f,num_Features)] + alpha_value*2*aux_sum_R[INDEX(k,f,num_Features)];                 
                } 
            }


            matrix_mul(L1,R1,init->v,slaves[id].upper_bound-slaves[id].lower_bound,num_Features,slaves[id].initial_row,user_portion);

        
            /*printf("Matrix L (process %d) na iteracao %d\n",id,i);
            printMatrix(L1,user_portion, num_Features); 
            printf("Matrix R na iteracao %d\n",i);       
            printMatrix(R1,num_Items, num_Features); */
            


        }
    }
    //printf("cheguei ao final\n");
    //fflush(stdout);

    int user_portion;
    user_portion= slaves[id].final_row-slaves[id].initial_row;   

    /*printf("Matrix L (process %d) \n",id);
    printMatrix(L1,user_portion, num_Features); 
    printf("Matrix R \n");       
    printMatrix(R1,num_Items, num_Features); */
    free(R2);
    free(L2);
    free(aux_sum_R);

    create_output(init->v,user_portion , num_Items, num_Features, L1, R1, slaves[id].upper_bound-slaves[id].lower_bound,result_per_process,slaves[id].initial_row);

    for (int i=0; i<p; i++) { 
        displs[i] = slaves[i].initial_row; 
        scounts[i] = slaves[i].final_row-slaves[i].initial_row; 
       
    } 
    MPI_Gatherv(&result_per_process[0], user_portion, MPI_INT, &all_results[0], scounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

    if(id==0){
        for(int i=0;i<num_Users;i++){
            printf("%d\n",all_results[i]);
        }
        free(all_results);
    }

        
    
    free(result_per_process);
    free(init->v);
    free(init);    
    free(L1);
    free(slaves);
    free(R1);
    free(displs);
    free(scounts);
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize ();
    return 0;
}