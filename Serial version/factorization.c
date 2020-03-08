#include <stdio.h>
#include <stdlib.h>

#define RAND01 ((double) random() / (double) RAND_MAX)

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

void random_fill_LR(double*** L, double*** R, int nU, int nI, int nF)
{   
    srandom(0);
    for(int i = 0; i < nU; i++)
        for(int j = 0; j < nF; j++)
            (*L)[i][j] = RAND01 / (double) nF;

    for(int i = 0; i < nF; i++)
        for(int j = 0; j < nI; j++)
            (*R)[i][j] = RAND01 / (double) nF;
}

void matrix_mul(double ***firstMatrix, double ***secondMatrix, double ***matrix3,int nU, int nI, int nF ){


    int rowFirst= nU;//sizeof(matrix_1)/sizeof(matrix_1[0]);
    int columnFirst =nF;// sizeof(matrix_1[0])/sizeof(matrix_1[0][0]);
    
    int columnSecond = nI;//sizeof(matrix_2[0])/sizeof(matrix_2[0][0]);
        
            
        
    double result=0.0;
	for(int i = 0; i < rowFirst; ++i)
	{
		for(int j = 0; j < columnSecond; ++j)
		{
			for(int k=0; k<columnFirst; ++k)
			{
               
                result=((*firstMatrix)[i][k] * (*secondMatrix)[k][j]);
				(*matrix3)[i][j] = (*matrix3)[i][j]+ result ;
                                
                                
			}
                        
		}
	}     
}


void recalculate_Matrix(double*** L, double*** R,double*** pre_L, double*** pre_R,double ***A,double*** B, double*** pre_B,int nU, int nI, int nF,int iter, double alpha){
    for(int i=0;i<nU;i++){
        for(int j=0; j<nI;j++){
            if((*A)[i][j]!=0.0){
                //printf("%f \n",(*A)[i][j]);
                for(int feature=0; feature < nF; feature ++){
                    double sum_L=0;
                    double sum_R=0;
                                        
                    for(int col=0;col<nI;col++){
                        //printf("im here\n");
                        sum_L=sum_L+(2*((*A)[i][col]-(*pre_B)[i][col])*(-(*pre_R)[feature][col]));
                                                
                    }
                    for(int line=0;line<nU;line++){
                        sum_R=sum_R+(2*((*A)[line][j]-(*pre_B)[line][j])*(-(*pre_L)[line][feature]));
                    }
                    (*L)[i][feature]=(*pre_L)[i][feature]-alpha*sum_L;
                    (*R)[feature][j]=(*pre_R)[feature][j]-alpha*sum_R;
                }
            }
        }
                //printf("\n");
    }   
}