#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "input.h"

struct _input_values{
  int iter;
  double alpha;
  int features;
  int rows;
  int columns;
  int non_zeros;
  double** matrix;
};

input_values* read_input(char* filename)
{
    FILE *fpIn;
    input_values* init ;
    int r = 0, c = 0;
    double ele = 0;
  
    init = (input_values*)malloc(sizeof(input_values));
    if (init == NULL){
      printf("ERROR: Out of memory.\n");
		  exit(1);
    }
    
    fpIn = fopen(filename, "r");
    if (fpIn == NULL){
      printf("ERROR: input file not open.");
    }

    fscanf(fpIn, "%d", &(init->iter) ); // 1st line
    fscanf(fpIn, "%lf", &(init->alpha)); // 2nd line
    fscanf(fpIn, "%d", &(init->features)); // 3rd line
    fscanf(fpIn, "%d %d %d", &init->rows, &init->columns, &init->non_zeros); // 4th line

    init->matrix = (double **)malloc(init->rows * sizeof(double *)); 
    if (init->matrix == NULL){
      printf("ERROR: Out of memory.\n");
		  exit(1);
    }
    for (int i = 0 ; i < init->rows ; i++) {
         init->matrix[i] = (double *)malloc(init->columns * sizeof(double));
    }
    
    while (fscanf(fpIn,"%d %d %lf", &r, &c, &ele) != EOF){ // remaining lines
        init->matrix[r][c] = ele;
    }

    return init;
}