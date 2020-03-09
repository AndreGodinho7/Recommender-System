#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "input.h"
#include "factorization.h"


/******************************************************************************
* read_input()
*
* Arguments: filename - character pointer to input filename
*
* Returns: pointer to an element of the type input_values
*
* Side-Effects: if any allocation fails exits with exit(1)
*               if the input does not open exits with exit(1)
*
* Description: reads input information of the instance to run
*
*****************************************************************************/

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
      exit(1);
    }

    fscanf(fpIn, "%d", &(init->iter) ); // 1st line
    fscanf(fpIn, "%lf", &(init->alpha)); // 2nd line
    fscanf(fpIn, "%d", &(init->nF)); // 3rd line
    fscanf(fpIn, "%d %d %d", &init->nU, &init->nI, &init->non_zeros); // 4th line

    init->matrix = MatrixInit(init->nU, init->nI);
    
    while (fscanf(fpIn,"%d %d %lf", &r, &c, &ele) != EOF){ // remaining lines
        init->matrix[r][c] = ele;
    }

    return init;
}