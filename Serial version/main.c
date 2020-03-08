#include <stdio.h>
#include <stdlib.h>

#include"input.h"

int main(int argc, char* argv[])
{
    input_values* init;
 
    if (argc != 2){
        printf("ERROR: inserted more than 1 input file.\n");
        exit(0);
    };

    init = read_input(argv[1]);
    


    return 0;
}