
#ifndef input
#define input


struct _non_zero{
  int process;
  int row;
  int column;
  double A;
  double B;
  int process;
};

typedef struct _non_zero non_zero;

struct _division_mpi{
  int lower_bound;
  int upper_bound;
  int initial_L;
  int final_L;
  int num_elements_L;
};

typedef struct _division_mpi division_mpi;


struct _input_values{
  int iter;
  double alpha;
  int nF;
  int nU;
  int nI;
  int num_zeros;
  non_zero* v;
  double* matrix;
};

typedef struct _input_values input_values;

input_values* read_input(char* filename);

#endif
