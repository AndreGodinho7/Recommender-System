#ifndef input
#define input

struct _input_values{
  int iter;
  double alpha;
  int nF;
  int nU;
  int nI;
  int non_zeros;
  double** matrix;
};

typedef struct _input_values input_values;

input_values* read_input(char* filename);

#endif
