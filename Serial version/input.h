
#ifndef input
#define input


struct _non_zero{
  int row;
  int column;
  double A;
  double B;
};

typedef struct _non_zero non_zero;

struct _input_values{
  int iter;
  double alpha;
  int nF;
  int nU;
  int nI;
  int num_zeros;
  non_zero* v;
  double** matrix;
};

typedef struct _input_values input_values;

input_values* read_input(char* filename);

#endif
