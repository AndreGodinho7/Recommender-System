
#ifndef input
#define input


typedef struct _non_zero{
  int row;
  int column;
  double value;
}_non_zero;


struct _input_values{
  int iter;
  double alpha;
  int nF;
  int nU;
  int nI;
  int non_zeros;
  double** matrix;
  _non_zero *v;
};


typedef struct _input_values input_values;

input_values* read_input(char* filename);

#endif
