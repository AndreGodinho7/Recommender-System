
#ifndef input
#define input


struct _non_zero{
  int row;
  int column;
  double value;
};
struct _internal_product{
  int row;
  int column;
  double value;
};


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

typedef struct _internal_product internal_product;
typedef struct _input_values input_values;
typedef struct _non_zero _non_zero;

input_values* read_input(char* filename);

#endif
