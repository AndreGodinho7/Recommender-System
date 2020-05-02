#ifndef factorization
#define factorization

void printMatrix(double** matrix, int rows, int columns);

double** MatrixInit(int rows, int columns);

void random_fill_LR(double** L_, double** R_, int nU, int nI, int nF);

double** transpose(double** matrix, int rows, int columns);

void matrix_mul(double **firstMatrix, double **secondMatrix, non_zero* v,int num_zeros, int nF);
void matrix_mul_mpi(double **firstMatrix, double **secondMatrix, non_zero* v,int lower, int up, int nF);

void zero_LR(double** L, double** R, int nU, int nI, int nF);

//void recalculate_Matrix(double*** L, double*** R,double*** pre_L, double*** pre_R,double ***A,double*** B, double*** pre_B,double*** L_calc,double*** R_calc,int nU, int nI, int nF,int iter, double alpha, _non_zero *v, int non_zero);
void recalculate_Matrix(double** L, double** R,double** pre_L, double** pre_R, int nU, int nI, int nF, double alpha, non_zero *v, int non_zero);

void copy_matrix(double** original, double** copied,int rows, int columns);

void create_output(non_zero *v, int nU, int nI, int nF, double** L, double** R, int num_zeros);


#endif