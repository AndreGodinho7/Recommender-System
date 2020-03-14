#ifndef factorization
#define factorization

void printMatrix(double** matrix, int rows, int columns);

double** MatrixInit(int rows, int columns);

void random_fill_LR(double** L, double** R, int nU, int nI, int nF);

void matrix_mul(double **firstMatrix, double **secondMatrix,double **mult_matrix,int nU, int nI, int nF);

//void recalculate_Matrix(double*** L, double*** R,double*** pre_L, double*** pre_R,double ***A,double*** B, double*** pre_B,double*** L_calc,double*** R_calc,int nU, int nI, int nF,int iter, double alpha, _non_zero *v, int non_zero);
void recalculate_Matrix(double** L, double** R,double** pre_L, double** pre_R,double **A,double** B, double** pre_B, int nU, int nI, int nF,int iter, double alpha, _non_zero *v, int non_zero);


#endif