#ifndef factorization
#define factorization

void printMatrix(double* matrix, int rows, int columns);

int getProcessUpBoundary(non_zero* v, int num_zeros,int p);

void copyPartMatrix(double* from, double* to,int rows, int columns);

void mark_process_in_nonzero(int num_zeros, non_zero *v, int NUM_PROCESSES);

double* MatrixInit(int rows, int columns);

void random_fill_LR(double* L_, double* R_, int nU, int nI, int nF,int initial_row,int final_row);

double* transpose(double* matrix, int rows, int columns);

void matrix_mul(double *firstMatrix, double *secondMatrix, non_zero* v,int num_zeros, int nF,int start,int num_users);

void zero_LR(double* L, double* R, int nU, int nI, int nF);

//void recalculate_Matrix(double** L, double** R,double** pre_L, double** pre_R,double **A,double** B, double** pre_B,double** L_calc,double** R_calc,int nU, int nI, int nF,int iter, double alpha, _non_zero *v, int non_zero);
void recalculate_Matrix(double* L, double* R,double* pre_L, double* pre_R, int nU, int nI, int nF, double alpha, non_zero *v, int non_zero,int id,int p, division_mpi *slaves);

void copy_matrix(double* original, double* copied,int rows, int columns);

void create_output(non_zero *v, int nU, int nI, int nF, double* L, double* R, int num_zeros,int *result,int start);

int find_upper_bound(int lower_row,int upper_row,int lower_bound,non_zero *v,int num_zeros);



#endif