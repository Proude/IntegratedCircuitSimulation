#include <gsl/gsl_matrix.h>
#include "csparse.h"

struct entry_s {
	char *key;
	size_t value;
	struct entry_s *next;
};

typedef struct entry_s entry_t;

struct hashtable_s {
	int size;
	struct entry_s **table;	
};

typedef struct hashtable_s hashtable_t; 
 
typedef struct elements
{
  char *name;
  char type;
  char *plus_d_c;
  char *minus_g_b;
  char *source_e;
  char *body;
  char *model_name;
  double L;
  double W;
  double area;
  double value;
  char *trans_type;
  double i1;
  double i2;
  double td1;
  double tc1;
  double td2;
  double tc2;
  double per;
  double *pwl;
  struct elements *nextPtr;
} elements;

// Function declarations
void add_to_list(char *name, char type, char *plus_d_c, char *minus_g_b, char *source_e, char *body, double value, char *model_name, double L,
		 double W, double area, char *trans_type, double i1, double i2, double td1, double tc1, double td2, double tc2, double per, double *pwl, int size_of_pwl);
void print_list(hashtable_t *hashtable);
void list_init();
int read_file(char *file_name, hashtable_t *hashtable, int *plot_list, int *checkForCholesky, char *typeOfDC, double *startValue2, double *endValue2, double *stepValue2, int *checkForDC, char *nameOfDC,int *checkForIter,int *checkForSpdIter,double *checkForItol,int *CheckForSparse,int *CheckForSpdSparce,int *CheckForIterSparce, int *CheckForSpdIterSparce, double *time_step, double *fin_time, int *numpwl,int *CheckForMethodTr,int *CheckForMethodBe);
hashtable_t *ht_create( int size );
int ht_hash( hashtable_t *hashtable, char *key );
entry_t *ht_newpair( char *key, size_t value );
void ht_set( hashtable_t *hashtable, char *key, size_t value );
size_t ht_get( hashtable_t *hashtable, char *key );
void create_MNA_matrix(hashtable_t *hashtable, double *MNA_matrix, double *MNA_right_vector);
void LUDecomposition(gsl_matrix_view A, gsl_vector_view b, gsl_vector *x);
void init(gsl_matrix_view *A, gsl_vector_view *b, double *MNA_matrix, double *MNA_right_vector);
void DC_init(gsl_matrix_view *A, gsl_vector_view *b, double *MNA_matrix, double *MNA_right_vector,char type, char var_name[], double val, hashtable_t *hashtable);
void DC_simulation(gsl_matrix_view A, gsl_vector_view b, gsl_vector *x);
void Cholesky_init(gsl_matrix_view *A, gsl_vector_view *b, double *MNA_matrix, double *MNA_right_vector);
void CholeskyDecomposition(gsl_matrix_view A, gsl_vector_view b, gsl_vector *x);
void BiCG(double *A, double *b, double *x, double itol);
double normC(double *mat);
void solveForBiCG(double *M, double *z, double *r);
void calculateQ(double *q, double *A, double *p);
void CG_Method(double *MNA_Matrix, double *b, double *x, double itol, int x_length);
void Calculate_Init_r(double *MNA_Matrix,double *MNA_right_vector,double *x,double *r,int r_length);
double norm2(double *vector,int length);
void Calculate_z(double *M,double *r,double *z,int length);
double dot_Product(double *x,double *y, int length);
void Multiply_Add(double *x, double *y, double *res, double val, int length); //res(vector) = x(vector) + val * y(vector)
void Multiply_Reduce(double *x, double *y, double *res, double val, int length); //res(vector) = x(vector) - val * y(vector)
void Product_Ap(double *A, double *p, double *q, int length);
int Calculate_nonZeros(hashtable_t *hashtable);
void create_Sparse_matrix(hashtable_t *hashtable,cs *Sparse_MNA_matrix,double *MNA_right_vector, int nonZeros);
void CholeskySparseDecomposition(cs *C, double *x, double *b);
void LUSparseDecomposition(cs *C, double *x, double *b);
void CG_Method_Sparse(cs *Sparse_MNA_matrix, double *right_vector, double *x, double itol, int x_length);
void BiCG_Sparse(cs *A, double *b, double *x, double itol);
void removeChar(char *str, char garbage);
void init_starting_value_LU(double *MNA_matrix, double *MNA_right_vector, double *xLU, double step_time, double fin_time);
void create_G_matrix(hashtable_t *hashtable, double *MNA_matrix, double *C_matrix,double *trans_right_vector, double step_time, double fin_time, int h_val, int numpwl);
void init_starting_value_Cholesky(double *MNA_matrix, double *MNA_right_vector, double *xLU, double step_time, double fin_time);
void trapezoidal_LU(double *G_matrix, double *C_matrix, double *x_final, int h_val, double *trans_right_matrix);
void trapezoidal_Cholesky(double *G_matrix, double *C_matrix, double *x_final, int h_val, double *trans_right_matrix);
void trapezoidal_CG(double *G_matrix, double *C_matrix, double *x_final, int h_val, double *trans_right_matrix, double itol);
void trapezoidal_BiCG(double *G_matrix, double *C_matrix, double *x_final, int h_val, double *trans_right_matrix, double itol);
void BE_method(double *BE,double *right_vector, double *G, double *C, double *xtkprev,double *etk, int h_val);

size_t hash_id;
size_t mnac;
elements *root;
int plot_count_words;