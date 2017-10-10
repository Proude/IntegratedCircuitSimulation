#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "hash_table.h"
#include "csparse.h"

void Cholesky_init(gsl_matrix_view *A, gsl_vector_view *b, double *MNA_matrix, double *MNA_right_vector)
{
  *A = gsl_matrix_view_array(MNA_matrix, mnac+hash_id, mnac+hash_id);
  *b = gsl_vector_view_array(MNA_right_vector, mnac+hash_id);
}

void CholeskyDecomposition(gsl_matrix_view A, gsl_vector_view b, gsl_vector *x){
  
  gsl_linalg_cholesky_decomp(&A.matrix);
  gsl_linalg_cholesky_solve(&A.matrix, &b.vector, x);
}

void CholeskySparseDecomposition(cs *C, double *x, double *b)
{
  css *S = NULL;
  csn *N = NULL;
  int i;
  
  S = cs_schol(1, C);
  N = cs_chol(C, S);
  
  if (!cs_ipvec(S->pinv, b, x, mnac+hash_id))
  {
    printf("An error occured during cholesky decomposition. Program exits!\n");
    exit(0);
  }
  if (!cs_lsolve(N->L, x))
  {
    printf("An error occured during cholesky decomposition. Program exits!\n");
    exit(0);
  }
  if (!cs_ltsolve(N->L, x))
  {
    printf("An error occured during cholesky decomposition. Program exits!\n");
    exit(0);
  }
  if (!cs_pvec(S->pinv, x, b, mnac+hash_id))
  {
    printf("An error occured during cholesky decomposition. Program exits!\n");
    exit(0);
  }
  
  for (i = 0; i < mnac + hash_id; i++)
    printf("%lf\n", x[i]);
  
  cs_sfree(S);
  cs_nfree(N);
}