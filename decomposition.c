#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "hash_table.h"

void init(gsl_matrix_view *A, gsl_vector_view *b, double *MNA_matrix, double *MNA_right_vector)
{
  *A = gsl_matrix_view_array(MNA_matrix, mnac+hash_id, mnac+hash_id);
  *b = gsl_vector_view_array(MNA_right_vector, mnac+hash_id);
}

void LUDecomposition(gsl_matrix_view A, gsl_vector_view b, gsl_vector *x)
{ 
  int s;
  gsl_permutation *p = gsl_permutation_alloc(mnac+hash_id);

  gsl_linalg_LU_decomp (&A.matrix, p, &s);
  gsl_linalg_LU_solve (&A.matrix, p, &b.vector, x);
  gsl_permutation_free(p);
}

void LUSparseDecomposition(cs *C, double *x, double *b)
{
  int i;
  css *S = NULL;
  csn *N = NULL;
  
  S = cs_sqr(2, C, 0);
  N = cs_lu(C, S, 1);
  
  if (!cs_ipvec(N->pinv, b, x, mnac + hash_id))
  {
    printf("An error occured during LU decomposition. Program exits!\n");
    exit(0);
  }
  if (!cs_lsolve(N->L, x))
  {
    printf("An error occured during LU decomposition. Program exits!\n");
    exit(0);
  }
  if (!cs_usolve(N->U, x))
  {
    printf("An error occured during LU decomposition. Program exits!\n");
    exit(0);
  }
  if (!cs_ipvec(S->q, x, b, mnac + hash_id))
  {
    printf("An error occured during LU decomposition. Program exits!\n");
    exit(0);
  }
  
  for (i = 0; i < mnac + hash_id; i++)
    printf("%lf\n", x[i]);
  
  cs_sfree(S);
  cs_nfree(N);
}