#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "hash_table.h"
#include "csparse.h"

void init_starting_value_LU(double *MNA_matrix, double *MNA_right_vector, double *xLU, double step_time, double fin_time)
{
  int i;
  size_t iiter = mnac + hash_id;
  size_t jiter = (fin_time / step_time) + 1;
  double b1[iiter];
  gsl_matrix_view A;
  gsl_vector_view b;
  gsl_vector *x_temp;

  for (i = 0; i < iiter; i++)
    b1[i] = MNA_right_vector[i * jiter];
  
  init(&A, &b, MNA_matrix, b1);

  x_temp = gsl_vector_alloc(iiter);
  LUDecomposition(A, b, x_temp);
  
  printf("To x mesa einai:\n");
  for (i = 0; i < iiter; i++)
  {
    xLU[i] = x_temp->data[i];
    printf("%lf\n", xLU[i]);
  }
}

void init_starting_value_Cholesky(double *MNA_matrix, double *MNA_right_vector, double *xLU, double step_time, double fin_time)
{
  int i;
  size_t iiter = mnac + hash_id;
  size_t jiter = (fin_time / step_time) + 1;
  gsl_matrix_view A;
  gsl_vector_view b;
  gsl_vector *x_temp;
  double b1[iiter];
  
  for (i = 0; i < iiter; i++)
    b1[i] = MNA_right_vector[i * jiter];
  
  Cholesky_init(&A, &b, MNA_matrix, b1);
  x_temp = gsl_vector_alloc(iiter);
  
  CholeskyDecomposition(A, b, x_temp);
  
  for (i = 0; i < iiter; i++)
    xLU[i] = x_temp->data[i];
}