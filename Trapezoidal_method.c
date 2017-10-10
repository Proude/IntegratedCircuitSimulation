#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "hash_table.h"

void trapezoidal_LU(double *G_matrix, double *C_matrix, double *x_final, int h_val, double *trans_right_matrix)
{
  int i, j, z;
  double *A = NULL;
  double *A_temp = NULL;
  double *A_temp2 = NULL;
  double *b = NULL;
  double *b_temp = NULL;
  double *b_temp2 = NULL;
  //double *x_temp = NULL;
  gsl_matrix_view A_view;
  gsl_vector_view b_view;
  gsl_vector *x_temp_view;
  
  A = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
  A_temp = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
  A_temp2 = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
  b = (double *)calloc((mnac + hash_id), sizeof(double));
  b_temp = (double *)calloc((mnac + hash_id), sizeof(double));
  b_temp2 = (double *)calloc((mnac + hash_id), sizeof(double));
  //x_temp = (double *)calloc((mnac + hash_id), sizeof(double));
  
  // Compute the A matrix to prepare for solving
  for (i = 0; i < ((mnac + hash_id) * (mnac + hash_id)); i++)
  {
    A[i] = G_matrix[i] + ((2 / h_val) * C_matrix[i]);
    A_temp[i] = G_matrix[i] - ((2 / h_val) * C_matrix[i]);
  }
  
  for (j = 1; j < (h_val + 1); j++)
  {
    for (i = 0; i < (mnac + hash_id); i++)
      b[i] = trans_right_matrix[i * (h_val + 1) + j] + trans_right_matrix[i * (h_val + 1) + (j - 1)];
    
    for (i = 0; i < (mnac + hash_id); i++)
      b_temp[i] = 0;
    
    for (i = 0; i < (mnac + hash_id); i++)
      for (z = 0; z < (mnac + hash_id); z++)
	b_temp[i] += A_temp[(i * (mnac + hash_id)) + z]*x_final[(z * (h_val + 1)) + (j - 1)];
      
    for (i = 0; i < (mnac + hash_id); i++)
      b[i] = b[i] - b_temp[i];
    
    for (i = 0; i < (mnac + hash_id) * (mnac + hash_id); i++)
      A_temp2[i] = A[i];
    for (i = 0; i < (mnac + hash_id); i++)
      b_temp2[i] = b[i];
    
    init(&A_view, &b_view, A_temp2, b_temp2);
    x_temp_view = gsl_vector_alloc(mnac + hash_id);
    LUDecomposition(A_view, b_view, x_temp_view);

    printf("%d iteration:\n", j);
    for (i = 0; i < mnac + hash_id; i++)
    {
      x_final[i * (h_val + 1) + (j + 1)] = x_temp_view->data[i];
      printf("%lf \n", x_final[i * (h_val + 1) + (j + 1)]);
    }
    gsl_vector_free(x_temp_view);
  }
  
}

void trapezoidal_Cholesky(double *G_matrix, double *C_matrix, double *x_final, int h_val, double *trans_right_matrix)
{
  int i, j, z;
  double *A = NULL;
  double *A_temp = NULL;
  double *A_temp2 = NULL;
  double *b = NULL;
  double *b_temp = NULL;
  double *b_temp2 = NULL;
  //double *x_temp = NULL;
  gsl_matrix_view A_view;
  gsl_vector_view b_view;
  gsl_vector *x_temp_view;
  
  A = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
  A_temp = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
  A_temp2 = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
  b = (double *)calloc((mnac + hash_id), sizeof(double));
  b_temp = (double *)calloc((mnac + hash_id), sizeof(double));
  b_temp2 = (double *)calloc((mnac + hash_id), sizeof(double));
  //x_temp = (double *)calloc((mnac + hash_id), sizeof(double));
  
  // Compute the A matrix to prepare for solving
  for (i = 0; i < ((mnac + hash_id) * (mnac + hash_id)); i++)
  {
    A[i] = G_matrix[i] + ((2 / h_val) * C_matrix[i]);
    A_temp[i] = G_matrix[i] - ((2 / h_val) * C_matrix[i]);
  }
  
  for (j = 1; j < (h_val + 1); j++)
  {
    for (i = 0; i < (mnac + hash_id); i++)
      b[i] = trans_right_matrix[i * (h_val + 1) + j] + trans_right_matrix[i * (h_val + 1) + (j - 1)];
    
    for (i = 0; i < (mnac + hash_id); i++)
      b_temp[i] = 0;
    
    for (i = 0; i < (mnac + hash_id); i++)
      for (z = 0; z < (mnac + hash_id); z++)
	b_temp[i] += A_temp[(i * (mnac + hash_id)) + z]*x_final[(z * (h_val + 1)) + (j - 1)];
      
    for (i = 0; i < (mnac + hash_id); i++)
      b[i] = b[i] - b_temp[i];
    
    for (i = 0; i < (mnac + hash_id) * (mnac + hash_id); i++)
      A_temp2[i] = A[i];
    for (i = 0; i < (mnac + hash_id); i++)
      b_temp2[i] = b[i];
    
    Cholesky_init(&A_view, &b_view, A_temp2, b_temp2);
    x_temp_view = gsl_vector_alloc(mnac + hash_id);
    CholeskyDecomposition(A_view, b_view, x_temp_view);
  
    for (i = 0; i < mnac + hash_id; i++)
      x_final[i * (h_val + 1) + (j + 1)] = x_temp_view->data[i];
    
    printf("%d iteration:\n", j);
    for (i = 0; i < mnac + hash_id; i++)
      printf("%lf \n", x_final[i * (h_val + 1) + (j + 1)]);
    gsl_vector_free(x_temp_view);
  }
  
}

void trapezoidal_CG(double *G_matrix, double *C_matrix, double *x_final, int h_val, double *trans_right_matrix, double itol)
{
  int i, j, z;
  double *A = NULL;
  double *A_temp = NULL;
  //double *A_temp2 = NULL;
  double *b = NULL;
  double *b_temp = NULL;
  //double *b_temp2 = NULL;
  double *x_temp = NULL;
  
  A = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
  A_temp = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
  //A_temp2 = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
  b = (double *)calloc((mnac + hash_id), sizeof(double));
  b_temp = (double *)calloc((mnac + hash_id), sizeof(double));
  //b_temp2 = (double *)calloc((mnac + hash_id), sizeof(double));
  x_temp = (double *)calloc((mnac + hash_id), sizeof(double));
  
  // Compute the A matrix to prepare for solving
  for (i = 0; i < ((mnac + hash_id) * (mnac + hash_id)); i++)
  {
    A[i] = G_matrix[i] + ((2 / h_val) * C_matrix[i]);
    A_temp[i] = G_matrix[i] - ((2 / h_val) * C_matrix[i]);
  }
  
  for (j = 1; j < (h_val + 1); j++)
  {
    for (i = 0; i < (mnac + hash_id); i++)
      b[i] = trans_right_matrix[i * (h_val + 1) + j] + trans_right_matrix[i * (h_val + 1) + (j - 1)];
    
    for (i = 0; i < (mnac + hash_id); i++)
      b_temp[i] = 0;
    
    for (i = 0; i < (mnac + hash_id); i++)
      for (z = 0; z < (mnac + hash_id); z++)
	b_temp[i] += A_temp[(i * (mnac + hash_id)) + z]*x_final[(z * (h_val + 1)) + (j - 1)];
      
    for (i = 0; i < (mnac + hash_id); i++)
      b[i] = b[i] - b_temp[i];
    
    //for (i = 0; i < (mnac + hash_id) * (mnac + hash_id); i++)
    //  A_temp2[i] = A[i];
    //for (i = 0; i < (mnac + hash_id); i++)
    //  b_temp2[i] = b[i];
    
    CG_Method(A, b, x_temp, itol, mnac + hash_id);
  
    for (i = 0; i < mnac + hash_id; i++)
      x_final[i * (h_val + 1) + (j + 1)] = x_temp[i];
    printf("%d iteration:\n", j);
    for (i = 0; i < mnac + hash_id; i++)
      printf("%lf \n", x_final[i * (h_val + 1) + (j + 1)]);
  }
  
}

void trapezoidal_BiCG(double *G_matrix, double *C_matrix, double *x_final, int h_val, double *trans_right_matrix, double itol)
{
  int i, j, z;
  double *A = NULL;
  double *A_temp = NULL;
  //double *A_temp2 = NULL;
  double *b = NULL;
  double *b_temp = NULL;
  //double *b_temp2 = NULL;
  double *x_temp = NULL;
  
  A = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
  A_temp = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
  //A_temp2 = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
  b = (double *)calloc((mnac + hash_id), sizeof(double));
  b_temp = (double *)calloc((mnac + hash_id), sizeof(double));
  //b_temp2 = (double *)calloc((mnac + hash_id), sizeof(double));
  x_temp = (double *)calloc((mnac + hash_id), sizeof(double));
  
  // Compute the A matrix to prepare for solving
  for (i = 0; i < ((mnac + hash_id) * (mnac + hash_id)); i++)
  {
    A[i] = G_matrix[i] + ((2 / h_val) * C_matrix[i]);
    A_temp[i] = G_matrix[i] - ((2 / h_val) * C_matrix[i]);
  }
  
  for (j = 1; j < (h_val + 1); j++)
  {
    for (i = 0; i < (mnac + hash_id); i++)
      b[i] = trans_right_matrix[i * (h_val + 1) + j] + trans_right_matrix[i * (h_val + 1) + (j - 1)];
    
    for (i = 0; i < (mnac + hash_id); i++)
      b_temp[i] = 0;
    
    for (i = 0; i < (mnac + hash_id); i++)
      for (z = 0; z < (mnac + hash_id); z++)
	b_temp[i] += A_temp[(i * (mnac + hash_id)) + z]*x_final[(z * (h_val + 1)) + (j - 1)];
      
    for (i = 0; i < (mnac + hash_id); i++)
      b[i] = b[i] - b_temp[i];
    
    //for (i = 0; i < (mnac + hash_id) * (mnac + hash_id); i++)
    //  A_temp2[i] = A[i];
    //for (i = 0; i < (mnac + hash_id); i++)
    //  b_temp2[i] = b[i];
    
    BiCG(A, b, x_temp, itol);
  
    for (i = 0; i < mnac + hash_id; i++)
      x_final[i * (h_val + 1) + (j + 1)] = x_temp[i];
    printf("%d iteration:\n", j);
    for (i = 0; i < mnac + hash_id; i++)
      printf("%lf \n", x_final[i * (h_val + 1) + (j + 1)]);

  }
  
}