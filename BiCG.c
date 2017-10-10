#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "hash_table.h"

void BiCG(double *A, double *b, double *x, double itol)
{
  double *r = malloc(sizeof(double)*(hash_id + mnac));
  double *rt = malloc(sizeof(double)*(hash_id + mnac));
  double *p = malloc(sizeof(double)*(hash_id + mnac));
  double *pt = malloc(sizeof(double)*(hash_id + mnac));
  double *M = malloc(sizeof(double)*(hash_id + mnac));
  double *Mt = malloc(sizeof(double)*(hash_id + mnac));
  double *q = malloc(sizeof(double)*(hash_id + mnac));
  double *qt = malloc(sizeof(double)*(hash_id + mnac));
  double *z = malloc(sizeof(double)*(hash_id + mnac));
  double *zt = malloc(sizeof(double)*(hash_id + mnac));
  double *At= malloc(sizeof(double)*((hash_id + mnac)*(hash_id + mnac)));
  int i = 0;
  int j = 0;
  int iter = 0;
  int n = hash_id + mnac;
  double rho = 0.0;
  double rho1 = 1.0;
  double beta, alpha, omega;
  
  // Initialiazations
  for (i = 0; i < hash_id + mnac; i++)
  {
    for (j = 0; j < hash_id + mnac; j++)
    {
      At[i * (mnac + hash_id) + j] = A[j * (mnac + hash_id) + i];
    }
  }
  for (i = 0; i < hash_id + mnac; i++)
  {
    z[i] = 0;
    zt[i] = 0;
    q[i] = 0;
    qt[i] = 0;
    if (fabs(A[i*(hash_id + mnac) + i]) < 0.00001)
    {
      M[i] = 1.0;
      Mt[i] = 1.0;
    }
    else
    {
      M[i] = A[i*(hash_id + mnac) + i];
      Mt[i] = A[i*(hash_id + mnac) + i];
    }
    x[i] = 0;
    r[i] = b[i];
    rt[i] = b[i];
    p[i] = 0;
    pt[i] = 0;
  }
  while (normC(r)/normC(b) > itol && iter < n)
  {
    iter = iter + 1;
    
    solveForBiCG(M, z, r);
    solveForBiCG(Mt, zt, rt);
    
    // rho = zt*rt
    rho = 0.0;
    for (i = 0; i < mnac + hash_id; i++)
      rho += z[i]*rt[i];
    if (fabs(rho) < 1e-14)
      break;
    if (iter == 1)
    {
      // p = z
      // pt = zt
      for (i = 0; i < mnac + hash_id; i++)
      {
	p[i] = z[i];
	pt[i] = zt[i];
      }
    }
    else
    {
      beta = rho/rho1;
      // p = z + beta * p
      // pt = zt + beta * pt
      for(i = 0;i < mnac + hash_id; i++)
      {
	p[i] = z[i] + beta * p[i];
	pt[i] = zt[i] + beta * pt[i];
      }
    }
    rho1 = rho;
    // q = A*p
    calculateQ(q, A, p);
    // qt = A*pt
    calculateQ(qt, At, pt);
    omega = 0.0;
    // omega = pt * q
    for (i = 0; i < mnac+hash_id; i++)
    {
      omega += pt[i]*q[i];
    }
    if (fabs(omega)< 1e-14)
      break;
    alpha = rho/omega;
    // x = x + alpha * p
    // r = r - alpha * q
    // rt = rt - alpha * qt
    for(i = 0; i < mnac + hash_id; i++)
    {
      x[i] += alpha * p[i];
      r[i] -= (alpha * q[i]);
      rt[i] -= (alpha * qt[i]);
    }
  }
}

double normC(double *mat)
{
  int i = 0;
  double sum = 0.0;
  for (i = 0; i < hash_id + mnac; i++)
  {
    sum += (mat[i]*mat[i]);
  }
  return sqrt(sum);
}

void solveForBiCG(double *M, double *z, double *r)
{
  int i = 0;
  for (i = 0; i < mnac + hash_id; i++)
  {
    z[i] = r[i]/M[i];
  }
}
  
void calculateQ(double *q, double *A, double *p)
{
  int i, j;
  for (i = 0; i < mnac + hash_id; i++)
  {
    q[i] = 0;
  }
  for (i = 0; i < mnac + hash_id; i++)
  {
    for (j = 0; j < mnac + hash_id; j++)
    {
      q[i] += A[i*(hash_id + mnac) + j] * p[j];
    }
  }
}

