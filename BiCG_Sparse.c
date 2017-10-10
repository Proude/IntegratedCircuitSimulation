#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hash_table.h"
#include "csparse.h"

void BiCG_Sparse(cs *A, double *b, double *x, double itol)
{
  int iter = 0;
  int i = 0;
  int n = mnac + hash_id;
  double rho = 0.0;
  double rho1 = 0.0;
  double beta, alpha, omega;
  double *z = malloc((mnac + hash_id)*sizeof(double));
  double *zt = malloc((mnac + hash_id) * sizeof(double));
  double *q = malloc((mnac + hash_id) * sizeof(double));
  double *qt = malloc((mnac + hash_id) * sizeof(double));
  double *r = malloc((mnac + hash_id) * sizeof(double));
  double *rt = malloc((mnac + hash_id) * sizeof(double));
  double *p = malloc((mnac + hash_id) * sizeof(double));
  double *pt = malloc((mnac + hash_id) * sizeof(double));
  double *M = malloc((mnac + hash_id) * sizeof(double));
  double *Mt = malloc((mnac + hash_id) * sizeof(double));
  
  //initializations
  for (i = 0; i < hash_id + mnac; i++)
  {
    z[i] = 0;
    zt[i] = 0;
    q[i] = 0;
    qt[i] = 0;
    if (fabs(A->x[A->p[i]]) < 0.00001)
    {
      M[i] = 1.0;
      Mt[i] = 1.0;
    }
    else
    {
      M[i] = A->x[A->p[i]];
      Mt[i] = A->x[A->p[i]];
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
    for (i = 0; i < mnac + hash_id; i++)
    {
      q[i] = 0;
      qt[i] = 0;
    }
    // q = A*p
    if (!cs_gaxpy(A, p, q))
    {
      printf("Multiplication Failed! Program exits\n");
      exit(0);
    }
    // qt = A*pt
    if (!cs_gaxpy(A, pt, qt))
    {
      printf("Multiplication Failed! Program exits\n");
      exit(0);
    }
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
  
  for (i = 0; i < mnac + hash_id; i++)
  {
    printf("%lf\n", x[i]);
  }
    
}