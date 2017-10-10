#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hash_table.h"
#include "csparse.h"

void CG_Method_Sparse(cs *Sparse_MNA_matrix, double *b, double *x, double itol, int x_length){
  int i,j;
  double rho,rho1,beta,alpha;
  double *M = (double *) malloc(x_length*sizeof(double));
  double *z = (double *) malloc(x_length*sizeof(double));
  double *p = (double *) malloc(x_length*sizeof(double));
  double *q = (double *) calloc(x_length,sizeof(double));
  for(i=0;i<x_length;i++){
    M[i] = Sparse_MNA_matrix->x[Sparse_MNA_matrix->p[i]];
    //printf("%f\n",M[i]);
  }
  for(i=0;i<x_length;i++){
    x[i] = 0;
  }
  int iter = 0;
  double *r = (double *) calloc(x_length,sizeof(double));
  for (i=0;i<x_length;i++){
    r[i]=b[i];
    //printf("%f\n",r[i]);
  }
  //printf("%f\n",norm2(r,x_length)/norm2(b,x_length));
  while(((norm2(r,x_length)/norm2(b,x_length))>itol) && (iter<x_length)){
    iter += 1;
    Calculate_z(M,r,z,x_length);
    rho = dot_Product(r,z,x_length);
    if(iter == 1){
      memcpy(p,z,x_length*sizeof(double));
    }
    else{
      beta = rho/rho1;
      Multiply_Add(z,p,p,beta,x_length);
    }
    rho1 = rho;
    for(i=0;i<x_length;i++){
      q[i]=0;
    }
    int res = cs_gaxpy(Sparse_MNA_matrix, p, q);
    alpha = rho / dot_Product(p,q,x_length);
    Multiply_Add(x,p,x,alpha,x_length);
    Multiply_Reduce(r,q,r,alpha,x_length);
  }

  for (i = 0; i < mnac+hash_id; i++){
    printf("%lf\n", x[i]);
  } 
}