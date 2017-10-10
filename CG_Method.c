#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hash_table.h"


void CG_Method(double *MNA_Matrix, double *b, double *x, double itol, int x_length){
  int i =0,j=0;
  double rho,rho1,beta,alpha;
  double *M = (double *) malloc(x_length*sizeof(double));
  double *z = (double *) malloc(x_length*sizeof(double));
  double *p = (double *) malloc(x_length*sizeof(double));
  double *q = (double *) calloc(x_length,sizeof(double));
  for(i=0;i<x_length;i++){
    if(MNA_Matrix[i*x_length+i]==0.0){
      M[i] = 1;
    }
    else{
      M[i] = MNA_Matrix[i*x_length+i];
    }
  }
  for(i=0;i<x_length;i++){
    x[i] = 0;
  }
  int iter = 0;
  double *r = (double *) calloc(x_length,sizeof(double));
  Calculate_Init_r(MNA_Matrix,b,x,r,x_length);
  //printf("%f\n",norm2(r,3)/norm2(b,3));
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
    Product_Ap(MNA_Matrix,p,q,x_length);
    alpha = rho / dot_Product(p,q,x_length);
    Multiply_Add(x,p,x,alpha,x_length);
    Multiply_Reduce(r,q,r,alpha,x_length);
  }
  //free(M);
  //free(z);
  //free(p);
  //free(q);
}

void Calculate_Init_r(double *MNA_Matrix,double *MNA_right_vector,double *x,double *r,int r_length){
  int i,j;
  for (i=0;i<r_length;i++){
    r[i]=0;
  }
  for (i=0;i<r_length;i++){
    r[i] = MNA_right_vector[i];
    for(j=0;j<r_length;j++){
       r[i] += MNA_Matrix[i*r_length+j]*x[j];
    }
  }
}

double norm2(double *vector,int length){
  double res;
  int i;
  res =0;
  for(i=0;i<length;i++){
    res += vector[i]*vector[i];
  }
  return sqrt(res);
}

void Calculate_z(double *M,double *r,double *z,int length){
  int i=0;
  for(i=0;i<length;i++){
    z[i] = r[i] / M[i];
  }
}

double dot_Product(double *x,double *y, int length){
  double res = 0;
  int i=0;
  for(i=0;i<length;i++){
    res += x[i]*y[i];
  }
  return res;
}

void Multiply_Add(double *x, double *y, double *res, double val, int length){
  int i;
  for(i=0;i<length;i++){
    res[i] = x[i] + val * y[i];
  }
}

void Multiply_Reduce(double *x, double *y, double *res, double val, int length){
  int i;
  for(i=0;i<length;i++){
    res[i] = x[i] - val * y[i];
  }
}

void Product_Ap(double *A, double *p, double *q, int length){
  //double *temp = (double *) malloc(length*sizeof(double));
  int i,j;
  for(i=0;i<length;i++){
    q[i]=0;
  }
  for(i=0;i<length;i++){
    for(j=0;j<length;j++){
      q[i] += A[i*length + j]*p[j];
    }
  }
}
