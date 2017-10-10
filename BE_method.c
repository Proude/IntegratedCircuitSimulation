#include <stdio.h>
#include <stdlib.h>
#include "hash_table.h"

void multiply_h(double *C_new, double *C, int h);
void add_matrix(double *BE, double *C, double *G);
void multiply_vector_matrix(double *right_vector, double *C, double *x, double *et);

void BE_method(double *BE,double *right_vector, double *G, double *C, double *xtkprev,double *etk, int h_val){
  double *C_new;
  C_new = calloc((mnac+hash_id)*(mnac+hash_id),sizeof(double));
  multiply_h(C_new,C,h_val);
  add_matrix(BE,G,C_new);
  multiply_vector_matrix(right_vector,C_new, xtkprev, etk);
}

void multiply_h(double *C_new, double *C, int h){
  int i;
  for(i=0;i<(mnac+hash_id)*(mnac+hash_id);i++){
    C_new[i] = C[i] / h;
    //printf("%f\t %f\n",C_new[i],C[i]);
  }
}

void add_matrix(double *BE, double *G, double *C){
  int i,j;
  //printf("G + D/h\n");
  for(i=0;i<mnac+hash_id;i++){
    for(j=0;j<mnac+hash_id;j++){
      BE[i*(mnac+hash_id) + j] = C[i*(mnac+hash_id) + j] + G[i*(mnac+hash_id) + j];
      //printf("%f\t",BE[i*(mnac+hash_id) + j]);
    }
    //printf("\n");
  }
}

void multiply_vector_matrix(double *right_vector, double *C, double *x, double *et){
  int i,j;
  double res;
  for(i=0;i<mnac+hash_id;i++){
    res = 0;
    for(j=0;j<mnac+hash_id;j++){
      res += C[i*(mnac+hash_id) + j] * x[j];
      //printf("%f\t",res);
    }
    //printf("\n");
    right_vector[i] = res + et[i];
    //printf("e=%f riv=%.10f\n", et[i],right_vector[i]);
  }
}