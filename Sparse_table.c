#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "hash_table.h"
#include "csparse.h"

void create_Sparse_matrix(hashtable_t *hashtable,cs *Sparse_MNA_matrix,double *MNA_right_vector, int nonZeros){
  
  int counter=0;
  size_t group2_counter=0;
  
  printf("Creating sparse matrix %ld size\n", mnac + hash_id);
  
  elements *check;
  for(check = root; check != NULL; check = check->nextPtr){
    if(check->type == 'r'){
      if(ht_get(hashtable, check->plus_d_c)==0){
	Sparse_MNA_matrix->i[counter] = ht_get(hashtable, check->minus_g_b)-1;
	Sparse_MNA_matrix->p[counter] = ht_get(hashtable, check->minus_g_b)-1;
	Sparse_MNA_matrix->x[counter] = 1 / check->value;
	counter++;
      }
      else if(ht_get(hashtable, check->minus_g_b)==0){
	Sparse_MNA_matrix->i[counter] = ht_get(hashtable, check->plus_d_c)-1;
	Sparse_MNA_matrix->p[counter] = ht_get(hashtable, check->plus_d_c)-1;
	Sparse_MNA_matrix->x[counter] = 1 / check->value;
	counter++;
      }
      else{
	Sparse_MNA_matrix->i[counter] = ht_get(hashtable, check->minus_g_b)-1;
	Sparse_MNA_matrix->p[counter] = ht_get(hashtable, check->minus_g_b)-1;
	Sparse_MNA_matrix->x[counter] = 1 / check->value;
	counter++;
	Sparse_MNA_matrix->i[counter] = ht_get(hashtable, check->plus_d_c)-1;
	Sparse_MNA_matrix->p[counter] = ht_get(hashtable, check->plus_d_c)-1;
	Sparse_MNA_matrix->x[counter] = 1 / check->value;
	counter++;
	Sparse_MNA_matrix->i[counter] = ht_get(hashtable, check->minus_g_b)-1;
	Sparse_MNA_matrix->p[counter] = ht_get(hashtable, check->plus_d_c)-1;
	Sparse_MNA_matrix->x[counter] = -1 / check->value;
	counter++;
	Sparse_MNA_matrix->i[counter] = ht_get(hashtable, check->plus_d_c)-1;
	Sparse_MNA_matrix->p[counter] = ht_get(hashtable, check->minus_g_b)-1;
	Sparse_MNA_matrix->x[counter] = -1 / check->value;
	counter++;
      }
    }
    else if(check->type == 'v'){
      if(ht_get(hashtable, check->plus_d_c)==0){
	Sparse_MNA_matrix->i[counter] = hash_id+group2_counter;
	Sparse_MNA_matrix->p[counter] =  ht_get(hashtable, check->minus_g_b)-1;
	Sparse_MNA_matrix->x[counter] = -1;
	counter++;
	Sparse_MNA_matrix->i[counter] = ht_get(hashtable, check->minus_g_b)-1;
	Sparse_MNA_matrix->p[counter] =  hash_id+group2_counter;
	Sparse_MNA_matrix->x[counter] = -1;
	counter++;
      }
      else if(ht_get(hashtable, check->minus_g_b)==0){
	Sparse_MNA_matrix->i[counter] = hash_id+group2_counter;
	Sparse_MNA_matrix->p[counter] =  ht_get(hashtable, check->plus_d_c)-1;
	Sparse_MNA_matrix->x[counter] = 1;
	counter++;
	Sparse_MNA_matrix->i[counter] = ht_get(hashtable, check->plus_d_c)-1;
	Sparse_MNA_matrix->p[counter] =  hash_id+group2_counter;
	Sparse_MNA_matrix->x[counter] = 1;
	counter++;
      }
      else{
	Sparse_MNA_matrix->i[counter] = hash_id+group2_counter;
	Sparse_MNA_matrix->p[counter] =  ht_get(hashtable, check->minus_g_b)-1;
	Sparse_MNA_matrix->x[counter] = -1;
	counter++;
	Sparse_MNA_matrix->i[counter] = ht_get(hashtable, check->minus_g_b)-1;
	Sparse_MNA_matrix->p[counter] =  hash_id+group2_counter;
	Sparse_MNA_matrix->x[counter] = -1;
	counter++;
	Sparse_MNA_matrix->i[counter] = hash_id+group2_counter;
	Sparse_MNA_matrix->p[counter] =  ht_get(hashtable, check->plus_d_c)-1;
	Sparse_MNA_matrix->x[counter] = 1;
	counter++;
	Sparse_MNA_matrix->i[counter] = ht_get(hashtable, check->plus_d_c)-1;
	Sparse_MNA_matrix->p[counter] =  hash_id+group2_counter;
	Sparse_MNA_matrix->x[counter] = 1;
	counter++;
      }
      MNA_right_vector[hash_id+group2_counter] += check->value;
      group2_counter++;
    }
    else if(check->type == 'l'){
      if(ht_get(hashtable, check->plus_d_c)==0){
	Sparse_MNA_matrix->i[counter] = hash_id+group2_counter;
	Sparse_MNA_matrix->p[counter] =  ht_get(hashtable, check->minus_g_b)-1;
	Sparse_MNA_matrix->x[counter] = -1;
	counter++;
	Sparse_MNA_matrix->i[counter] = ht_get(hashtable, check->minus_g_b)-1;
	Sparse_MNA_matrix->p[counter] =  hash_id+group2_counter;
	Sparse_MNA_matrix->x[counter] = -1;
	counter++;
      }
      else if(ht_get(hashtable, check->minus_g_b)==0){
	Sparse_MNA_matrix->i[counter] = hash_id+group2_counter;
	Sparse_MNA_matrix->p[counter] =  ht_get(hashtable, check->plus_d_c)-1;
	Sparse_MNA_matrix->x[counter] = 1;
	counter++;
	Sparse_MNA_matrix->i[counter] = ht_get(hashtable, check->plus_d_c)-1;
	Sparse_MNA_matrix->p[counter] =  hash_id+group2_counter;
	Sparse_MNA_matrix->x[counter] = 1;
	counter++;
      }
      else{
	Sparse_MNA_matrix->i[counter] = hash_id+group2_counter;
	Sparse_MNA_matrix->p[counter] =  ht_get(hashtable, check->minus_g_b)-1;
	Sparse_MNA_matrix->x[counter] = -1;
	counter++;
	Sparse_MNA_matrix->i[counter] = ht_get(hashtable, check->minus_g_b)-1;
	Sparse_MNA_matrix->p[counter] =  hash_id+group2_counter;
	Sparse_MNA_matrix->x[counter] = -1;
	counter++;
	Sparse_MNA_matrix->i[counter] = hash_id+group2_counter;
	Sparse_MNA_matrix->p[counter] =  ht_get(hashtable, check->plus_d_c)-1;
	Sparse_MNA_matrix->x[counter] = 1;
	counter++;
	Sparse_MNA_matrix->i[counter] = ht_get(hashtable, check->plus_d_c)-1;
	Sparse_MNA_matrix->p[counter] =  hash_id+group2_counter;
	Sparse_MNA_matrix->x[counter] = 1;
	counter++;
      }
      group2_counter++;
    }
    else if(check->type == 'i'){
      if(ht_get(hashtable, check->plus_d_c)==0){
	MNA_right_vector[ht_get(hashtable, check->minus_g_b)-1] += check->value;
      }
      else if(atoi(check->minus_g_b)==0){
	MNA_right_vector[ht_get(hashtable, check->plus_d_c)-1] += -check->value;
      }
      else{
	MNA_right_vector[ht_get(hashtable, check->minus_g_b)-1] += check->value;
	MNA_right_vector[ht_get(hashtable, check->plus_d_c)-1] += -check->value;
      }
    }
  }
  Sparse_MNA_matrix->nz = nonZeros;
}