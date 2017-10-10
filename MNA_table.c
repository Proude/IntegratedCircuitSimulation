#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "hash_table.h"

void create_MNA_matrix(hashtable_t *hashtable, double *MNA_matrix, double *MNA_right_vector){
  
  size_t i,j;
  size_t group2_counter=0;
  
  printf("Creating matrix %ld size\n", mnac + hash_id);
  for(i=0; i < mnac+hash_id; i++){
    for(j=0; j < mnac+hash_id; j++){
      MNA_matrix[i * (mnac+hash_id) + j] = 0;
    }
  }
  for(i=0; i < mnac+hash_id; i++){
    MNA_right_vector[i] = 0;
  }
  elements *check;
  for(check = root; check != NULL; check = check->nextPtr){
    if(check->type == 'r'){
      if(ht_get(hashtable, check->plus_d_c)==0){
	MNA_matrix[(ht_get(hashtable, check->minus_g_b)-1) * (mnac+hash_id) + ht_get(hashtable, check->minus_g_b)-1] += 1 / check->value;
      }
      else if(ht_get(hashtable, check->minus_g_b)==0){
	MNA_matrix[(ht_get(hashtable, check->plus_d_c)-1) * (mnac+hash_id) + ht_get(hashtable, check->plus_d_c)-1] += 1 / check->value;
      }
      else{
	MNA_matrix[(ht_get(hashtable, check->minus_g_b)-1) * (mnac+hash_id) + ht_get(hashtable, check->minus_g_b)-1] += 1 / check->value;
	MNA_matrix[(ht_get(hashtable, check->plus_d_c)-1) * (mnac+hash_id) + ht_get(hashtable, check->plus_d_c)-1] += 1 / check->value;
	MNA_matrix[(ht_get(hashtable, check->minus_g_b)-1) * (mnac+hash_id) + ht_get(hashtable, check->plus_d_c)-1] -= 1 / check->value;
	MNA_matrix[(ht_get(hashtable, check->plus_d_c)-1) * (mnac+hash_id) + ht_get(hashtable, check->minus_g_b)-1] -= 1 / check->value;
      }
    }
    else if(check->type == 'v'){
      if(ht_get(hashtable, check->plus_d_c)==0){
	MNA_matrix[(hash_id+group2_counter) * (mnac+hash_id) + ht_get(hashtable, check->minus_g_b)-1] = -1;
	MNA_matrix[(ht_get(hashtable, check->minus_g_b)-1) * (mnac+hash_id) + hash_id+group2_counter] = -1;
      }
      else if(ht_get(hashtable, check->minus_g_b)==0){
	MNA_matrix[(hash_id+group2_counter) * (mnac+hash_id) + ht_get(hashtable, check->plus_d_c)-1] = 1;
	MNA_matrix[(ht_get(hashtable, check->plus_d_c)-1) * (mnac+hash_id) + hash_id+group2_counter] = 1;
      }
      else{
	MNA_matrix[(hash_id+group2_counter) * (mnac+hash_id) + ht_get(hashtable, check->minus_g_b)-1] = -1;
	MNA_matrix[(ht_get(hashtable, check->minus_g_b)-1) * (mnac+hash_id) + hash_id+group2_counter] = -1;
	MNA_matrix[(hash_id+group2_counter) * (mnac+hash_id) + ht_get(hashtable, check->plus_d_c)-1] = 1;
	MNA_matrix[(ht_get(hashtable, check->plus_d_c)-1) * (mnac+hash_id) + hash_id+group2_counter] = 1;
      }
      MNA_right_vector[hash_id+group2_counter] += check->value;
      group2_counter++;
    }
    else if(check->type == 'l'){
      if(ht_get(hashtable, check->plus_d_c)==0){
	MNA_matrix[(hash_id+group2_counter) * (mnac+hash_id) + ht_get(hashtable, check->minus_g_b)-1] = -1;
	MNA_matrix[(ht_get(hashtable, check->minus_g_b)-1) * (mnac+hash_id) + hash_id+group2_counter] = -1;
      }
      else if(ht_get(hashtable, check->minus_g_b)==0){
	MNA_matrix[(hash_id+group2_counter) * (mnac+hash_id) + ht_get(hashtable, check->plus_d_c)-1] = 1;
	MNA_matrix[(ht_get(hashtable, check->plus_d_c)-1) * (mnac+hash_id) + hash_id+group2_counter] = 1;
      }
      else{
	MNA_matrix[(hash_id+group2_counter) * (mnac+hash_id) + ht_get(hashtable, check->minus_g_b)-1] = -1;
	MNA_matrix[(ht_get(hashtable, check->minus_g_b)-1) * (mnac+hash_id) + hash_id+group2_counter] = -1;
	MNA_matrix[(hash_id+group2_counter) * (mnac+hash_id) + ht_get(hashtable, check->plus_d_c)-1] = 1;
	MNA_matrix[(ht_get(hashtable, check->plus_d_c)-1) * (mnac+hash_id) + hash_id+group2_counter] = 1;
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
  /*printf("Matrix A\n");
  for(i=0; i < mnac+hash_id; i++){
    printf("\n");
    for(j=0; j < mnac+hash_id; j++){
      printf("%f\t",MNA_matrix[i * (mnac+hash_id) + j]);
    }
  }
  printf("\nMatrix b\n");
  printf("\n");
  for(i=0;i < mnac+hash_id; i++){
    printf("%f\n",MNA_right_vector[i]);
  }*/
}