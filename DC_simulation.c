#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <string.h>
#include "hash_table.h"

void DC_init(gsl_matrix_view *A, gsl_vector_view *b, double *MNA_matrix, double *MNA_right_vector,char type, char var_name[], double val, hashtable_t *hashtable){
  elements *check;
  int group2_counter = 0;
  for(check = root; check != NULL; check = check->nextPtr){
    if((strcmp(check->name,var_name)==0)&&(check->type == type)){
      if(type == 'v'){
	MNA_right_vector[hash_id+group2_counter] = val;
      }
      else if(type =='i'){
	if(ht_get(hashtable, check->plus_d_c)==0){
	  MNA_right_vector[ht_get(hashtable, check->minus_g_b)-1] = val;
	}
	else if(atoi(check->minus_g_b)==0){
	  MNA_right_vector[ht_get(hashtable, check->plus_d_c)-1] = -val;
	}
	else{
	  MNA_right_vector[ht_get(hashtable, check->minus_g_b)-1] = val;
	  MNA_right_vector[ht_get(hashtable, check->plus_d_c)-1] = -val;
	}
      }
    }
    if((check->type == 'v')||(check->type=='l')){
      group2_counter++;
    }
  }
  *b = gsl_vector_view_array(MNA_right_vector, mnac+hash_id);
  *A = gsl_matrix_view_array(MNA_matrix, mnac+hash_id, mnac+hash_id);
}

void DC_simulation(gsl_matrix_view A, gsl_vector_view b, gsl_vector *x){
  int s;
  gsl_permutation *p = gsl_permutation_alloc (mnac+hash_id);

  gsl_linalg_LU_decomp (&A.matrix, p, &s);
  gsl_linalg_LU_solve (&A.matrix, p, &b.vector, x);
  gsl_permutation_free(p);
}
