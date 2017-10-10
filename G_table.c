#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "hash_table.h"
#include <math.h>

#define PI 3.14159265358979323846


void create_G_matrix(hashtable_t *hashtable, double *MNA_matrix, double *C_matrix,double *trans_right_vector, double step_time, double fin_time, int h_val, int numpwl){
  
  size_t i,j;
  size_t group2_counter=0;
  double time;
  double temp_h;
  
  printf("Creating matrix %ld size\n", mnac + hash_id);
  for(i=0; i < mnac+hash_id; i++){
    for(j=0; j < mnac+hash_id; j++){
      MNA_matrix[i * (mnac+hash_id) + j] = 0;
    }
  }
  for(i=0; i < mnac+hash_id; i++){
    for(j=0; j < mnac+hash_id; j++){
      C_matrix[i * (mnac+hash_id) + j] = 0;
    }
  }
  /*for(i=0; i < mnac+hash_id; i++){
    MNA_right_vector[i] = 0;
  }*/
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
      time = 0;
      for(temp_h=0; temp_h<h_val; temp_h++){
	int pos = (hash_id+group2_counter)*h_val + temp_h;
	double res = 0;
	//printf("time is %f\t",time);
	if((check->trans_type != NULL) && (strcmp(check->trans_type,"exp") == 0)){
	  if(time<=check->td1){
	    //printf("%f\n",check->i1);
	    trans_right_vector[pos] = check->i1;
	  }
	  else if((time > check->td1) && (time <= check->td2)){
	    res = check->i1 + (check->i2 - check->i1) * (1 - exp(-(time-check->td1)/check->tc1));
	    trans_right_vector[pos] = res;
	    //printf("%f\n",res);
	  }
	  else if(time > check->td2){
	    res = check->i1 + (check->i2 - check->i1) * (exp(-(time-check->td2)/check->tc2) - exp(-(time-check->td1)/check->tc1));
	    trans_right_vector[pos] = res;
	    //printf("%f\n",res);
	  }
	  else{
	    trans_right_vector[pos] = 5;
	  }
	}
	else if((check->trans_type != NULL) && (strcmp(check->trans_type,"sin") == 0)){
	  if(time <= check->tc1){
	    res = check->i1 + check->i2 * sin(2*PI*(check->tc2 / 360));
	    trans_right_vector[pos] = res;
	  }
	  else if(time > check->tc1){
	    res = check->i1 + check->i2 * sin(2*PI*(check->td1)*(time-check->tc1) + 2*PI*(check->tc2 / 360)) * exp(-(time-check->tc1)*check->td2);
	    trans_right_vector[pos] = res;
	  }
	  else{
	    trans_right_vector[pos] = 2.2;
	  }
	}
	else if((check->trans_type != NULL) && (strcmp(check->trans_type,"pulse") == 0)){
	  if(time<=check->td1){
	    res = check->i1;
	    trans_right_vector[pos] = res;
	  }
	  else if((time > check->td1) && (time <= check->td1 + check->tc1)){
	    res = ((check->i2 - check->i1) / check->tc1) * (time - check->td1) + check->i1;
	    trans_right_vector[pos] = res;
	  }
	  else if((time > check->td1 + check->tc1) && (time <= check->td1 + check->tc1 + check->tc2)){
	    res = check->i2;
	    trans_right_vector[pos] = res;
	  }
	  else if((time > check->td1 + check->tc1 + check->tc2) && (time <= check->td1 + check->tc1 + check->tc2 + check->td2) ){
	    res = ((check->i1 - check->i2) / check->td2) * (time - (check->td1 + check->tc1 + check->tc2)) + check->i2;
	    trans_right_vector[pos] = res;
	  }
	  else if((time > check->td1 + check->tc1 + check->tc2 + check->td2)){
	    res = check->i1;
	    trans_right_vector[pos] = res;
	  }
	  else{
	    trans_right_vector[pos] = 3.3;
	  }
	}
	else if((check->trans_type != NULL) && (strcmp(check->trans_type,"pwl") == 0)){
	  //printf("size of pwl is %d\t",numpwl);
	  //printf("first element of pwl is %f\n",check->pwl[2]);
	  if(time < check->pwl[0]){
	    res = check->pwl[1];
	  }
	  for(i=2; i<numpwl; i=i+2){
	    if((time < check->pwl[i]) && (time >= check->pwl[i-2])){
	      res = (((check->pwl[i+1]-check->pwl[i-1]) / (check->pwl[i] - check->pwl[i-2])) *(time - check->pwl[i-2]) + check->pwl[i-1]);
	    }
	  }
	  if(time >= check->pwl[numpwl-2]){
	    res = check->pwl[numpwl-1];
	  }
	  trans_right_vector[pos] = res;
	}
	else{
	  trans_right_vector[pos] = check->value;
	}
	time = time + step_time;
	//temp_h++;
      }
      //MNA_right_vector[hash_id+group2_counter] += check->value;
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
      C_matrix[(hash_id+group2_counter) * (mnac+hash_id) + hash_id+group2_counter] -= check->value;
      group2_counter++;
    }
    if(check->type == 'c'){
      if(ht_get(hashtable, check->plus_d_c)==0){
	C_matrix[(ht_get(hashtable, check->minus_g_b)-1) * (mnac+hash_id) + ht_get(hashtable, check->minus_g_b)-1] += check->value;
      }
      else if(ht_get(hashtable, check->minus_g_b)==0){
	C_matrix[(ht_get(hashtable, check->plus_d_c)-1) * (mnac+hash_id) + ht_get(hashtable, check->plus_d_c)-1] += check->value;
      }
      else{
	C_matrix[(ht_get(hashtable, check->minus_g_b)-1) * (mnac+hash_id) + ht_get(hashtable, check->minus_g_b)-1] += check->value;
	C_matrix[(ht_get(hashtable, check->plus_d_c)-1) * (mnac+hash_id) + ht_get(hashtable, check->plus_d_c)-1] += check->value;
	C_matrix[(ht_get(hashtable, check->minus_g_b)-1) * (mnac+hash_id) + ht_get(hashtable, check->plus_d_c)-1] -= check->value;
	C_matrix[(ht_get(hashtable, check->plus_d_c)-1) * (mnac+hash_id) + ht_get(hashtable, check->minus_g_b)-1] -= check->value;
      }
    }
    else if(check->type == 'i'){
      if(ht_get(hashtable, check->plus_d_c)==0){
	//MNA_right_vector[ht_get(hashtable, check->minus_g_b)-1] += check->value;
	time = 0;
      for(temp_h=0; temp_h<h_val; temp_h++){
	int pos = (ht_get(hashtable, check->minus_g_b)-1)*h_val + temp_h;
	double res = 0;
	//printf("time is %f\t",time);
	if((check->trans_type != NULL) && (strcmp(check->trans_type,"exp") == 0)){
	  if(time<=check->td1){
	    //printf("%f\n",check->i1);
	    trans_right_vector[pos] -= check->i1;
	  }
	  else if((time > check->td1) && (time <= check->td2)){
	    res = check->i1 + (check->i2 - check->i1) * (1 - exp(-(time-check->td1)/check->tc1));
	    trans_right_vector[pos] -= res;
	    //printf("%f\n",res);
	  }
	  else if(time > check->td2){
	    res = check->i1 + (check->i2 - check->i1) * (exp(-(time-check->td2)/check->tc2) - exp(-(time-check->td1)/check->tc1));
	    trans_right_vector[pos] -= res;
	    //printf("%f\n",res);
	  }
	  else{
	    trans_right_vector[pos] = 5;
	  }
	}
	else if((check->trans_type != NULL) && (strcmp(check->trans_type,"sin") == 0)){
	  if(time <= check->tc1){
	    res = check->i1 + check->i2 * sin(2*PI*(check->tc2 / 360));
	    trans_right_vector[pos] -= res;
	  }
	  else if(time > check->tc1){
	    res = check->i1 + check->i2 * sin(2*PI*(check->td1)*(time-check->tc1) + 2*PI*(check->tc2 / 360)) * exp(-(time-check->tc1)*check->td2);
	    trans_right_vector[pos] -= res;
	  }
	  else{
	    trans_right_vector[pos] = 2.2;
	  }
	}
	else if((check->trans_type != NULL) && (strcmp(check->trans_type,"pulse") == 0)){
	  if(time<=check->td1){
	    res = check->i1;
	    trans_right_vector[pos] -= res;
	  }
	  else if((time > check->td1) && (time <= check->td1 + check->tc1)){
	    res = ((check->i2 - check->i1) / check->tc1) * (time - check->td1) + check->i1;
	    trans_right_vector[pos] -= res;
	  }
	  else if((time > check->td1 + check->tc1) && (time <= check->td1 + check->tc1 + check->tc2)){
	    res = check->i2;
	    trans_right_vector[pos] -= res;
	  }
	  else if((time > check->td1 + check->tc1 + check->tc2) && (time <= check->td1 + check->tc1 + check->tc2 + check->td2) ){
	    res = ((check->i1 - check->i2) / check->td2) * (time - (check->td1 + check->tc1 + check->tc2)) + check->i2;
	    trans_right_vector[pos] -= res;
	  }
	  else if((time > check->td1 + check->tc1 + check->tc2 + check->td2)){
	    res = check->i1;
	    trans_right_vector[pos] -= res;
	  }
	  else{
	    trans_right_vector[pos] = 3.3;
	  }
	}
	else if((check->trans_type != NULL) && (strcmp(check->trans_type,"pwl") == 0)){
	  //printf("size of pwl is %d\t",numpwl);
	  //printf("first element of pwl is %f\n",check->pwl[2]);
	  if(time < check->pwl[0]){
	    res = check->pwl[1];
	  }
	  for(i=2; i<numpwl; i=i+2){
	    if((time < check->pwl[i]) && (time >= check->pwl[i-2])){
	      res = (((check->pwl[i+1]-check->pwl[i-1]) / (check->pwl[i] - check->pwl[i-2])) *(time - check->pwl[i-2]) + check->pwl[i-1]);
	    }
	  }
	  if(time >= check->pwl[numpwl-2]){
	    res = check->pwl[numpwl-1];
	  }
	  trans_right_vector[pos] -= res;
	}
	else{
	  trans_right_vector[pos] -= check->value;
	}
	time = time + step_time;
      }
      }
      else if(atoi(check->minus_g_b)==0){
	//MNA_right_vector[ht_get(hashtable, check->plus_d_c)-1] += -check->value;
	time = 0;
      for(temp_h=0; temp_h<h_val; temp_h++){
	int pos = (ht_get(hashtable, check->plus_d_c)-1)*h_val + temp_h;
	double res = 0;
	//printf("time is %f\t",time);
	if((check->trans_type != NULL) && (strcmp(check->trans_type,"exp") == 0)){
	  if(time<=check->td1){
	    //printf("%f\n",check->i1);
	    trans_right_vector[pos] = check->i1;
	  }
	  else if((time > check->td1) && (time <= check->td2)){
	    res = check->i1 + (check->i2 - check->i1) * (1 - exp(-(time-check->td1)/check->tc1));
	    trans_right_vector[pos] += res;
	    //printf("%f\n",res);
	  }
	  else if(time > check->td2){
	    res = check->i1 + (check->i2 - check->i1) * (exp(-(time-check->td2)/check->tc2) - exp(-(time-check->td1)/check->tc1));
	    trans_right_vector[pos] += res;
	    //printf("%f\n",res);
	  }
	  else{
	    trans_right_vector[pos] = 5;
	  }
	}
	else if((check->trans_type != NULL) && (strcmp(check->trans_type,"sin") == 0)){
	  if(time <= check->tc1){
	    res = check->i1 + check->i2 * sin(2*PI*(check->tc2 / 360));
	    trans_right_vector[pos] += res;
	  }
	  else if(time > check->tc1){
	    res = check->i1 + check->i2 * sin(2*PI*(check->td1)*(time-check->tc1) + 2*PI*(check->tc2 / 360)) * exp(-(time-check->tc1)*check->td2);
	    trans_right_vector[pos] += res;
	  }
	  else{
	    trans_right_vector[pos] = 2.2;
	  }
	}
	else if((check->trans_type != NULL) && (strcmp(check->trans_type,"pulse") == 0)){
	  if(time<=check->td1){
	    res = check->i1;
	    trans_right_vector[pos] += res;
	  }
	  else if((time > check->td1) && (time <= check->td1 + check->tc1)){
	    res = ((check->i2 - check->i1) / check->tc1) * (time - check->td1) + check->i1;
	    trans_right_vector[pos] += res;
	  }
	  else if((time > check->td1 + check->tc1) && (time <= check->td1 + check->tc1 + check->tc2)){
	    res = check->i2;
	    trans_right_vector[pos] += res;
	  }
	  else if((time > check->td1 + check->tc1 + check->tc2) && (time <= check->td1 + check->tc1 + check->tc2 + check->td2) ){
	    res = ((check->i1 - check->i2) / check->td2) * (time - (check->td1 + check->tc1 + check->tc2)) + check->i2;
	    trans_right_vector[pos] += res;
	  }
	  else if((time > check->td1 + check->tc1 + check->tc2 + check->td2)){
	    res = check->i1;
	    trans_right_vector[pos] += res;
	  }
	  else{
	    trans_right_vector[pos] = 3.3;
	  }
	}
	else if((check->trans_type != NULL) && (strcmp(check->trans_type,"pwl") == 0)){
	  //printf("size of pwl is %d\t",numpwl);
	  //printf("first element of pwl is %f\n",check->pwl[2]);
	  if(time < check->pwl[0]){
	    res = check->pwl[1];
	  }
	  for(i=2; i<numpwl; i=i+2){
	    if((time < check->pwl[i]) && (time >= check->pwl[i-2])){
	      res = (((check->pwl[i+1]-check->pwl[i-1]) / (check->pwl[i] - check->pwl[i-2])) *(time - check->pwl[i-2]) + check->pwl[i-1]);
	    }
	  }
	  if(time >= check->pwl[numpwl-2]){
	    res = check->pwl[numpwl-1];
	  }
	  trans_right_vector[pos] += res;
	}
	else{
	  trans_right_vector[pos] += check->value;
	}
	time = time + step_time;
      }
      }
      else{
	//MNA_right_vector[ht_get(hashtable, check->minus_g_b)-1] += check->value;
	//MNA_right_vector[ht_get(hashtable, check->plus_d_c)-1] += -check->value;time = 0;
	time = 0;
      for(temp_h=0; temp_h<h_val; temp_h++){
	int pos = (ht_get(hashtable, check->plus_d_c)-1)*h_val + temp_h;
	double res = 0;
	//printf("time is %f\t",time);
	if((check->trans_type != NULL) && (strcmp(check->trans_type,"exp") == 0)){
	  if(time<=check->td1){
	    //printf("%f\n",check->i1);
	    trans_right_vector[pos] = check->i1;
	  }
	  else if((time > check->td1) && (time <= check->td2)){
	    res = check->i1 + (check->i2 - check->i1) * (1 - exp(-(time-check->td1)/check->tc1));
	    trans_right_vector[pos] += res;
	    //printf("%f\n",res);
	  }
	  else if(time > check->td2){
	    res = check->i1 + (check->i2 - check->i1) * (exp(-(time-check->td2)/check->tc2) - exp(-(time-check->td1)/check->tc1));
	    trans_right_vector[pos] += res;
	    //printf("%f\n",res);
	  }
	  else{
	    trans_right_vector[pos] = 5;
	  }
	}
	else if((check->trans_type != NULL) && (strcmp(check->trans_type,"sin") == 0)){
	  if(time <= check->tc1){
	    res = check->i1 + check->i2 * sin(2*PI*(check->tc2 / 360));
	    trans_right_vector[pos] += res;
	  }
	  else if(time > check->tc1){
	    res = check->i1 + check->i2 * sin(2*PI*(check->td1)*(time-check->tc1) + 2*PI*(check->tc2 / 360)) * exp(-(time-check->tc1)*check->td2);
	    trans_right_vector[pos] += res;
	  }
	  else{
	    trans_right_vector[pos] = 2.2;
	  }
	}
	else if((check->trans_type != NULL) && (strcmp(check->trans_type,"pulse") == 0)){
	  if(time<=check->td1){
	    res = check->i1;
	    trans_right_vector[pos] += res;
	  }
	  else if((time > check->td1) && (time <= check->td1 + check->tc1)){
	    res = ((check->i2 - check->i1) / check->tc1) * (time - check->td1) + check->i1;
	    trans_right_vector[pos] += res;
	  }
	  else if((time > check->td1 + check->tc1) && (time <= check->td1 + check->tc1 + check->tc2)){
	    res = check->i2;
	    trans_right_vector[pos] += res;
	  }
	  else if((time > check->td1 + check->tc1 + check->tc2) && (time <= check->td1 + check->tc1 + check->tc2 + check->td2) ){
	    res = ((check->i1 - check->i2) / check->td2) * (time - (check->td1 + check->tc1 + check->tc2)) + check->i2;
	    trans_right_vector[pos] += res;
	  }
	  else if((time > check->td1 + check->tc1 + check->tc2 + check->td2)){
	    res = check->i1;
	    trans_right_vector[pos] += res;
	  }
	  else{
	    trans_right_vector[pos] = 3.3;
	  }
	}
	else if((check->trans_type != NULL) && (strcmp(check->trans_type,"pwl") == 0)){
	  //printf("size of pwl is %d\t",numpwl);
	  //printf("first element of pwl is %f\n",check->pwl[2]);
	  if(time < check->pwl[0]){
	    res = check->pwl[1];
	  }
	  for(i=2; i<numpwl; i=i+2){
	    if((time < check->pwl[i]) && (time >= check->pwl[i-2])){
	      res = (((check->pwl[i+1]-check->pwl[i-1]) / (check->pwl[i] - check->pwl[i-2])) *(time - check->pwl[i-2]) + check->pwl[i-1]);
	    }
	  }
	  if(time >= check->pwl[numpwl-2]){
	    res = check->pwl[numpwl-1];
	  }
	  trans_right_vector[pos] += res;
	}
	else{
	  trans_right_vector[pos] += check->value;
	}
	time = time + step_time;
      }
      time = 0;
      for(temp_h=0; temp_h<h_val; temp_h++){
	int pos = (ht_get(hashtable, check->minus_g_b)-1)*h_val + temp_h;
	double res = 0;
	//printf("time is %f\t",time);
	if((check->trans_type != NULL) && (strcmp(check->trans_type,"exp") == 0)){
	  if(time<=check->td1){
	    //printf("%f\n",check->i1);
	    trans_right_vector[pos] -= check->i1;
	  }
	  else if((time > check->td1) && (time <= check->td2)){
	    res = check->i1 + (check->i2 - check->i1) * (1 - exp(-(time-check->td1)/check->tc1));
	    trans_right_vector[pos] -= res;
	    //printf("%f\n",res);
	  }
	  else if(time > check->td2){
	    res = check->i1 + (check->i2 - check->i1) * (exp(-(time-check->td2)/check->tc2) - exp(-(time-check->td1)/check->tc1));
	    trans_right_vector[pos] -= res;
	    //printf("%f\n",res);
	  }
	  else{
	    trans_right_vector[pos] = 5;
	  }
	}
	else if((check->trans_type != NULL) && (strcmp(check->trans_type,"sin") == 0)){
	  if(time <= check->tc1){
	    res = check->i1 + check->i2 * sin(2*PI*(check->tc2 / 360));
	    trans_right_vector[pos] -= res;
	  }
	  else if(time > check->tc1){
	    res = check->i1 + check->i2 * sin(2*PI*(check->td1)*(time-check->tc1) + 2*PI*(check->tc2 / 360)) * exp(-(time-check->tc1)*check->td2);
	    trans_right_vector[pos] -= res;
	  }
	  else{
	    trans_right_vector[pos] = 2.2;
	  }
	}
	else if((check->trans_type != NULL) && (strcmp(check->trans_type,"pulse") == 0)){
	  if(time<=check->td1){
	    res = check->i1;
	    trans_right_vector[pos] -= res;
	  }
	  else if((time > check->td1) && (time <= check->td1 + check->tc1)){
	    res = ((check->i2 - check->i1) / check->tc1) * (time - check->td1) + check->i1;
	    trans_right_vector[pos] -= res;
	  }
	  else if((time > check->td1 + check->tc1) && (time <= check->td1 + check->tc1 + check->tc2)){
	    res = check->i2;
	    trans_right_vector[pos] -= res;
	  }
	  else if((time > check->td1 + check->tc1 + check->tc2) && (time <= check->td1 + check->tc1 + check->tc2 + check->td2) ){
	    res = ((check->i1 - check->i2) / check->td2) * (time - (check->td1 + check->tc1 + check->tc2)) + check->i2;
	    trans_right_vector[pos] -= res;
	  }
	  else if((time > check->td1 + check->tc1 + check->tc2 + check->td2)){
	    res = check->i1;
	    trans_right_vector[pos] -= res;
	  }
	  else{
	    trans_right_vector[pos] = 3.3;
	  }
	}
	else if((check->trans_type != NULL) && (strcmp(check->trans_type,"pwl") == 0)){
	  //printf("size of pwl is %d\t",numpwl);
	  //printf("first element of pwl is %f\n",check->pwl[2]);
	  if(time < check->pwl[0]){
	    res = check->pwl[1];
	  }
	  for(i=2; i<numpwl; i=i+2){
	    if((time < check->pwl[i]) && (time >= check->pwl[i-2])){
	      res = (((check->pwl[i+1]-check->pwl[i-1]) / (check->pwl[i] - check->pwl[i-2])) *(time - check->pwl[i-2]) + check->pwl[i-1]);
	    }
	  }
	  if(time >= check->pwl[numpwl-2]){
	    res = check->pwl[numpwl-1];
	  }
	  trans_right_vector[pos] -= res;
	}
	else{
	  trans_right_vector[pos] -= check->value;
	}
	time = time + step_time;
      }
      }
    }
  }
  /*printf("Matrix G\n");
  for(i=0; i < mnac+hash_id; i++){
    printf("\n");
    for(j=0; j < mnac+hash_id; j++){
      printf("%f\t",MNA_matrix[i * (mnac+hash_id) + j]);
    }
  }*/
  /*printf("\n");
  printf("Matrix C\n");
  for(i=0; i < mnac+hash_id; i++){
    printf("\n");
    for(j=0; j < mnac+hash_id; j++){
      printf("%.10f\t",C_matrix[i * (mnac+hash_id) + j]);
    }
  }*/
  /*for(i=0; i < mnac+hash_id; i++){
    printf("\n");
    for(j=0; j < h_val;j++){
      printf("%f\t",trans_right_vector[(i*h_val)+j]);
    }
  }*/
  printf("\n");
  
  /*printf("\nMatrix b\n");
  printf("\n");
  for(i=0;i < mnac+hash_id; i++){
    printf("%f\n",MNA_right_vector[i]);
  }*/
}