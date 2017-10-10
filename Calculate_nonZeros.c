#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "hash_table.h"


int Calculate_nonZeros(hashtable_t *hashtable){
  int res=0;
  elements *check;
  for(check = root; check != NULL; check = check->nextPtr){
    if(check->type == 'r'){
      if((ht_get(hashtable, check->plus_d_c)==0)||(ht_get(hashtable, check->minus_g_b)==0)){
	res += 1;
      }
      else{
	res += 4;
      }
    }
     else if(check->type == 'v'){
       if((ht_get(hashtable, check->plus_d_c)==0)||(ht_get(hashtable, check->minus_g_b)==0)){
	 res += 2;
       }
       else{
	 res += 4;
       }
     }
     else if(check->type == 'l'){
       if((ht_get(hashtable, check->plus_d_c)==0)||(ht_get(hashtable, check->minus_g_b)==0)){
	 res += 2;
       }
       else{
	 res += 4;
       }
     }
  }
 
  return res;
}