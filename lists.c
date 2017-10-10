#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "hash_table.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
// List initialisation
///////////////////////////////////////////////////////////////////////////////////////////////////
void list_init()
{
  root = NULL;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
// This function is used to add elements to our list
//////////////////////////////////////////////////////////////////////////////////////////////////
void add_to_list(char *name, char type, char *plus_d_c, char *minus_g_b, char *source_e, char *body, double value, char *model_name, double L, double W, double area,
		 char *trans_type, double i1, double i2, double td1, double tc1, double td2, double tc2, double per, double *pwl, int size_of_pwl)
{
  int i;
  elements *temp;
  temp = (elements *)malloc(sizeof(elements));
  temp->name = (char *)malloc((strlen(name) + 1) * sizeof(char));
  //temp->type = (char *)malloc((strlen(name)+1)*sizeof(char));
  temp->plus_d_c = (char *)malloc((strlen(plus_d_c) + 1) * sizeof(char));
  temp->minus_g_b = (char *)malloc((strlen(minus_g_b) + 1) * sizeof(char));
  temp->source_e = (char *)malloc((strlen(source_e) + 1) * sizeof(char));
  temp->body = (char *)malloc((strlen(body) + 1)* sizeof(char));
  temp->model_name = (char *)malloc((strlen(model_name) + 1) * sizeof(char));
  if (trans_type != NULL)
    temp->trans_type = (char *)malloc(sizeof(trans_type)* sizeof(char));
  strcpy(temp->name, name);
  temp->type = type;
  strcpy(temp->plus_d_c, plus_d_c);
  strcpy(temp->minus_g_b, minus_g_b);
  strcpy(temp->source_e, source_e);
  strcpy(temp->body, body);
  strcpy(temp->model_name, model_name);
  temp->L = L;
  temp->W = W;
  temp->area = area;
  temp->value = value;
  if (trans_type != NULL)
    strcpy(temp->trans_type, trans_type);
  temp->i1 = i1;
  temp->i2 = i2;
  temp->td1 = td1;
  temp->tc1 = tc1;
  temp->td2 = td2;
  temp->tc2 = tc2;
  temp->per = per;
  if (size_of_pwl != 0)
  {
    temp->pwl = (double *)malloc(size_of_pwl * sizeof(double));
    memcpy(temp->pwl, pwl, size_of_pwl*sizeof(double));
  }
  else
    temp->pwl = NULL;
  temp->nextPtr = root;
  root = temp;
}
//////////////////////////////////////////////////////////////////////////////////////////////////
// This function prints our list
//////////////////////////////////////////////////////////////////////////////////////////////////
void print_list(hashtable_t *hashtable)
{
  elements *check;
  for (check = root; check != NULL; check = check->nextPtr)
  {
    if (check->type == 'v' || check->type == 'i' || check->type == 'r' || check->type == 'c' || check->type == 'l')
      //printf("%s %c %s %s %e\n", check->name, toupper(check->type), check->plus_d_c, check->minus_g_b, check->value);
      printf("%s %ld %s %ld\n", check->plus_d_c, ht_get(hashtable, check->plus_d_c), check->minus_g_b, ht_get(hashtable, check->minus_g_b));
    else if (check->type == 'd')
      //printf("%s %c %s %s %s %s %e\n", check->name, toupper(check->type), check->plus_d_c, check->minus_g_b, check->source_e, check->model_name, check->area);
      printf("%s %ld %s %ld %s %ld\n", check->plus_d_c, ht_get(hashtable, check->plus_d_c), check->minus_g_b, ht_get(hashtable, check->minus_g_b), check->source_e, ht_get(hashtable, check->source_e));
    else if (check->type == 'm')
      //printf("%s %c %s %s %s %s %s L=%e W=%e\n", check->name, toupper(check->type), check->plus_d_c, check->minus_g_b, check->source_e, check->body, check->model_name, check->L, check->W);
      printf("%s %ld %s %ld %s %ld\n", check->plus_d_c, ht_get(hashtable, check->plus_d_c), check->minus_g_b, ht_get(hashtable, check->minus_g_b), check->source_e, ht_get(hashtable, check->source_e), check->body, ht_get(hashtable, check->body));
    else if (check->type == 'q')
      //printf("%s %c %s %s %s %s\n", check->name, toupper(check->type), check->plus_d_c, check->minus_g_b, check->source_e, check->model_name);
      printf("%s %ld %s %ld %s %ld\n", check->plus_d_c, ht_get(hashtable, check->plus_d_c), check->minus_g_b, ht_get(hashtable, check->minus_g_b), check->source_e, ht_get(hashtable, check->source_e));
  }
}