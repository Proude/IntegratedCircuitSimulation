#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "hash_table.h"

////////////////////////////////////////////////////////////////////////
// This functions is used to read from a .txt file and insert elements
// into our list. We also check for txt correctness.
////////////////////////////////////////////////////////////////////////

int read_file(char *file_name, hashtable_t *hashtable, int *plot_list, int *checkForCholesky, char *typeOfDC, double *startValue2, double *endValue2, double *stepValue2, int *checkForDC, char *nameOfDC,int *checkForIter,int *checkForSpdIter,double *checkForItol,int *CheckForSparse,int *CheckForSpdSparce,int *CheckForIterSparce, int *CheckForSpdIterSparce, double *time_step, double *fin_time, int *numpwl,int *CheckForMethodTr,int *CheckForMethodBe)
{

  FILE *netlistFile = fopen(file_name,"r");
  char *line;
  size_t line_size = 0;
  char temp_type;
  char *temp_name;
  char *temp_node_pos_d_c;
  char *temp_node_neg_g_b;
  char *temp_source_e;
  char *temp_body;
  char *temp_model_name;
  char *temp_L_string;
  char *temp_W_string;
  double temp_area;
  double temp_L;
  double temp_W;
  double temp_value;

  int i;
  int check = 0;
  int count_words = 0, count_line = 0;
  char *result = NULL;
  char *result2 = NULL;
  char *result3 = NULL;
  char *temp_substring = NULL;
  char *linecounter = NULL;
  char *tempora = NULL;
  int temporary;

  char *sparce;
  char *spditer;
  char *option;
  char *spd;
  char *dc;
  char *plot;
  char *print;
  char *nodePrintPlot;
  char *test;
  char dot;
  char rchar;
  double startValue;
  double endValue;
  double increament;
  int helpTemp;
  int counter = 0;
  
  int SpdIterSparce =0;
  int CiterSparce =0;
  int CSparce =0;
  int CspdSparce =0;
  int Citer = 0 ; 
  int Cspditer = 0;
  double itol = 1e-3;
  
  char *trans_type = NULL;
  char *temp_rem = NULL;
  char *temp_trans = NULL;
  char *time_st = NULL;
  char *fin_ti = NULL;
  double exp_i1 = 0.0;
  double exp_i2 = 0.0;
  double exp_td1 = 0.0;
  double exp_tc1 = 0.0;
  double exp_td2 = 0.0;
  double exp_tc2 = 0.0;
  double pulse_per = 0.0;
  double *pwl;
  int num_of_pwl = 0;
  
  *checkForItol = itol;
  *checkForCholesky = 0;
  mnac = 0;
  hash_id = 1;
  while(getline(&line, &line_size, netlistFile ) != -1)
  {
    count_line++;
    if (line[0] == '\n') continue;	// If line is empty ignore.
    for(i = 0; line[i]; i++)
    {
      line[i] = tolower(line[i]);	// Non case-sensitive
    }
    
    linecounter = (char *)malloc(line_size*sizeof(char));
    temp_name = (char *)malloc(line_size*sizeof(char));
    temp_node_pos_d_c = (char *)malloc(line_size*sizeof(char));
    temp_node_neg_g_b = (char *)malloc(line_size*sizeof(char));
    temp_source_e = (char *)malloc(line_size*sizeof(char));
    temp_body = (char *)malloc(line_size*sizeof(char));
    temp_model_name = (char *)malloc(line_size*sizeof(char));
    temp_L_string = (char *)malloc(line_size*sizeof(char));
    temp_W_string = (char *)malloc(line_size*sizeof(char));
    //spyros
    spd = (char *)malloc(line_size*sizeof(char));
    sparce = (char *)malloc(line_size*sizeof(char));
    spditer = (char *)malloc(line_size*sizeof(char));
    nodePrintPlot = (char *)malloc(line_size*sizeof(char));
    print = (char *)malloc(line_size*sizeof(char));
    plot = (char *)malloc(line_size*sizeof(char));
    option = (char *)malloc(line_size*sizeof(char));
    plot = (char *)malloc(line_size*sizeof(char));
    dc = (char *)malloc(line_size*sizeof(char));
    test = (char *)malloc(line_size*sizeof(char));
    
    temp_trans = (char *)malloc(6*sizeof(char));
    time_st = (char *)malloc(10*sizeof(char));
    fin_ti = (char *)malloc(10*sizeof(char));
    pwl = (double *)malloc(0);
    
    strcpy(linecounter, line);
    result = strtok(linecounter, "\n");		//we remove '\n' char from our line 
    result2 = strtok(result, " ");
    temp_type = result[0];

      
    count_words = 0;
    while (result2 != NULL)	// Here we count the words we expect so as to check for correctness
    {
        if (temp_type == '.')
        {
            sscanf(line, "%s",test);
            if ((strcmp (test,".PLOT" )==0) || (strcmp (test,".plot" )==0) || (strcmp (test,".PRINT" )==0) || (strcmp (test,".print" )==0))
            {
                plot_count_words = count_words;
            }
        }
        count_words++;
        result2 = strtok(NULL, " ");
    }
    
    if (plot_count_words != 0 )
    {
      plot_list  = (int *)realloc(plot_list, plot_count_words * sizeof(int));
    }

    strcpy(linecounter, line);
    result = strtok(linecounter, "\n");		//we remove '\n' char from our line
    result2 = strtok(result, " ");
      
    while (result2 != NULL)
    {
        if (temp_type == '.')
        {
            sscanf(line, "%s",test);
            if (strcmp (test,".PLOT" )==0 || (strcmp (test,".plot" ))==0 || strcmp (test,".PRINT" )==0 || (strcmp (test,".print" ))==0)
            {
                if ((strcmp (result2,".PLOT" )!= 0) && ((strcmp (result2,".plot" ))!= 0) && (strcmp (result2,".PRINT" )!= 0) && ((strcmp (result2,".print" ))!= 0) )
                {
		  
		 if(result2[1] == '(')
		  {
		    temporary = (int)(result2[2]);
		    plot_list[counter] = temporary;
		    counter++;
		  }
		  else
		  {
		    result2 = ++result2;
		    helpTemp = atoi(result2);
		    plot_list[counter] = helpTemp;
		    counter++;
		  }
                }
            }
        }
        result2 = strtok(NULL, " ");
    }
    
      
      
    if (temp_type == '.') {
       sscanf(line, "%c %s",&temp_type,test);

       if ((strcmp (test,"OPTIONS" )== 0) || ((strcmp (test,"options" ))== 0) ) {
            sscanf(line, "%c %s %s %s %lf", &dot,option,spd,spditer,&itol);
	    
	    if (!strcmp(spd,"itol"))
	    {
	      *checkForItol = itol;
	      printf("to itol einai = %lf \n",itol);
	    }
            if (!strcmp(spd, "method=tr"))
	    {
	      *CheckForMethodTr =1;
	    }
	    if (!strcmp(spd, "method=be"))
	    {
	      *CheckForMethodBe =1;
	      *CheckForMethodTr = 0;
	    }
	    if (!strcmp(spd,"sparse"))
	    {
	      CSparce = 1;
	      *CheckForSparse=CSparce;
	    }
	    if (!strcmp(spd, "iter"))
	    {
	      if (!strcmp(spditer,"sparse"))
	      {
		CiterSparce=1;
		*CheckForIterSparce = CiterSparce;
		*checkForIter= Citer;
	      }
	      else
	      {
	      Citer=1;
	      *checkForIter= Citer;
	      }
	    }
	    if (!strcmp(spd, "spd"))
	    {
	      if(!strcmp(spditer,"iter"))
	      {
		sscanf(line, "%c %s %s %s %s", &dot,option,spd,spditer,sparce);
		if(!strcmp(sparce,"sparse"))
		{
		  SpdIterSparce=1;
		  *CheckForSpdIterSparce = SpdIterSparce;
		  *checkForSpdIter = Cspditer;
		}
		else
		{
		Cspditer = 1;
		*checkForSpdIter = Cspditer;
		}
	      }
	      else if(!strcmp(spditer,"sparse"))
	      {
		CspdSparce = 1;
		*CheckForSpdSparce = CspdSparce;
		*checkForCholesky = 1;
	      }
	      else
	      {
		*checkForCholesky = 1;
	      }
	      
	    }
       }
       if (strcmp (test,"DC" ) ==0 || (strcmp (test,"dc" ))==0) {
            sscanf(line, "%c %s %s %lf %lf %lf", &dot,dc,nameOfDC,&startValue,&endValue,&increament);
	    *checkForDC = 1;
	    *typeOfDC = nameOfDC[0];
	    *startValue2 = startValue;
	    printf("%s ", nameOfDC);
	    *endValue2 = endValue;
	    if (increament != 0)
	      *stepValue2 = increament;
	    else
	      *stepValue2 = 1;
       }
       if (strcmp (test,"PLOT" )==0 || (strcmp (test,"plot" ))==0) {
           sscanf(line, "%c %s %s", &rchar,plot,nodePrintPlot);
       }
       if (strcmp (test,"PRINT" )==0 || (strcmp (test,"print" ))==0) {
           sscanf(line, "%c %s %s", &rchar,print,nodePrintPlot);            
       }
       if (strcmp (test,"tran") == 0) {
	 sscanf(line, "%s %s %s", temp_trans, time_st, fin_ti);
	 *CheckForMethodTr = 1;
	 *time_step = atof(time_st);
	 *fin_time = atof(fin_ti);
       }
    }
      
    if(temp_type == 'v' || temp_type == 'i' || temp_type == 'r' || temp_type == 'c' || temp_type == 'l')
    {
      if (count_words == 4)
      {
	sscanf(line,"%c %s %s %s %lf",&temp_type, temp_name, temp_node_pos_d_c, temp_node_neg_g_b, &temp_value);
	if (strcmp(temp_node_pos_d_c, "0") == 0 || strcmp(temp_node_neg_g_b, "0") == 0)
	  check = 1;
	if (strcmp(temp_node_pos_d_c, "0") == 0)
	  ht_set(hashtable, temp_node_pos_d_c, 0);
	else
	{
	  ht_set(hashtable, temp_node_pos_d_c, hash_id);
	  if(ht_get(hashtable, temp_node_pos_d_c) == hash_id)
	    hash_id++;
	}
	if(strcmp(temp_node_neg_g_b, "0") == 0)
	  ht_set(hashtable, temp_node_neg_g_b, 0);
	else
	{
	  ht_set(hashtable, temp_node_neg_g_b, hash_id);
	  if(ht_get(hashtable, temp_node_neg_g_b) == hash_id)
	    hash_id++;
	}
	add_to_list(temp_name, temp_type, temp_node_pos_d_c, temp_node_neg_g_b, "", "", temp_value, "", -1.0, -1.0, 1.0, NULL, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, NULL, 0);
      }
      else
      {
	sscanf(line,"%c %s %s %s %lf", &temp_type, temp_name, temp_node_pos_d_c, temp_node_neg_g_b, &temp_value);
	if (strcmp(temp_node_pos_d_c, "0") == 0 || strcmp(temp_node_neg_g_b, "0") == 0)
	  check = 1;
	if (strcmp(temp_node_pos_d_c, "0") == 0)
	  ht_set(hashtable, temp_node_pos_d_c, 0);
	else
	{
	  ht_set(hashtable, temp_node_pos_d_c, hash_id);
	  if(ht_get(hashtable, temp_node_pos_d_c) == hash_id)
	    hash_id++;
	}
	if(strcmp(temp_node_neg_g_b, "0") == 0)
	  ht_set(hashtable, temp_node_neg_g_b, 0);
	else
	{
	  ht_set(hashtable, temp_node_neg_g_b, hash_id);
	  if(ht_get(hashtable, temp_node_neg_g_b) == hash_id)
	    hash_id++;
	}
	strcpy(linecounter, line);
	result = strtok(linecounter, "\n");		//we remove '\n' char from our line
	result2 = strtok(result, " ");
	while ((strcmp(result2, "exp") != 0) && (strcmp(result2, "sin") != 0) && (strcmp(result2, "pulse") != 0) && (strcmp(result2, "pwl") != 0))
	  result2 = strtok(NULL, " ");
	if (strcmp(result2, "pwl") != 0)
	{
	  trans_type = (char *)malloc(sizeof(result2)*sizeof(char));
	  strcpy(trans_type, result2);
	  result2 = strtok(NULL, " ");
	  //result3 = strtok(result2, "(");
	  exp_i1 = (double)atof(result2+1);
	  result2 = strtok(NULL, " ");
	  exp_i2 = (double)atof(result2);
	  result2 = strtok(NULL, " ");
	  exp_td1 = (double)atof(result2);
	  result2 = strtok(NULL, " ");
	  exp_tc1 = (double)atof(result2);
	  result2 = strtok(NULL, " ");
	  exp_td2 = (double)atof(result2);
	  result2 = strtok(NULL, " ");
	  exp_tc2 = (double)atof(result2);
	  result2 = strtok(NULL, " ");
	  if (strcmp(trans_type, "pulse") == 0)
	  {
	    pulse_per = (double)atof(result2);
	    add_to_list(temp_name, temp_type, temp_node_pos_d_c, temp_node_neg_g_b, "", "", temp_value, "", -1.0, -1.0, 1.0, trans_type, exp_i1, exp_i2, exp_td1, exp_tc1, exp_td2, exp_tc2, pulse_per, NULL, 0);
	  }
	  else
	  {
	    add_to_list(temp_name, temp_type, temp_node_pos_d_c, temp_node_neg_g_b, "", "", temp_value, "", -1.0, -1.0, 1.0, trans_type, exp_i1, exp_i2, exp_td1, exp_tc1, exp_td2, exp_tc2, 0.0, NULL, 0);
	  }
	  free(trans_type);
	}
	else
	{
	  trans_type = (char *)malloc(sizeof(result2)*sizeof(char));
	  strcpy(trans_type, result2);
	  result2 = strtok(NULL, " ");
	  num_of_pwl = 0;
	  while (result2 != NULL)
	  {
	    num_of_pwl += 2;
	    pwl = (double *)realloc(pwl, num_of_pwl * sizeof(double));
	    pwl[num_of_pwl - 2] = (double)atof(result2+1);
	    result2 = strtok(NULL, " ");
	    pwl[num_of_pwl - 1] = atof(result2);
	    result2 = strtok(NULL, " ");
	  }
	  *numpwl = num_of_pwl;
	  add_to_list(temp_name, temp_type, temp_node_pos_d_c, temp_node_neg_g_b, "", "", temp_value, "", -1.0, -1.0, 1.0, trans_type, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pwl, num_of_pwl);
	  free(trans_type);
	}
      }
      if(temp_type == 'v' || temp_type == 'l')
	mnac++;
    }
    else if(temp_type == 'd')
    {
      if (count_words == 4)
      {
	sscanf(line,"%c %s %s %s %s",&temp_type, temp_name, temp_node_pos_d_c, temp_node_neg_g_b, temp_model_name);
	add_to_list(temp_name, temp_type, temp_node_pos_d_c, temp_node_neg_g_b, "", "", 0, temp_model_name, -1, -1, 1, NULL, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, NULL, 0);
      }
      else if (count_words == 5)
      {
	sscanf(line,"%c %s %s %s %s %lf",&temp_type, temp_name, temp_node_pos_d_c, temp_node_neg_g_b, temp_model_name, &temp_area);
	add_to_list(temp_name, temp_type, temp_node_pos_d_c, temp_node_neg_g_b, "", "", 0, temp_model_name, -1, -1, temp_area, NULL, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, NULL, 0);
      }
      else
      {
	printf("Syntactic error in line %d\n", count_line);
	return 0;
      }
      if (strcmp(temp_node_pos_d_c, "0") == 0|| strcmp(temp_node_neg_g_b, "0") == 0)
	  check = 1;
      if (strcmp(temp_node_pos_d_c, "0") == 0)
	  ht_set(hashtable, temp_node_pos_d_c, 0);
	else
	{
	  ht_set(hashtable, temp_node_pos_d_c, hash_id);
	    if(ht_get(hashtable, temp_node_pos_d_c) == hash_id)
	    hash_id++;
	}
	if(strcmp(temp_node_neg_g_b, "0") == 0)
	  ht_set(hashtable, temp_node_neg_g_b, 0);
	else
	{
	  ht_set(hashtable, temp_node_neg_g_b, hash_id);
	  if(ht_get(hashtable, temp_node_neg_g_b) == hash_id)
	    hash_id++;
	}
    }
    else if(temp_type == 'm')
    {
      if (count_words == 8)
      {
	sscanf(line,"%c %s %s %s %s %s %s %s %s",&temp_type, temp_name, temp_node_pos_d_c, temp_node_neg_g_b, temp_source_e, temp_body, temp_model_name, temp_L_string, temp_W_string);
	temp_substring = strtok(temp_L_string, "=");	// Next instructions are used to convert (e.g) "L=10e-2" to "10e-2" so as we can use it as double
	temp_substring = strtok(NULL, "=");
	temp_L = atof(temp_substring);
	temp_substring = strtok(temp_W_string, "=");
	temp_substring = strtok(NULL, "=");
	temp_W = atof(temp_substring);
	if (strcmp(temp_node_pos_d_c, "0") == 0 || strcmp(temp_node_neg_g_b, "0") == 0 || strcmp(temp_source_e, "0") == 0 || strcmp(temp_body, "0") == 0)
	  check = 1;
	if (strcmp(temp_node_pos_d_c, "0") == 0)
	  ht_set(hashtable, temp_node_pos_d_c, 0);
	else
	{
	  ht_set(hashtable, temp_node_pos_d_c, hash_id);
	  if(ht_get(hashtable, temp_node_pos_d_c) == hash_id)
	    hash_id++;
	}
	if(strcmp(temp_node_neg_g_b, "0") == 0)
	  ht_set(hashtable, temp_node_neg_g_b, 0);
	else
	{
	  ht_set(hashtable, temp_node_neg_g_b, hash_id);
	  if(ht_get(hashtable, temp_node_neg_g_b) == hash_id)
	    hash_id++;
	}
	if(strcmp(temp_source_e, "0") == 0)
	  ht_set(hashtable, temp_source_e, 0);
	else
	{
	  ht_set(hashtable, temp_source_e, hash_id);
	  if(ht_get(hashtable, temp_source_e) == hash_id)
	    hash_id++;
	}
	if(strcmp(temp_body, "0") == 0)
	  ht_set(hashtable, temp_body, 0);
	else
	{
	  ht_set(hashtable, temp_body, hash_id);
	  if(ht_get(hashtable, temp_body) == hash_id)
	    hash_id++;
	}
	add_to_list(temp_name, temp_type, temp_node_pos_d_c, temp_node_neg_g_b, temp_source_e, temp_body, 0, temp_model_name, temp_L, temp_W, 1, NULL, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, NULL, 0);
      }
      else
      {
	printf("Syntactic error in line %d\n", count_line);
	return 0;
      }
    }
    else if(temp_type == 'q')
    {
      if (count_words == 5)
      {
	sscanf(line,"%c %s %s %s %s %s",&temp_type, temp_name, temp_node_pos_d_c, temp_node_neg_g_b, temp_source_e, temp_model_name);
	add_to_list(temp_name, temp_type, temp_node_pos_d_c, temp_node_neg_g_b, temp_source_e, "", 0, temp_model_name, -1, -1, 1, NULL, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, NULL, 0);
      }
      else if (count_words == 6)
      {
	sscanf(line,"%c %s %s %s %s %s %lf",&temp_type, temp_name, temp_node_pos_d_c, temp_node_neg_g_b, temp_source_e, temp_model_name, &temp_area);
	add_to_list(temp_name, temp_type, temp_node_pos_d_c, temp_node_neg_g_b, temp_source_e, "", 0, temp_model_name, -1, -1, temp_area, NULL, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, NULL, 0);
      }
      else
      {
	printf("Syntactic error in line %d\n", count_line);
	return 0;
      }
      if (strcmp(temp_node_pos_d_c, "0") == 0 || strcmp(temp_node_neg_g_b, "0") == 0 || strcmp(temp_source_e, "0") == 0)
	  check = 1;
      if (strcmp(temp_node_pos_d_c, "0") == 0)
	  ht_set(hashtable, temp_node_pos_d_c, 0);
	else
	{
	  ht_set(hashtable, temp_node_pos_d_c, hash_id);
	  if(ht_get(hashtable, temp_node_pos_d_c) == hash_id)
	    hash_id++;
	}
	if(strcmp(temp_node_neg_g_b, "0") == 0)
	  ht_set(hashtable, temp_node_neg_g_b, 0);
	else
	{
	  ht_set(hashtable, temp_node_neg_g_b, hash_id);
	  if(ht_get(hashtable, temp_node_neg_g_b) == hash_id)
	    hash_id++;
	}
	if(strcmp(temp_source_e, "0") == 0)
	  ht_set(hashtable, temp_source_e, 0);
	else
	{
	  ht_set(hashtable, temp_source_e, hash_id);
	  if(ht_get(hashtable, temp_source_e) == hash_id)
	    hash_id++;
	}
    }
    else if(temp_type == '*')
      continue;
    
    free(temp_L_string);
    free(linecounter);
    free(temp_name);
    free(temp_node_pos_d_c);
    free(temp_node_neg_g_b);
    free(temp_source_e);
    free(temp_body);
    free(temp_model_name);
    free(temp_W_string);

  }
  hash_id = hash_id - 1;
  fclose(netlistFile);
  return check;
}

void removeChar(char *str, char garbage) {

    char *src, *dst;
    for (src = dst = str; *src != '\0'; src++) {
        *dst = *src;
        if (*dst != garbage) dst++;
    }
    *dst = '\0';
}