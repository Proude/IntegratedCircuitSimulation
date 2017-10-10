#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "hash_table.h"
#include "csparse.h"

int main(int argc, char *argv[])
{
  double *MNA_matrix = NULL;
  double *MNA_right_vector = NULL;
  double *MNA_matrix2 = NULL;
  double *MNA_right_vector2 = NULL;
  double *MNA_right_vector_sparse = NULL;
  double *G_matrix = NULL;
  double *G_matrix2 = NULL;
  double *C_matrix = NULL;
  double *C_matrix2 = NULL;
  double *trans_right_matrix = NULL;
  double *trans_right_matrix2 = NULL;
  double *xBiCG = NULL;
  double *xCG = NULL;
  double *xSparse = NULL;
  double *xLU = NULL;
  double *x_final;
  double *BE_table = NULL;
  double *BE_right_vector = NULL;
  double *etk = NULL;
  double *xti = NULL;
  int *plot_list = malloc(0*sizeof(int));
  int checkForCholesky = 0;
  int checkForIter = 0;
  int checkForSpdIter = 0;
  int checkForDC = 0;
  int CheckForSparse = 0;
  int CheckForSpdSparce = 0;
  int CheckForIterSparce = 0;
  int CheckForSpdIterSparce = 0;
  int CheckForMethodTr = 0;
  int CheckForMethodBe = 0;
  int i, j;
  int res;
  int nonZeros;
  double time_step;
  double fin_time;
  double temp_type;
  int h_val = 0;
  int numpwl;
  char typeOfDC;
  char *nameOfDC = malloc(256*sizeof(char));
  double startValue, endValue, stepValue,checkForItol;
  double itol,iter;
  gsl_matrix_view A;
  gsl_vector_view b;
  gsl_vector *x;
  gsl_vector *xtk;
  cs *Sparse_MNA_matrix;
  cs *Sparse_MNA_matrix_compress;
  hashtable_t *hashtable = ht_create( 65536 );
  
  list_init();
  if (argv[1] == NULL)
  {
    printf("No input file\n");
    return 0;
  }
  res = read_file(argv[1], hashtable, plot_list, &checkForCholesky, &typeOfDC, &startValue, &endValue, &stepValue, &checkForDC,
		  nameOfDC,&checkForIter,&checkForSpdIter,&checkForItol,&CheckForSparse,&CheckForSpdSparce,&CheckForIterSparce, &CheckForSpdIterSparce, &time_step, &fin_time, &numpwl,&CheckForMethodTr,&CheckForMethodBe);
  for (temp_type = 0; temp_type <= fin_time; temp_type += time_step)
    h_val++;
  if (res == 0)
  {
    printf("fail");
    return 0;
  }
  if (CheckForMethodBe || CheckForMethodTr)
  {
    if (CheckForSparse || CheckForSpdSparce || CheckForIterSparce || CheckForSpdIterSparce)
    {
      /*xSparse = malloc((mnac + hash_id)*sizeof(double));
  
      nonZeros = Calculate_nonZeros(hashtable);
      MNA_right_vector_sparse = calloc((mnac+hash_id),sizeof(double));
      Sparse_MNA_matrix = cs_spalloc(mnac+hash_id,mnac+hash_id,nonZeros,1,1);
  
      create_Sparse_matrix(hashtable, Sparse_MNA_matrix,MNA_right_vector_sparse, nonZeros);
      Sparse_MNA_matrix_compress = cs_compress(Sparse_MNA_matrix);
      cs_dupl(Sparse_MNA_matrix_compress);
      
      if (CheckForSparse)
      {
	printf("LU method Sparse\n");
	LUSparseDecomposition(Sparse_MNA_matrix_compress, xSparse, MNA_right_vector_sparse);
      }
      if (CheckForSpdSparce)
      {
	printf("Cholesky method Sparse\n");
	CholeskySparseDecomposition(Sparse_MNA_matrix_compress, xSparse, MNA_right_vector_sparse);
      }
      if (CheckForSpdIterSparce)
      {
	printf("CG Method Sparse\n");
	CG_Method_Sparse(Sparse_MNA_matrix_compress, MNA_right_vector_sparse, xSparse, itol, mnac + hash_id);
      }
      if (CheckForIterSparce)
      {
	printf("BiCG Method Sparse\n");
	BiCG_Sparse(Sparse_MNA_matrix_compress, MNA_right_vector_sparse, xSparse, itol);
      }*/
    }
    else
    {
      G_matrix = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
      C_matrix = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
      trans_right_matrix = (double *)calloc((mnac + hash_id) * (h_val + 1), sizeof(double));
      xLU = (double *)calloc((mnac + hash_id), sizeof(double));

      create_G_matrix(hashtable, G_matrix, C_matrix, trans_right_matrix, time_step, fin_time, h_val + 1, numpwl);

      if (checkForCholesky == 0 && checkForIter == 0 && checkForSpdIter == 0)
      {
	G_matrix2 = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
	trans_right_matrix2 = (double *)calloc((mnac + hash_id) * (h_val + 1), sizeof(double));
	for (i = 0; i < (mnac + hash_id) * (mnac + hash_id); i++)
	  G_matrix2[i] = G_matrix[i];
	for (i = 0; i < (mnac + hash_id) * (h_val + 1); i++)
	  trans_right_matrix2[i] = trans_right_matrix[i]; 
	
	init_starting_value_LU(G_matrix2, trans_right_matrix2, xLU, time_step, fin_time);
	
	x_final = (double *)calloc((mnac + hash_id) * (h_val + 1), sizeof(double));
	for(i = 0; i < mnac + hash_id; i++)
	  x_final[i * (h_val + 1)] = xLU[i];
	if (CheckForMethodTr)
	{
	  trapezoidal_LU(G_matrix, C_matrix, x_final, h_val, trans_right_matrix);
	}
	if (CheckForMethodBe)
	{
	  BE_table = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
	  BE_right_vector = (double *)calloc((mnac + hash_id),sizeof(double));
	  etk = (double*)calloc(mnac+hash_id,sizeof(double));
	  
	  for(i=1;i<=h_val;i++){
	    for(j=0;j<mnac+hash_id;j++){
	      etk[j] = trans_right_matrix[j*(h_val+1) + i];
	    }
	    BE_method(BE_table, BE_right_vector, G_matrix, C_matrix, xLU, etk, h_val);
	    init(&A, &b, BE_table, BE_right_vector);
	    xtk = gsl_vector_alloc(mnac+hash_id);
	    LUDecomposition(A, b, xtk);
	    printf("%d step of BE with LU\n",i);
	    gsl_vector_fprintf(stdout,xtk,"%lf");
	    for (j = 0; j < mnac+hash_id; j++){
	      xLU[j] = xtk->data[j];
	    }
	  }
	}
      }
      if (checkForCholesky == 1)
      {
	G_matrix2 = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
	trans_right_matrix2 = (double *)calloc((mnac + hash_id) * (h_val + 1), sizeof(double));
	for (i = 0; i < (mnac + hash_id) * (mnac + hash_id); i++)
	  G_matrix2[i] = G_matrix[i];
	for (i = 0; i < (mnac + hash_id) * (h_val + 1); i++)
	  trans_right_matrix2[i] = trans_right_matrix[i];
	init_starting_value_Cholesky(G_matrix2, trans_right_matrix2, xLU, time_step, fin_time);
	x_final = (double *)calloc((mnac + hash_id) * (h_val + 1), sizeof(double));
	for(i = 0; i < mnac + hash_id; i++)
	  x_final[i * (h_val + 1)] = xLU[i];
	if (CheckForMethodTr)
	{
	  trapezoidal_Cholesky(G_matrix, C_matrix, x_final, h_val, trans_right_matrix);
	}
	if (CheckForMethodBe)
	{
	  BE_table = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
	  BE_right_vector = (double *)calloc((mnac + hash_id),sizeof(double));
	  etk = (double*)calloc(mnac+hash_id,sizeof(double));
	  
	  for(i=1;i<=h_val;i++){
	    for(j=0;j<mnac+hash_id;j++){
	      etk[j] = trans_right_matrix[j*(h_val+1) + i];
	    }
	    BE_method(BE_table, BE_right_vector, G_matrix, C_matrix, xLU, etk, h_val);
	    Cholesky_init(&A, &b, BE_table, BE_right_vector);
	    xtk = gsl_vector_alloc(mnac+hash_id);
	    CholeskyDecomposition(A, b, xtk);
	    printf("%d step of BE with Cholesky\n",i);
	    gsl_vector_fprintf(stdout,xtk,"%lf");
	    for (j = 0; j < mnac+hash_id; j++){
	      xLU[j] = xtk->data[j];
	    }
	  }
	}
      }
      if (checkForSpdIter == 1)
      {
	/*G_matrix2 = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
	C_matrix2 = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
	
	for (i = 0; i < mnac + hash_id; i++)
	{
	  G_matrix2[i] = G_matrix[i];
	  C_matrix2[i] = C_matrix[i];
	}*/
	printf("CG method\n");
	xCG = malloc((mnac + hash_id)*sizeof(double));
	CG_Method(G_matrix, trans_right_matrix, xCG, itol, mnac+hash_id);
	x_final = (double *)calloc((mnac + hash_id) * (h_val + 1), sizeof(double));
	for(i = 0; i < mnac + hash_id; i++)
	  x_final[i * (h_val + 1)] = xLU[i];
	if (CheckForMethodTr)
	{
	  trapezoidal_CG(G_matrix, C_matrix, x_final, h_val, trans_right_matrix, itol);
	}
	if (CheckForMethodBe)
	{
	  BE_table = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
	  BE_right_vector = (double *)calloc((mnac + hash_id),sizeof(double));
	  etk = (double*)calloc(mnac+hash_id,sizeof(double));
	  xti = (double *)calloc(mnac + hash_id, sizeof(double));
	  
	  for(i=1;i<=h_val;i++){
	    for(j=0;j<mnac+hash_id;j++){
	      etk[j] = trans_right_matrix[j*(h_val+1) + i];
	    }
	    BE_method(BE_table, BE_right_vector, G_matrix, C_matrix, xCG, etk, h_val);
	    printf("%d step of BE with CG\n",i);
	    CG_Method(BE_table, BE_right_vector, xti, itol, mnac + hash_id);
	    for (j = 0; j < mnac+hash_id; j++){
	      xCG[j] = xti[j];
	    }
	  }
	  
	}
      }
      if (checkForIter == 1)
      {
	/*G_matrix2 = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
	C_matrix2 = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
	
	for (i = 0; i < mnac + hash_id; i++)
	{
	  G_matrix2[i] = G_matrix[i];
	  C_matrix2[i] = C_matrix[i];
	}
	*/
	printf("BiCG method\n");
	xBiCG = malloc((mnac + hash_id)*sizeof(double));
	BiCG(G_matrix, trans_right_matrix, xBiCG, itol);
	x_final = (double *)calloc((mnac + hash_id) * (h_val + 1), sizeof(double));
	for(i = 0; i < mnac + hash_id; i++)
	  x_final[i * (h_val + 1)] = xLU[i];
	if (CheckForMethodTr)
	{
	  trapezoidal_BiCG(G_matrix, C_matrix, x_final, h_val, trans_right_matrix, itol);
	}
	if (CheckForMethodBe)
	{
	  BE_table = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
	  BE_right_vector = (double *)calloc((mnac + hash_id),sizeof(double));
	  etk = (double*)calloc(mnac+hash_id,sizeof(double));
	  xti = (double *)calloc(mnac + hash_id, sizeof(double));
	  
	  for(i=1;i<=h_val;i++){
	    for(j=0;j<mnac+hash_id;j++){
	      etk[j] = trans_right_matrix[j*(h_val+1) + i];
	    }
	    BE_method(BE_table, BE_right_vector, G_matrix, C_matrix, xBiCG, etk, h_val);
	    printf("%d step of BE with BiCG\n",i);
	    BiCG(BE_table, BE_right_vector, xti, itol);
	    for (j = 0; j < mnac+hash_id; j++){
	      xBiCG[j] = xti[j];
	    }
	  }
	}
      }
    }
  }
  else
  {
    if (CheckForSparse || CheckForSpdSparce || CheckForIterSparce || CheckForSpdIterSparce)
    {
      xSparse = malloc((mnac + hash_id)*sizeof(double));
  
      nonZeros = Calculate_nonZeros(hashtable);
      MNA_right_vector_sparse = calloc((mnac+hash_id),sizeof(double));
      Sparse_MNA_matrix = cs_spalloc(mnac+hash_id,mnac+hash_id,nonZeros,1,1);
  
      create_Sparse_matrix(hashtable, Sparse_MNA_matrix,MNA_right_vector_sparse, nonZeros);
      Sparse_MNA_matrix_compress = cs_compress(Sparse_MNA_matrix);
      cs_dupl(Sparse_MNA_matrix_compress);
    
      if (CheckForSparse)
      {
	printf("LU method Sparse\n");
	LUSparseDecomposition(Sparse_MNA_matrix_compress, xSparse, MNA_right_vector_sparse);
      }
      if (CheckForSpdSparce)
      {
	printf("Cholesky method Sparse\n");
	CholeskySparseDecomposition(Sparse_MNA_matrix_compress, xSparse, MNA_right_vector_sparse);
      }
      if (CheckForSpdIterSparce)
      {
	printf("CG Method Sparse\n");
	CG_Method_Sparse(Sparse_MNA_matrix_compress, MNA_right_vector_sparse, xSparse, itol, mnac + hash_id);
      }
      if (CheckForIterSparce)
      {
	printf("BiCG Method Sparse\n");
	BiCG_Sparse(Sparse_MNA_matrix_compress, MNA_right_vector_sparse, xSparse, itol);
      }
      cs_spfree(Sparse_MNA_matrix);
    }
  /*else if (CheckForMethodBe || CheckForMethodTr)
  {
    
    G_matrix = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
    C_matrix = (double *)calloc((mnac + hash_id) * (mnac + hash_id), sizeof(double));
    
    trans_right_matrix = (double *)calloc((mnac + hash_id) * (h_val + 1), sizeof(double));
    create_G_matrix(hashtable, G_matrix, C_matrix, trans_right_matrix, time_step, fin_time, h_val + 1, numpwl);
    if (CheckForMethodBe == 1)
    {
    }
    if (CheckForMethodTr == 1)
    {
      printf("check\n");
      init_starting_value(G_matrix, trans_right_matrix, time_step, fin_time);
    }
  }*/
    else
    {
      MNA_matrix = calloc((mnac+hash_id)*(mnac+hash_id), sizeof(double));
      MNA_right_vector = calloc((mnac+hash_id), sizeof(double));
      MNA_matrix2 = calloc((mnac+hash_id)*(mnac+hash_id), sizeof(double));
      MNA_right_vector2 = calloc((mnac+hash_id), sizeof(double));
    
      create_MNA_matrix(hashtable, MNA_matrix, MNA_right_vector);
  
      for (i = 0; i < (mnac + hash_id)*(mnac + hash_id); i++)
	MNA_matrix2[i] = MNA_matrix[i];
      for (i = 0; i < mnac + hash_id; i++)
      {
	printf("%lf\n", MNA_right_vector[i]);
	MNA_right_vector2[i] = MNA_right_vector[i];
      }  
      FILE *outputfile = fopen("output.txt", "w");
      if (checkForCholesky == 0 && checkForIter == 0 && checkForSpdIter == 0)
      {
	printf("LU method\n");
	init(&A, &b, MNA_matrix, MNA_right_vector);
	x = gsl_vector_alloc(mnac+hash_id);
	LUDecomposition(A, b, x);
	gsl_vector_fprintf(stdout, x, "%lf");
      }
      if (checkForCholesky == 1)
      {
	printf("cholesky method\n");
	Cholesky_init(&A, &b, MNA_matrix, MNA_right_vector);
	x = gsl_vector_alloc(mnac+hash_id);
	CholeskyDecomposition(A, b, x);
	gsl_vector_fprintf(stdout, x, "%lf");
      }
      FILE *files[plot_count_words];    
      char tempstr[50];
      for(i=0; i<plot_count_words;i++){
	sprintf(tempstr,"%d",plot_list[i]);
	files[i] = fopen(tempstr,"w");
      }
      if (checkForDC == 1)
      {
	for (iter = startValue; iter <= endValue; iter+=stepValue)
	{
	  DC_init(&A, &b, MNA_matrix, MNA_right_vector, typeOfDC, nameOfDC+1, i, hashtable);
	  x = gsl_vector_alloc(mnac+hash_id);
	  DC_simulation(A, b, x);
	  for(j=0; j<plot_count_words;j++){
	    fprintf(files[j],"%lf \n",gsl_vector_get(x,plot_list[j]-1));
	  }
	  gsl_vector_fprintf(outputfile, x, "%lf");
	}
      }
      if (checkForSpdIter == 1)
      {
	printf("CG_method\n");
	xCG = malloc((mnac + hash_id)*sizeof(double));
	CG_Method(MNA_matrix2, MNA_right_vector2, xCG, itol, mnac+hash_id);
	for (i = 0; i < mnac+hash_id; i++)
	  printf("%lf\n", xCG[i]);
      }
      if (checkForIter == 1)
      {
	printf("biCG method\n");
	xBiCG = malloc((mnac + hash_id)*sizeof(double));
	BiCG(MNA_matrix2, MNA_right_vector2, xBiCG, itol);
	for (i = 0; i < mnac + hash_id; i++)
	  printf("%lf\n", xBiCG[i]);
      }
    //free(MNA_matrix);
    //free(MNA_right_vector);
    }
  }
  //CholeskySparseDecomposition(Sparse_MNA_matrix_compress, xSparse, MNA_right_vector_sparse);
  //CG_Method_Sparse(Sparse_MNA_matrix_compress, MNA_right_vector_sparse, xSparse, itol, mnac+hash_id);
  //BiCG_Sparse(Sparse_MNA_matrix_compress, MNA_right_vector_sparse, xSparse, itol);
  //cs_spfree(Sparse_MNA_matrix);
}

