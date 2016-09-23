#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#define pi 3.14

/***********************
 * N: tamanho do sistema linear
 * k: numero da diagonal, 0 = diagonal principal, 1 = acima/abaixo da diagonal, 2 = ...
 * kMax: numero de bandas do sistema linear
 * diag: vetor para armazenar os valores da diagonal. Deve ser alocado por quem chama a função.
 ***********************/
int generateRandomDiagonal( unsigned int N, unsigned int k, unsigned int kMax, double *diag )
{
  if ( !diag || N < 3 || kMax > N/2 || k < 0 || k > kMax )
    return (-1);

  /* garante valor dominante para diagonal principal */
  double fator = (k == 0) ? ((double)(kMax-1)) : (0.0);

  double invRandMax = 1.0 / (double)RAND_MAX;

  for (int i=0; i < N-k; ++i)
  {
    diag[i] = fator + (double)rand() * invRandMax;
  }

  return (0);
}

double f(double x)
{
	return((4*pow(pi, 2))*(sin(2*pi*x)+sin(2*pi*(pi-x))));
}


double *generateResultVector(unsigned int N)
{
	double *result = malloc(N*sizeof(double));
	for (int i=0; i < N; ++i)
	{
		result[i] = f(i*pi/N);
	}
	return result;
}

double *multiply_matrix_array(double **A, double *x, int n, int k)
{
	double *result = malloc(n*sizeof(double));
	for(int i=0; i<n; ++i)
        {
        	for(int j=i-k; j<=i+k; ++j)
                {
                        if (j >= 0)&&(j < n)
                        {
                                result[i] += A[i+abs(j)][i] * x[j];
                        }
                }
        }
	return result
}

double multiply_arrays(double *a, double* b, int n)
{
	double result;
	for(int i=0; i<n; ++i)
	{
		result += a[i] * b[i];
	}
	return result;
}

double *conjugatedGradient(double **A, double *x, double *b, int n, int k)
{
	double *result, *r, *Ax, *Ar;
	double s;
	r = malloc(n*sizeof(double));
	s = malloc(n*sizeof(double));

	Ax = multiply_matrix_array(A, x);

	for(int i=0; i<n; ++i)
	{
		r[i] = b[i] - Ax[i];
	}

	Ar = multiply_matrix_array(A, r);
	s = multiply_arrays(r, r)/multiply_arrays(r, Ar);

	for(int i=0; i<n; ++i)
	{
		result[i] = x[i] + (s * r[i]);
	}

	return result;
}

int main (int argc, char *argv[])
{
	int n, k, i;
	double t = 0.0;
	double *b, *x;
	char *output = malloc(256*sizeof(char));

 	n = atoi(argv[1]);
	k = atoi(argv[2]);

	if (argc == 9)
	{ //todos os parametros
		i = atoi(argv[4]);
		t = atof(argv[6]);
		strcpy (output, argv[8]);
	}
	else if (argc == 7)
	{
		if (strcmp(argv[3],"-i")==0)
			 i = atoi(argv[4]);
		else t = atof(argv[4]);
		strcpy(output, argv[6]);
	}
	else
	{
		strcpy(output, argv[4]);
	}


	x = malloc(n*sizeof(double));
	for (int i=0; i<n; ++i)
		x[i] = 0.0;

	b = generateResultVector(n);

	double ** A = malloc ((k+1)*sizeof(double *));
	srand(20162);
	for (int count=0; count<=k; ++count)
	{
		A[count] = (double *)malloc(n*sizeof(double));
		generateRandomDiagonal (n, count, k, A[count]);
	}
	
	for(int i = 0; i < n; ++i)
	{
		printf("b[%d] = %f\n",i,  b[i]);
	}

	for(int i=0; i<(k+1); ++i)
	{
		for(int j=0; j<n; ++j)
		{
			printf("A[%d][%d] = %f", i, j, A[i][j]);
		}
		printf("\n");
	}

	do
	{
		x = conjugatedGradient(A, x, b, n, k);
	}while()||()
	return 0;
}
