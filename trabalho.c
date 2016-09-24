#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#define pi 3.14

double *tm, *tr, *res, *err;
int count = 0;

double timestamp(void)
{
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
}

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
			
                        for(int l=0; l<=k; ++l)
			{
				if((j >= 0)&&(j < n))
                                	result[i] += A[l][i] * x[j];
                        }
                }
        }
	return result;
}

double multiply_arrays(double *a, double *b, int n)
{
	double result;
	for(int i=0; i<n; ++i)
	{
		result += a[i] * b[i];
	}
	return result;
}

double euclidean_norm(double *v, unsigned int n)
{
	double norm = 0;
	for(int i=0; i<n; ++i)
		norm += pow(v[i], 2);
	return sqrt(norm);
}

double *conjugatedGradient(double **A, double *x, double *b, int n, int k)
{
	double *result, *r, *Ax, *Ar;
	double begin, end;
	double s;

	r = malloc(n*sizeof(double));
	result = malloc(n*sizeof(double));
	
	begin = timestamp();
	Ax = multiply_matrix_array(A, x, n, k);
	for(int i=0; i<n; ++i)
	{
		r[i] = b[i] - Ax[i];
	}
	end = timestamp();

	tr[count] = end - begin;
	res[count] = euclidean_norm(r, n);

	if(count > 0)
		err[count] = res[count] - res[count-1];
	else
		err[count] = 0.0;

	Ar = multiply_matrix_array(A, r, n, k);
	s = multiply_arrays(r, r, n)/multiply_arrays(r, Ar, n);
	for(int i=0; i<n; ++i)
	{
		result[i] = x[i] + (s * r[i]);
	}
	return result;
}

double min(double *v, int n)
{
	double min = v[0];
	for(int i=1; i<n; ++i)
	{
		if(v[i] < min)
			min = v[i];
	}
	return min;
}

double max(double *v, int n)
{
	double max = v[0];
	for(int i=1; i<n; ++i)
	{
		if(v[i] > max)
			max = v[i];
	}
	return max;
}

double avg(double *v, int n)
{
	double avg = 0.0;
	for(int i=0; i<n; ++i)
		avg += v[i];
	return avg/n;
}

void save_file(char *filename, double *x,  int i, int n)
{
	FILE *fp;

	fp = fopen(filename, "w+");
        fprintf (fp, "###########\n");
	fprintf (fp, "# Tempo método CG: %f %f %f\n", min(tm, n), avg(tm, n), max(tm, n));
	fprintf (fp, "# Tempo resíduo: %f %f %f", min(tr, n), avg(tr, n), max(tr,n));
	fprintf (fp, "#\n");

	fprintf (fp, "# Norma Euclidiana do Resíduo e Erro aproximado\n");
	for(int j=1; j<=i; ++j)
		fprintf(fp, "# i=%d: %f %f\n", j, res[j], err[j]);

	fprintf (fp, "###########\n");
	fprintf (fp, "%d", n);
	for(int j=0; j<n; ++j)
		fprintf(fp, "%.14g ", x[j]);

	fclose(fp);

}

int main (int argc, char *argv[])
{
	int n, k, i;
	double t = 0.0;
	double *b, *x;
	double begin, end;
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
		else
		{
			t = atof(argv[4]);
			i = n;
		}
		strcpy(output, argv[6]);
	}
	else
	{
		strcpy(output, argv[4]);
	}

	err = malloc(i*sizeof(double));
	res = malloc(i*sizeof(double));
	tm = malloc(i*sizeof(double));
	tr = malloc(i*sizeof(double));

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
	
	if (t==0.0)
	{
		while(count < i)
		{
			begin = timestamp();
			x = conjugatedGradient(A, x, b, n, k);
			end = timestamp();
			tm[count] = end-begin;
			++count;
		}
	}
	else
	{
		do
		{	
			begin = timestamp();
			x = conjugatedGradient(A, x, b, n, k);
			end = timestamp();
			tm[count] = end-begin;
			++count;
		}while((count < i)||(err[count] <= t));
	}
	save_file(output, x, i, n);

	return 0;
}

