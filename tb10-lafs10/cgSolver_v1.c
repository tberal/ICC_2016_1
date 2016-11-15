#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#define M_PI 3.14159265358979323846

double *tm, *tr, *res, *err;
int count = 0;

// pega o tempo atual
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

// f(x) = 4 * pi² * (sin(2 * pi * x) + sin(2 * pi * (pi - x)))
double f(double x)
{
	return((4*pow(M_PI, 2))*(sin(2*M_PI*x)+sin(2*M_PI*(M_PI-x))));
}

// gera um vetor com os valores de f(x)
double *generateResultVector(unsigned int N)
{
	double *result = malloc(N*sizeof(double));
	for (int i=0; i < N; ++i)
	{
		result[i] = f(i*M_PI/N);
	}
	return result;
}

/*
	Função para multiplicar elementos de uma matriz NxK+1 por um vetor de N elementos	

*/
double *multiply_matrix_array(double **A, double *x, int n, int k)
{
	int line, offset;
	double *result = malloc(n*sizeof(double));
  for(int i=0; i<n; ++i)
    result[i] = 0.0;

	for(int i=0; i<n; ++i)
        {
        	for(int j=i-k; j<=i+k; ++j)
                {
			if ((j>=0)&&(j<n))
			{
				line = abs(i - j);
				offset = j + k;
				result[i] += A[line][i] * x[j];
			}	
                }
        }
	return result;
}

// multiplicação de vetores
double multiply_arrays(double *a, double *b, int n)
{
	double result = 0.0;
  for(int i=0; i<n; ++i)
	{
		result += a[i] * b[i];
	}
	return result;
}

// calcula a norma euclideana
double euclidean_norm(double *v, unsigned int n)
{
	double norm = 0.0;
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
		err[count] = fabs(res[count] - res[count-1]);
	else
		err[count] = fabs(res[count]);

	Ar = multiply_matrix_array(A, r, n, k);
	s = multiply_arrays(r, r, n)/multiply_arrays(r, Ar, n);
	for(int i=0; i<n; ++i)
	{
		result[i] = x[i] + (s * r[i]);
	}
  free (Ax);
  free (Ar);
  free (r);
	return result;
}

// encontra o menor valor de um vetor
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

// encontra o maior valor de um vetor
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

// calcula a média dos valor de um vetor v de tamanho n
double avg(double *v, int n)
{
	double avg = 0.0;
	for(int i=0; i<n; ++i)
		avg += v[i];
	return avg/n;
}

//geração do arquivo
void save_file(char *filename, double *x, int n)
{
	FILE *fp;

	fp = fopen(filename, "w+");
  fprintf (fp, "###########\n");
	fprintf (fp, "# Tempo método CG: %f %f %f\n", min(tm, count), avg(tm, count), max(tm, count));
	fprintf (fp, "# Tempo resíduo: %f %f %f\n", min(tr, count), avg(tr, count), max(tr,count));
	fprintf (fp, "#\n");

	fprintf (fp, "# Norma Euclidiana do Resíduo e Erro aproximado\n");
	for(int j=0; j<count; ++j)
		fprintf(fp, "# i=%d: %f %f\n", j+1, res[j], err[j]);

	fprintf (fp, "###########\n");
	fprintf (fp, "%d\n", n);
	for(int j=0; j<n; ++j)
		fprintf(fp, "%.14g ", x[j]);

	fclose(fp);

}

//programa principal
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
	else if(argc == 5)
	{
		strcpy(output, argv[4]);
		i = n;
	}
	else
	{
		printf("Parâmetros insuficientes\n");
		exit(1);
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
	for (int i=0; i<=k; ++i)
	{
		A[i] = (double *)malloc(n*sizeof(double));
		generateRandomDiagonal (n, i, k, A[i]);
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
			if(err[count] <= t)
			{
				++count;
				break;
			}
			++count;
		}while(count < i);
	}
	save_file(output, x, n);
  for (int i=0; i < (k+1); ++i)
  {
    free (A[i]);
  }
  free (A);
  free (x);
  free (b);
  free (err);
  free (res);
  free (tm);
  free (tr);
	return 0;
}

