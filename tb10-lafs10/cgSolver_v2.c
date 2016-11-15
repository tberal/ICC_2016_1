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
double *generateResultVector(unsigned int N, int size)
{
	double *result = malloc(size*sizeof(double));
	for (int i=0; i < N; ++i)
	{
		result[i] = f(i*M_PI/N);
	}
        for (int i=N; i< size; ++i)
            result[i] = 0.0;
	return result;
}

/*
	Função para multiplicar elementos de uma matriz NxK+1 por um vetor de N elementos

*/
double *multiply_matrix_array(double *A, double *x, int k, int size)
{
	int line, offset;
	double *result = malloc(size*sizeof(double));
        double *test = malloc(size*sizeof(double));

        for(int i=0; i<size; ++i)
            test[i] = 0.0;

  	for(int i=0; i<size; ++i)
    	    result[i] = 0.0;

	// Código original para acesso da matriz: Ineficiente pois acessa
	// elementos da matriz coluna por coluna, assim não aproveitando bem a cache
	for(int i=0; i<size; ++i)
        {
        	for(int j=i-k; j<=i+k; ++j)
                {
			if ((j>=0)&&(j<size))
			{
				line = abs(i - j);
				test[i] += A[line*size+j] * x[j];
			}
                }
        }

	// Otimização da leitura da matriz: acessa elementos da matriz por linha
	// ao invés de coluna, afim de melhor aproveitar dados em cache
	// Otimização de pipeline: Lasso desenrolado para melhor aproveitar
	// a paralelização do processador
	for(int i=0; i<=k; ++i)
	{
		for(int j=0; j<size; j++)
		{
			if(i==0)
			{
				result[j] += A[i*size+j] * x[j];
			}
			else
			{
				if(j-i >= 0)
					result[j] += A[i*size+j-i]* x[j-i];
				if(j+i < size)
					result[j] += A[i*size+j+i] * x[j+i];
			}
		}
	}
        for(int i=0; i<size; ++i)
            printf("test[%d]: %f result[%d]: %f\n", i, test[i], i, result[i]);

	return result;
}

// multiplicação de vetores
double multiply_arrays(double *a, double *b, int size)
{
	double result = 0.0;
        // Otimização: Loop desenrolado para melhor aproveitar o
        // paralelismo do processador
        for(int i=0; i<size; i+=4)
	{
		result += a[i] * b[i];
		result += a[i+1] * b[i+1];
		result += a[i+2] * b[i+2];
		result += a[i+3] * b[i+3];
	}
	return result;
}

// calcula a norma euclideana
double euclidean_norm(double *v, unsigned int size)
{
	double norm = 0.0;
	for(int i=0; i<size; i+=4)
        {
		norm += v[i] * v[i];
		norm += v[i+1] * v[i+1];
		norm += v[i+2] * v[i+2];
		norm += v[i+3] * v[i+3];
        }
	return sqrt(norm);
}

double *conjugatedGradient(double *A, double *x, double *b, int k, int size)
{
	double *result, *r, *Ax, *Ar;
	double begin, end;
	double s;

	r = malloc(size*sizeof(double));
	result = malloc(size*sizeof(double));

	begin = timestamp();
	Ax = multiply_matrix_array(A, x, k, size);
	for(int i=0; i<size; i+=4)
	{
		r[i] = b[i] - Ax[i];
		r[i+1] = b[i+1] - Ax[i+1];
		r[i+2] = b[i+2] - Ax[i+2];
		r[i+3] = b[i+3] - Ax[i+3];
	}
	end = timestamp();

	tr[count] = end - begin;
	res[count] = euclidean_norm(r, size);

	if(count > 0)
		err[count] = fabs(res[count] - res[count-1]);
	else
		err[count] = fabs(res[count]);

	Ar = multiply_matrix_array(A, r, k, size);
	s = multiply_arrays(r, r, size)/multiply_arrays(r, Ar, size);
	for(int i=0; i<size; i+=4)
	{
		result[i] = x[i] + (s * r[i]);
		result[i+1] = x[i+1] + (s * r[i+1]);
		result[i+2] = x[i+2] + (s * r[i+2]);
		result[i+3] = x[i+3] + (s * r[i+3]);
	}
	free (Ax);
	free (Ar);
	free (r);
	Ax = NULL;
	Ar =NULL;
	r = NULL;

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

        int size;
        if (n%5 !=0)
            size = n+((n/5+1)*5-n);
        else
            size = n;

	x = malloc(size*sizeof(double));
	for (int i=0; i<size; ++i)
		x[i] = 0.0;

	b = generateResultVector(n, size);

	double * A = malloc((k+1)*size*sizeof(double));
	double * aux = malloc(n*sizeof(double));

	srand(20162);

	for (int i=0; i<=k; ++i)
	{
		generateRandomDiagonal (n, i, k, aux);
		for (int j=0; j<n; ++j)
		{
			A[i*size+j] = aux[j];
			aux[j] = 0.0;
		}
                for (int j=n; j<size; ++j)
                {
                    A[i*size+j] = 0.0;
                }
	}

        free(aux);

	if (t==0.0)
	{
		while(count < i)
		{
			begin = timestamp();
			x = conjugatedGradient(A, x, b, k, size);
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
			x = conjugatedGradient(A, x, b, k, size);
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

        free (A);
        free (x);
        free (b);
        free (err);
        free (res);
        free (tm);
        free (tr);

	return 0;
}

