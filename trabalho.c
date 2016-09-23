#include <stdlib.h>
#include <string.h>
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
	return((4*pow(pi, 2))*(sin(2*pi*x)+(2*pi*(pi-x))));
}


double *generateResultVector(unsigned int N)
{
	double *b = malloc(N*sizeof(double));
	for (int i=0; i < N; ++i)
	{
		b[i] = f(i*pi/N);
	}
	return b;
}



int main (int argc, char *argv[])
{
	int n, k, i;
	double t = 0.0;
	double *b;
	char *output = malloc(256*sizeof(char));
	if (argc == 9)
	{ //todos os parametros
 		n = atoi(argv[1]);
		k = atoi(argv[2]);
		i = atoi(argv[4]);
		t = atof(argv[6]);
		strcpy (output, argv[8]);
	}
	else if (argc == 7)
	{
		if (strcmp(argv[3],"-i")==0)
			 i = atof(argv[4]);
		else t = atof(argv[4]);
		strcpy(output, argv[6]);
	}
	else
	{
		strcpy(output, argv[4]);
	}

	b = generateResultVector(n);
	double ** m = malloc ((k+1)*sizeof(double));
	for (int count=0; count<k; ++count)
	{
		m[count] = malloc(n*sizeof(double));
		generateRandomDiagonal (n, count, k, m[count]);
	}
	return 0;
}
