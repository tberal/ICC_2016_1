#include <stdlib.h>

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

int main (int argc, char *argv[])
{
	int n = k = i = 0;
	double t = 0.0;
	char *output;
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
		if (argv[3] == '-i')
			 i = argv[4]
		else t = argv[4];
		strcpy(output, argv[6];
	}
	else
	{
		strcpy(output, argv[4];
	}

	double ** m = malloc ((k+1)*sizeof(double));
	for (int i=0; i<k; ++i)
	{
		m[i] = malloc(n*sizeof(double));
		generateRandomDiagonal (n, i, k, m[i]);
	}
	return 0;
}
