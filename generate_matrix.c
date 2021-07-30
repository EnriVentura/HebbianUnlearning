#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

int N, P;
long double alpha, p = 0.5;

void generate_csi(int P, int **csi)
{ //generate memories as a bernoulli process with probability p
	int i, j;
	for (i = 0; i < P; i++)
	{
		for (j = 0; j < N; j++)
		{
			if ((lrand48() / (double)RAND_MAX) < p)
			{
				csi[i][j] = 1;
			}
			else
			{
				csi[i][j] = -1;
			}
		}
	}
}

void generate_J(int P, int **csi, double **J)
{ //assembly connectivity matrix from memories to be symmetric with no autapses

	int i, j, k, s;

	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			J[i][j] = 0;
		}
	}

	for (i = 0; i < N; i++)
	{
		for (j = i + 1; j < N; j++)
		{
			s = 0;
			for (k = 0; k < P; k++)
			{
				s = s + (csi[k][i] * csi[k][j]);
			}
			J[i][j] = (double)s / (double)N;
			J[j][i] = J[i][j];
		}
		J[i][i] = 0;
	}
}

int main(int argc, char *argv[]){

	if (argc != 4)
	{
		printf("Please use 3 input parameters \n");
		printf("Usage: ./generate_matrix.exe N alpha delta_seed\n");
		exit(-9);
	}

	N = atoi(argv[1]);
	alpha = atof(argv[2]);
	int delta_seed = atoi(argv[3]);
	int seed = time(0) + delta_seed;
	srand48(seed);

	FILE *fout1;
	FILE *fout2;
	char string1[150], string2[150];
	sprintf(string1, "J_matrix_N%d_alpha%Lg.dat", N, alpha);
	sprintf(string2, "J_patts_N%d_alpha%Lg.dat", N, alpha);
	fout1 = fopen(string1, "w");
	fout2 = fopen(string2, "w");
	

	P = (int)(alpha*(double)N);

	int **csi;
	double **J;

	csi = (int **)malloc(P * sizeof(int *));
	J = (double **)malloc(N * sizeof(double *));
	for (int i = 0; i < N; i++)
	{
		J[i] = (double *)malloc(N * sizeof(double));
	}
	for (int i = 0; i < P; i++)
	{
		csi[i] = (int *)malloc(N * sizeof(int));
	}

	if (J == NULL)
	{
		printf("malloc of J matrix failed.\n");
		exit(EXIT_FAILURE);
	}
	else if (csi == NULL)
	{
		printf("malloc of csi matrix failed.\n");
		exit(EXIT_FAILURE);
	}
	generate_csi(P,csi);
	generate_J(P,csi,J);

	for(int i = 0; i < N; i++){
		for(int j = 0; j < N; j++){
			fprintf(fout1, "%lf\t", J[i][j]);
		}
		fprintf(fout1, "\n");
	}
	for(int i = 0; i < P;i++){
		for(int j = 0; j < N;j++){
			fprintf(fout2, "%d\t", csi[i][j]);
		}
		fprintf(fout2, "\n");
	}

	fclose(fout1);
	fclose(fout2);


}