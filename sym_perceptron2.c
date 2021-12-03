#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

int N, P, N_samp, max_iter, T_max;
long double lambda, alpha, den, c, p = 0.5, p_in = 0;

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

long double get_norm(double **J, int i)
{
	int j;

	long double norm = 0;
	for (j = 0; j < N; j++)
	{
		norm += J[i][j] * J[i][j];
	}
	norm = sqrt(norm);
	return norm;
}

int main(int argc, char *argv[])
{

	if (argc != 7)
	{
		printf("Please use 5 input parameters \n");
		printf("Usage: ./sym_perceptron2.exe N alpha lambda maxiter N_samp delta_seed\n");
		exit(-9);
	}

	N = atoi(argv[1]);
	alpha = atof(argv[2]);
	lambda = atof(argv[3]);
	T_max = atoi(argv[4]);
	N_samp = atoi(argv[5]);
	int delta_seed = atoi(argv[6]);
	int seed = time(0) + delta_seed;
	srand48(seed);

	double cc_max;
	int t, flagg;

	//double alpha_array[1] = {0.4};//{0.3, 0.4, 0.5};
	//double c_array[12] = {1.2959, 1.1265, 0.9875, 0.8695, 0.7681, 0.6795, 0.616, 0.601, 0.5306, 0.4671, 0.4095, 0.3565};
	double c_055[14] = {0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86}; //this was first attempt to go beyond gardner {0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79};
	int len_c_055 = 14;
	double c_05[14] = {0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94}; //this was first attempt to go beyond gardner{0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93};//this stopped at gardner's limit{0.7681, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0};
	int len_c_05 = 14;
	double c_04[14] = {1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 1.12, 1.13, 1.14}; //this was first attempt to go beyond gardner{0.95, 0.96, 0.97, 0.98, 0.99, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09};//this stopped at gardner's limit{0.9875, 0.9, 0.75, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0};
	int len_c_04 = 14;
	double c_03[20] = {1.26, 1.27, 1.28, 1.29, 1.30, 1.31, 1.32, 1.33, 1.34, 1.35, 1.36, 1.37, 1.38, 1.39, 1.40, 1.41, 1.42, 1.43, 1.44}; //this was first attempt to go beyond gardner{1.2, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.30, 1.31, 1.32, 1.33, 1.34, 1.35, 1.36, 1.37, 1.38, 1.39};this stopped at gardner's limit {1.285, 1.1, 0.9, 0.7, 0.5, 0.4, 0.3, 0.2, 0.1, 0};//{1.2959, 1.1, 0.9, 0.7, 0.5, 0.4, 0.3, 0.2, 0.1, 0};
	int len_c_03 = 20;

	char string[150];
	sprintf(string, "sym_perceptron_pSAT_N%d_alpha%Lg_lambda%Lg_maxiter%d_Nsamp%dseed%d.dat", N, alpha, lambda, T_max, N_samp, seed);
	FILE *fout1;
	fout1 = fopen(string, "w");
	printf("\nalpha = %Lg\n", alpha);

	P = (int)(alpha * (double)N);

	int **csi, **mask;
	double **J;
	csi = (int **)malloc(P * sizeof(int *));
	J = (double **)malloc(N * sizeof(double *));
	mask = (int **)malloc(P * sizeof(int *));
	for (int i = 0; i < N; i++)
	{
		J[i] = (double *)malloc(N * sizeof(double));
	}
	for (int i = 0; i < P; i++)
	{
		csi[i] = (int *)malloc(N * sizeof(int));
		mask[i] = (int *)malloc(N * sizeof(int));
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

	if (alpha == 0.3)
	{
		cc_max = len_c_03;
	}
	else if (alpha == 0.4)
	{
		cc_max = len_c_04;
	}
	else if (alpha == 0.5)
	{
		cc_max = len_c_05;
	}
	else if (alpha == 0.55)
	{
		cc_max = len_c_055;
	}
	else
	{
		printf("ABORT: alpha is not right");
		exit(9);
	}

	for (int cc = 0; cc < cc_max; cc++)
	{
		if (alpha == 0.3)
		{
			c = c_03[cc];
		}
		else if (alpha == 0.4)
		{
			c = c_04[cc];
			;
		}
		else if (alpha == 0.5)
		{
			c = c_05[cc];
		}
		else if (alpha == 0.55)
		{
			c = c_055[cc];
		}
		double freq = 0;
		for (int samp = 0; samp < N_samp; samp++)
		{
			generate_csi(P, csi);
			flagg = 0;
			t = 0;

			for (int i=0; i<N ; i++)
			{
				J[i][i]=0;
			}
			for (int i = 0; i < N; i++)
			{
				for (int j = i+1; j < N; j++)
				{
					J[i][j] = ((lrand48() / (double)RAND_MAX) - 0.5);
					J[j][i] = J[i][j];
				}
			}

			for (t = 0; t < T_max; t++)
			{
				int endcounter = 0; //this avoids iterating after all masks are 0!
				for (int j = 0; j < P; j++)
				{
					for (int i = 0; i < N; i++)
					{
						mask[j][i] = 0;
						double stab = 0;
						double norm = get_norm(J, i);
						for (int kk = 0; kk < N; kk++)
						{
							stab += J[i][kk] * csi[j][kk] * csi[j][i] / norm;
						}

						if (stab >= c)
						{
							mask[j][i] = 0;
							endcounter++;
						}
						else
						{
							mask[j][i] = 1;
						}
					}
				}
				if (endcounter == N * P)
				{
					freq = freq + 1 / (double)N_samp;
					break;
				}
				for (int i = 0; i < N; i++)
				{
					for (int j = i + 1; j < N; j++)
					{
						for (int k = 0; k < P; k++)
						{
							J[i][j] += lambda * (mask[k][i] + mask[k][j]) * csi[k][i] * csi[k][j];
						}
						J[j][i] = J[i][j];
					}
					J[i][i] = 0;
				}
			}
		}
		fprintf(fout1, "%Lf\t%lf\n", c, freq);
	}
}