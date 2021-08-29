#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

int N, P, N_samp, max_iter;
long double lambda, alpha, c, p = 0.5;
int *sigma1, *sigma2;

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

int *generate_initial(int **csi, int pattern, double p_in)
{ //generate an initial configuration with an average distance from pattern given by p_in
	int i, j;
	for (i = 0; i < N; i++)
	{
		if ((lrand48() / (double)RAND_MAX) < p_in)
		{
			sigma1[i] = (-1) * csi[pattern][i];
		}
		else
		{
			sigma1[i] = csi[pattern][i];
		}
	}
	return sigma1;
}

int *generate_rand_initial()
{ //random shooting generator

	int i, j;
	for (i = 0; i < N; i++)
	{
		if ((lrand48() / (double)RAND_MAX) < p)
		{
			sigma2[i] = 1;
		}
		else
		{
			sigma2[i] = -1;
		}
	}
	return sigma2;
}

void async_dynamics(int *sigma, double **J)
{ //runs asyncronous hopfield dynamics
	int i, j, k, l, flag = 0, time = 0;
	double field;

	while (time < 10000000)
	{
		flag = 0;
		i = (int)((lrand48() / (double)RAND_MAX) * (double)N); //even if i=N there is no problem
		for (int site = i; site < N; site++)
		{
			field = 0;
			for (int j = 0; j < N; j++)
			{
				field += J[site][j] * sigma[j];
			}
			if (field * sigma[site] < 0)
			{
				sigma[site] = -sigma[site];
				flag = 1;
				break;
			}
		}
		if (flag == 0)
		{
			for (int site = 0; site < i; site++)
			{
				field = 0;
				for (int j = 0; j < N; j++)
				{
					field += J[site][j] * sigma[j];
				}
				if (field * sigma[site] < 0)
				{
					sigma[site] = -sigma[site];
					flag = 1;
					break;
				}
			}
		}
		if (flag == 0)
		{
			break;
		}
		time++;
		//printf("time = %d\n", time);
		if (time == 99999)
		{
			printf("ABORT async_dynamics did not converge\n");
			exit(-9);
		}
	}
}

double overlap(int **csi, int *vec, int pattern)
{ //computes overlap between pattern and a given vector

	int i;
	double m = 0;

	for (i = 0; i < N; i++)
	{
		m = m + (double)csi[pattern][i] * (double)vec[i] / (double)N;
	}
	return m;
}

double overlap_patterns(int **csi, int *sigma, int i, int mu)
{

	double m = 0;

	for (int j = 0; j < N; j++)
	{
		m = m + (double)csi[mu][i] * csi[mu][j] * (double)sigma[i] * sigma[j] / (double)N;
	}
	return m;
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

	if (argc != 8)
	{
		printf("Please use 7 input parameters \n");
		printf("Usage: ./perceptron_algo.exe N alpha lambda max_stability max_number_of_iterations n_samples delta_seed\n");
		exit(-9);
	}

	N = atoi(argv[1]);
	alpha = atof(argv[2]);
	lambda = atof(argv[3]);
	c = atof(argv[4]);
	max_iter = atoi(argv[5]);
	N_samp = atoi(argv[6]);
	int delta_seed = atoi(argv[7]);
	int seed = time(0) + delta_seed;
	srand48(seed);

	char string[150];

	sprintf(string, "sym_perceptron_basin_N%d_alpha%Lg_lambda%Lg_maxstab%Lg_Nsamp%d_seed%d.dat", N, alpha, lambda, c, N_samp, seed);

	FILE *fout1;
	fout1 = fopen(string, "w");

	int i;

	P = (int)(alpha * (double)N);

	int **csi, *sigma, *sigma_new, **mask;
	double **J;

	sigma1 = (int *)malloc(N * sizeof(int)); //this is an ugly way of getting the malloc out of generate_initial(int **csi, int pattern, double p_in)
	if (sigma1 == NULL)
	{
		printf("malloc of sigma array failed.\n");
	}
	sigma2 = (int *)malloc(N * sizeof(int)); //this is an ugly way of getting the malloc out of generate_random_initial()
	if (sigma2 == NULL)
	{
		printf("malloc of sigma array failed.\n");
	}

	sigma = (int *)malloc(N * sizeof(int));
	sigma_new = (int *)malloc(N * sizeof(int));
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
	else if (sigma == NULL)
	{
		printf("malloc of sigma array failed.\n");
		exit(EXIT_FAILURE);
	}

	double **over;
	over = (double **)malloc(N_samp * sizeof(double *));
	for (int i = 0; i < N_samp; i++)
	{
		over[i] = (double *)malloc((int)(0.5 * N) * sizeof(double));
	}

	int initial_pattern = 0;

	for (int samp = 0; samp < N_samp; samp++)
	{
		printf("Sample #%d of %d\n", samp + 1, N_samp);
		generate_csi(P, csi);
		generate_J(P, csi, J);

		double k;
		int t, endcounter;
		for (t = 0; t < max_iter; t++)
		{
			endcounter = 0; //this avoids iterating after all masks are 0!
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
			if (t == max_iter-1){printf("algorithm did not converge in %d iterations", max_iter); exit (1);}
		}

		for (int l = 0; l < (int)(0.5 * N); l++) //this was l < (int)(0.5*N)+1
		{
			sigma_new = generate_initial(csi, initial_pattern, (1 - l * (double)(2 / (double)N)) * 0.5);
			async_dynamics(sigma_new, J);
			over[samp][l] = overlap(csi, sigma_new, initial_pattern);
		}
	}

	for (int l = 0; l < (int)(0.5 * N); l++)
	{
		double m = 0;
		double sigma_m = 0;
		for (i = 0; i < N_samp; i++)
		{
			m = m + over[i][l] / (double)N_samp;
			sigma_m = sigma_m + over[i][l] * over[i][l] / (double)N_samp;
		}
		fprintf(fout1, "%Lg\t%lf\t%lf\n", (long double)(l * 2 / (double)N), m, sqrt((sigma_m - m * m) / (double)N_samp));
		fflush(fout1);
	}
	fclose(fout1);
	return 0;
}