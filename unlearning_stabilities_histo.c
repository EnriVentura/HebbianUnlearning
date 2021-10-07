#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define thresh 0
#define p 0.5
#define EQ 20
#define T_MAX 1000000

char *NORM_TYPE;
int N, N_samp, n_samp_overlap = 10, delta_seed, n_bins = 1000;
;
long double alpha, strenghtN, D_maxstrenght, stab_inf = -1, bin_size = 0.005;
int *sigma1, *sigma2, *sigma_new;
double *field;
int **stabilities;
int delta_seed;
FILE *fout1;
FILE *fout2;
FILE *fout3;
FILE *fout4;
FILE *fout5;
FILE *fout6;
FILE *fout4bis;
FILE *fout5bis;
FILE *fout6bis;

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

double H(double **J, int *sigma)
{ //computes energy from Hopfield hamiltonian

	double H = 0;
	int i, j;

	for (i = 0; i < N; i++)
	{
		for (j = i + 1; j < N; j++)
		{
			H = H - J[i][j] * sigma[i] * sigma[j];
		}
	}

	return H;
}

int sync_dynamics(int *sigma, int *sigma_new, double **J, double *energ)
{ //runs syncronous hopfield dynamics until convergence

	int i, j, flag = 0, count;
	double field;
	int time = 0;

	while (flag < EQ && time < T_MAX)
	{

		//energ[time] = H(J, sigma);

		count = 0;
		for (i = 0; i < N; i++)
		{
			field = 0;
			for (j = 0; j < N; j++)
			{
				field = field + J[i][j] * (double)sigma[j];
			}
			if (field > 0)
			{
				sigma_new[i] = 1;
			}
			else
			{
				sigma_new[i] = -1;
			}
			count = count + sigma[i] * sigma_new[i];
		}
		if (count == N)
		{
			flag++;
		}
		for (i = 0; i < N; i++)
		{
			sigma[i] = sigma_new[i];
		}

		time++;
	}

	return time;
}

void async_dynamics(int *sigma, double **J)
{ //runs asyncronous hopfield dynamics until convergence
	int i, j, k, l, flag = 0, time = 0;
	double field;

	while (time < 100000)
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
			fprintf(fout1, "ABORT async_dynamics did not converge/n");
			exit(-9);
		}
	}
	for (int i = 0; i < N; i++)
	{
		double stab = 0;
		for (int k = 0; k < N; k++)
		{
			stab += J[i][k] * sigma[k] * sigma[i];
		}
		if (stab < 0)
		{
			printf("\nAIUTO!!!\n\n");
			fprintf(fout1, "\nAIUTO!!!\n\n");
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

void updateJ(double **J, int *sigma, double strenght)
{ //updates couplings according to Hebbian Unlearning
	int i, j;
	for (i = 0; i < N; i++)
	{
		for (j = i + 1; j < N; j++)
		{
			J[i][j] = J[i][j] - strenght * sigma[i] * sigma[j];
			J[j][i] = J[i][j];
		}
	}
}

//marco
void normalizeJ(char *type_of_norm, double **J)
{ //updates couplings according to Hebbian Unlearning
	if (strcmp(type_of_norm, "NO_NORM") == 0)
	{
	}
	if (strcmp(type_of_norm, "ROW_NORM") == 0)
	{
		int i, j;
		for (i = 0; i < N; i++)
		{
			long double norm = 0;
			for (j = 0; j < N; j++)
			{
				norm += J[i][j] * J[i][j];
			}
			norm = sqrt(norm);
			for (j = 0; j < N; j++)
			{
				J[i][j] /= norm;
			}
		}
	}
	if (strcmp(type_of_norm, "TOT_NORM") == 0)
	{
		int i, j;
		long double norm = 0;
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < N; j++)
			{
				norm += J[i][j] * J[i][j];
			}
		}
		norm = sqrt(norm);
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < N; j++)
			{
				J[i][j] /= norm;
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
		printf("Please use 6 input parameters\n");
		printf("Usage: ./get_J_unlearning.exe N alpha strenghtN normalization_tipe n_samples delta_seed \n");
		exit(-9);
	}

	N = atoi(argv[1]);
	alpha = atof(argv[2]);
	strenghtN = atof(argv[3]);
	NORM_TYPE = argv[4];
	N_samp = atoi(argv[5]);
	delta_seed = atoi(argv[6]);
	int seed = time(0) + delta_seed;
	srand48(seed);

	int i, j, l, t, on, i_max, i_min, pattern;
	int *sigma, **csi;
	double **J;

	double time1, time2, time_tot = 0, time_tott = 0;

	int P = (int)(alpha * (double)N);
	double m, m2;
	int ***freq;

	double Overlap;
	int N_samp_over = 10;

	int N_samp_over_percept = 200;
	int index_min, index_max, mu_min, mu_max, over_min, over_max;

	double Deps_measure[3];
	char stringin[150];
	sprintf(stringin, "Dpoints_N%d_a%Lg_e%Lg.dat", N, alpha, strenghtN);
	FILE *take_D;
	take_D = fopen(stringin, "r");

	fscanf(take_D, "%lf\t%lf\t%lf\n", &Deps_measure[0], &Deps_measure[1], &Deps_measure[2]);
	D_maxstrenght = Deps_measure[2] + 0.05;
	double strenght = strenghtN / (long double)N;
	int D, D_max = (int)(D_maxstrenght / strenght), D_in = (int)(Deps_measure[0] / strenght);
	// Initialization of the main configuration/interaction arrays

	char string1[100], string2[100], string3[100], string4[100], string4bis[100], string5[100], string5bis[100], string6[100], string6bis[100];
	if (strcmp(NORM_TYPE, "NO_NORM") == 0)
	{
		sprintf(string1, "unlearning_stabhisto_Din_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string2, "unlearning_stabhisto_Dopt_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string3, "unlearning_stabhisto_Dfin_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string4, "unlearning_overlaphisto_min_Din_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string4bis, "unlearning_overlaphisto_max_Din_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string5, "unlearning_overlaphisto_min_Dopt_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string5bis, "unlearning_overlaphisto_max_Dopt_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string6, "unlearning_overlaphisto_min_Dfin_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string6bis, "unlearning_overlaphisto_max_Dfin_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
	}
	else if (strcmp(NORM_TYPE, "ROW_NORM") == 0)
	{
		sprintf(string1, "unlearning_stabhisto_Din_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string2, "unlearning_stabhisto_Dopt_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string3, "unlearning_stabhisto_Dfin_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string4, "unlearning_overlaphisto_min_Din_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string4bis, "unlearning_overlaphisto_max_Din_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string5, "unlearning_overlaphisto_min_Dopt_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string5bis, "unlearning_overlaphisto_max_Dopt_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string6, "unlearning_overlaphisto_min_Dfin_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string6bis, "unlearning_overlaphisto_max_Dfin_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
	}
	else if (strcmp(NORM_TYPE, "TOT_NORM") == 0)
	{
		sprintf(string1, "unlearning_stabhisto_Din_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string2, "unlearning_stabhisto_Dopt_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string3, "unlearning_stabhisto_Dfin_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string4, "unlearning_overlaphisto_min_Din_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string4bis, "unlearning_overlaphisto_max_Din_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string5, "unlearning_overlaphisto_min_Dopt_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string5bis, "unlearning_overlaphisto_max_Dopt_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string6, "unlearning_overlaphisto_min_Dfin_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string6bis, "unlearning_overlaphisto_max_Dfin_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
	}
	else
	{
		printf("please select a norm type: NO_NORM, ROW_NORM, TOT_NORM ");
		exit(1);
	}

	fout1 = fopen(string1, "w");
	fout2 = fopen(string2, "w");
	fout3 = fopen(string3, "w");
	fout4 = fopen(string4, "w");
	fout5 = fopen(string5, "w");
	fout6 = fopen(string6, "w");
	fout4bis = fopen(string4bis, "w");
	fout5bis = fopen(string5bis, "w");
	fout6bis = fopen(string6bis, "w");

	sigma = (int *)malloc(N * sizeof(int));

	// Initialization of the main configuration/interaction arrays
	sigma1 = (int *)malloc(N * sizeof(int));
	if (sigma1 == NULL)
	{
		printf("malloc of sigma1 array failed.\n");
	}
	sigma2 = (int *)malloc(N * sizeof(int));
	if (sigma2 == NULL)
	{
		printf("malloc of sigma2 array failed.\n");
	}
	sigma_new = (int *)malloc(N * sizeof(int));
	if (sigma_new == NULL)
	{
		printf("malloc of sigma array failed.\n");
	}

	csi = (int **)malloc(P * sizeof(int *));
	J = (double **)malloc(N * sizeof(double *));
	for (i = 0; i < N; i++)
	{
		J[i] = (double *)malloc(N * sizeof(double));
	}
	for (i = 0; i < P; i++)
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
	else if (sigma == NULL)
	{
		printf("malloc of sigma array failed.\n");
		exit(EXIT_FAILURE);
	}

	int ***freq_min, ***freq_max;
	freq_min = (int ***)malloc(3 * sizeof(int **));
	freq_max = (int ***)malloc(3 * sizeof(int **));

	for (int i = 0; i < 3; i++)
	{
		freq_min[i] = (int **)malloc(((int)(0.5 * N) + 1) * sizeof(int *));
		freq_max[i] = (int **)malloc(((int)(0.5 * N) + 1) * sizeof(int *));
	}
	for (int j = 0; j < 3; j++)
	{
		for (int i = 0; i < (int)(0.5 * N) + 1; i++)
		{
			freq_min[j][i] = (int *)malloc((2 * N + 1) * sizeof(int));
			freq_max[j][i] = (int *)malloc((2 * N + 1) * sizeof(int));
		}
	}

	stabilities = (int **)malloc(sizeof(int *) * 3);
	for (int t = 0; t < 3; t++)
	{
		stabilities[t] = (int *)malloc(sizeof(int) * n_bins);
	}

	for (int samp = 0; samp < N_samp; samp++)
	{ //Repetition of the run over N_samp realizations of disorder

		generate_csi(P, csi);
		generate_J(P, csi, J);
		printf("alpha = %Lg\tsample #%d\n", alpha, samp + 1);

		int t = 0;

		for (D = 0; D < D_max; D++)
		{ //Dreaming..
			//printf("Deps/N = %Lg\n", (long double)(D*strenght));

			if (D > 0)
			{
				sigma = generate_rand_initial();
				async_dynamics(sigma, J);
				updateJ(J, sigma, strenght);
				normalizeJ(NORM_TYPE, J);
			}

			if (D > 0 && t < 3 && D % (int)(Deps_measure[t] / strenght) == 0)
			{

				double stab_min = 100000;
				double stab_max = -100000;
				for (int j = 0; j < P; j++)
				{
					for (int i = 0; i < N; i++)
					{
						double stab = 0;
						double norm = get_norm(J, i);
						for (int kk = 0; kk < N; kk++)
						{
							stab += J[i][kk] * csi[j][kk] * csi[j][i] / norm;
						}
						if (stab >= stab_max)
						{
							i_max = i;
							mu_max = j;
						}
						if (stab <= stab_min)
						{
							i_min = i;
							mu_min = j;
						}
						if (stab < stab_inf || stab > stab_inf + n_bins * bin_size)
						{
							fprintf(fout1, "ABORT stability outside binning interval \n");
							exit(9);
						}
						else
						{
							stabilities[t][(int)((stab - stab_inf) / bin_size)]++;
						}
					}
				}

				for (int l = 0; l < (int)(0.5 * N) + 1; l++) //this was l < (int)(0.5*N)+1
				{
					for (int k = 0; k < n_samp_overlap; k++)
					{
						sigma_new = generate_initial(csi, mu_min, (1 - l * (double)(2 / (double)N)) * 0.5);
						async_dynamics(sigma_new, J);
						freq_min[t][l][(int)((overlap(csi, sigma_new, mu_min) + 1) * (double)N)]++;

						sigma_new = generate_initial(csi, mu_max, (1 - l * (double)(2 / (double)N)) * 0.5);
						async_dynamics(sigma_new, J);
						freq_max[t][l][(int)((overlap(csi, sigma_new, mu_max) + 1) * (double)N)]++;
					}
				}
				t++;
			}
		}
	}//samp

	//printing
	for (int t = 0; t < 3; t++)
	{
		if (t == 0)
		{
			for (int i = 0; i < n_bins; i++)
			{
				fprintf(fout1, "binned_stability %Lg counter %d \n", stab_inf + i * bin_size, stabilities[t][i]);
			}
			for (int j = 0; j < (int)(0.5 * N) + 1; j++)
			{

				fprintf(fout4, "%Lg\t", (long double)(j * 2 / (double)N));
				fflush(fout4);
				fprintf(fout4bis, "%Lg\t", (long double)(j * 2 / (double)N));
				fflush(fout4bis);
				for (int i = 0; i < 2 * N + 1; i++)
				{
					fprintf(fout4, "%d\t", freq_min[t][j][i]);
					fflush(fout4);
					fprintf(fout4bis, "%d\t", freq_max[t][j][i]);
					fflush(fout4bis);
				}
				fprintf(fout4bis, "\n");

				fprintf(fout4, "\n");
				fflush(fout4);
			}
		}
		else if (t == 1)
		{
			for (int i = 0; i < n_bins; i++)
			{
				fprintf(fout2, "binned_stability %Lg counter %d \n", stab_inf + i * bin_size, stabilities[t][i]);
			}
			for (int j = 0; j < (int)(0.5 * N) + 1; j++)
			{
				fprintf(fout5, "%Lg\t", (long double)(j * 2 / (double)N));
				fflush(fout5);
				fprintf(fout5bis, "%Lg\t", (long double)(j * 2 / (double)N));
				fflush(fout5bis);
				for (int i = 0; i < 2 * N + 1; i++)
				{
					fprintf(fout5, "%d\t", freq_min[t][j][i]);
					fflush(fout5);
					fprintf(fout5bis, "%d\t", freq_max[t][j][i]);
					fflush(fout5bis);
				}
				fprintf(fout5bis, "\n");

				fprintf(fout5, "\n");
				fflush(fout5);
			}
		}
		else if (t == 2)
		{
			for (int i = 0; i < n_bins; i++)
			{
				fprintf(fout3, "binned_stability %Lg counter %d \n", stab_inf + i * bin_size, stabilities[t][i]);
			}
			for (int j = 0; j < (int)(0.5 * N) + 1; j++)
			{
				fprintf(fout6, "%Lg\t", (long double)(j * 2 / (double)N));
				fflush(fout6);
				fprintf(fout6bis, "%Lg\t", (long double)(j * 2 / (double)N));
				fflush(fout6bis);
				for (int i = 0; i < 2 * N + 1; i++)
				{
					fprintf(fout6, "%d\t", freq_min[t][j][i]);
					fflush(fout6);
					fprintf(fout6bis, "%d\t", freq_max[t][j][i]);
					fflush(fout6bis);
				}
				fprintf(fout6bis, "\n");

				fprintf(fout6, "\n");
				fflush(fout6);
			}
		}
	}

	fclose(fout1);
	fclose(fout2);
	fclose(fout3);
	fclose(fout4);
	fclose(fout5);
	fclose(fout6);
	fclose(fout4bis);
	fclose(fout5bis);
	fclose(fout6bis);
	return 0;
}