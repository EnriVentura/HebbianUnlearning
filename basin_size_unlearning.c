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
int N, N_samp;
long double alpha, strenghtN, D_maxstrenght;
int *sigma1;
int *sigma2;
FILE *fout1;
FILE *fout2;
FILE *fout3;

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

	int i;
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
	int i;
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
{
	//runs syncronous hopfield dynamics until convergence

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
{ //runs asyncronous hopfield dynamics
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
			printf("ABORT async_dynamics did not converge/n");
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
			fprintf(fout2, "\nAIUTO!!!\n\n");
			fprintf(fout3, "\nAIUTO!!!\n\n");

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
		printf("Please use 6 input parameters \n");
		printf("Usage: hebbian_unlearning_v2.exe N alpha strenghtN normalization_tipe n_samples delta_seed\n");
		exit(-9);
	}

	N = atoi(argv[1]);
	alpha = atof(argv[2]);
	strenghtN = atof(argv[3]);
	NORM_TYPE = argv[4];
	N_samp = atoi(argv[5]);
	int delta_seed = atoi(argv[6]);
	int seed = time(0) + delta_seed;
	srand48(seed);

	int i, j, l, t, on;
	int *sigma, *sigma_new, **csi;
	double **J, **J_av, **J_sigma, J_Av, J_Sigma;

	double time1, time2, time_tot = 0, time_tott = 0;

	printf("\n\nalpha = %Lg\n", alpha);
	printf("Strength = %Lg\n\n", strenghtN);

	int P = (int)(alpha * (double)N);
	double m, sigma_m;
	double ***over;

	double Overlap;
	int N_samp_over = 10;

	int N_samp_over_percept = 200;
	int index_min, index_max, mu_min, mu_max, over_min, over_max;

	double Deps_measure[3];
	char stringin[150];
	sprintf(stringin, "Dpoints_N%d_a%Lg_e%Lg.dat", N, alpha, strenghtN);
	FILE *take_D;
	take_D = fopen(stringin, "r");
	if (take_D == NULL)
	{
		printf("need files for D!");
		exit(1);
	}
	fscanf(take_D, "%lf\t%lf\t%lf\n", &Deps_measure[0], &Deps_measure[1], &Deps_measure[2]); // these 3 are D_in*eps/N, D_opt*eps/N, D_fin*eps/N
	fclose(take_D);

	D_maxstrenght = Deps_measure[2] + 0.05;
	long double strenght = strenghtN / N; //strenghtN è epsilon
	int D, D_max = (int)(D_maxstrenght / strenght);
	//N = 100, alpha = 0.3, eps = 0.01, {0.105, 0.235, 0.324}; alpha = 0.4, eps = 0.01, {0.231, 0.35, 0.426}; alpha = 0.5, eps = 0.01, {0.336, 0.445, 0.5}; alpha = 0.59, eps = 0.01, {0.529};
	// N = 150, alpha = 0.3, eps = 0.01, {0.119, 0.242, 0.330}; alpha = 0.4, eps = 0.01, {0.234, 0.353, 0.430}; alpha = 0.5, eps = 0.01, {0.371, 0.451, 0.507}; alpha = 0.59, eps = 0.01, {0.537};
	// N = 200, alpha = 0.3, eps = 0.01, {0.126, 0.249, 0.332}; alpha = 0.4, eps = 0.01, {0.241, 0.353, 0.429}; alpha = 0.5, eps = 0.01, {0.375, 0.456, 0.511}; alpha = 0.59, eps = 0.01. {0.534};
	//N = 300, alpha = 0.3, eps = 0.01, {0.136, 0.249 ,0.332}; N = 300, alpha = 0.4, eps = 0.01, {0.248, 0.395 ,0.428}; N = 300, alpha = 0.5, eps = 0.01, {0.38, 0.46 ,0.51}; N = 300, alpha = 0.6, eps = 0.01, {0.55}; N = 300, alpha = 0.59, eps = 0.01, {0.540};
	// N = 400, alpha = 0.3, eps = 0.01, {0.141, 0.25, 0.331}; alpha = 0.4, eps = 0.01, {0.249, 0.363, 0.428}; alpha = 0.5, eps = 0.01, {0.381, 0.463, 0.512}; alpha = 0.59, eps = 0.01, {0.542};

	// Initialization of the main configuration/interaction arrays

	sigma = (int *)malloc(N * sizeof(int));
	sigma_new = (int *)malloc(N * sizeof(int));

	sigma1 = (int *)malloc(N * sizeof(int)); //this is an ugly way of getting the malloc out of generate_initial(int **csi, int pattern, double p_in)
	if (sigma1 == NULL)
	{
		printf("malloc of sigma1 array failed.\n");
	}
	sigma2 = (int *)malloc(N * sizeof(int)); //this is an ugly way of getting the malloc out of generate_rand_initial()
	if (sigma2 == NULL)
	{
		printf("malloc of sigma2 array failed.\n");
	}

	csi = (int **)malloc(P * sizeof(int *));
	J = (double **)malloc(N * sizeof(double *));
	over = (double ***)malloc(N_samp * sizeof(double **));
	for (i = 0; i < N; i++)
	{
		J[i] = (double *)malloc(N * sizeof(double));
	}
	for (i = 0; i < P; i++)
	{
		csi[i] = (int *)malloc(N * sizeof(int));
	}
	for (i = 0; i < N_samp; i++)
	{
		over[i] = (double **)malloc(3 * sizeof(double *));
		for (int j = 0; j < 3; j++)
		{
			over[i][j] = (double *)malloc((int)(0.5 * N) * sizeof(double));
		}
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

	int initial_pattern = 0; //index of the test pattern with respect to the which we compute the overlap

	// Initialization of the output files
	char string[100], string2[100], string3[100];
	if (strcmp(NORM_TYPE, "NO_NORM") == 0)
	{
		sprintf(string, "unlearningV2NONORM_basin_Din_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string2, "unlearningV2NONORM_basin_Dopt_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string3, "unlearningV2NONORM_basin_Dfin_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
	}
	else if (strcmp(NORM_TYPE, "ROW_NORM") == 0)
	{
		sprintf(string, "unlearningV2NONORM_basin_Din_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string2, "unlearningV2NONORM_basin_Dopt_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string3, "unlearningV2NONORM_basin_Dfin_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
	}
	else if (strcmp(NORM_TYPE, "TOT_NORM") == 0)
	{
		sprintf(string, "unlearningV2NONORM_basin_Din_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string2, "unlearningV2NONORM_basin_Dopt_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
		sprintf(string3, "unlearningV2NONORM_basin_Dfin_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);
	}
	else
	{
		printf("please select a norm type: NO_NORM, ROW_NORM, TOT_NORM ");
		exit(1);
	}

	fout1 = fopen(string, "w");
	fout2 = fopen(string2, "w");
	fout3 = fopen(string3, "w");

	//fprintf(fout3, "Samples %d N %d alpha %Lg  D_maxstrenghtN %Lg D_maxstrenght %Lg norm %s \n", N_samp, N, alpha, strenghtN, D_maxstrenght, NORM_TYPE);

	for (i = 0; i < N_samp; i++)
	{ //Repetition of the run over N_samp realizations of disorder

		printf("Sample #%d of %d\n", i + 1, N_samp);
		generate_csi(P, csi);
		generate_J(P, csi, J);

		int t = 0;

		for (D = 0; D < D_max; D++)
		{ //Dreaming..
			//printf("D = %d\n", D);
			if (D > 0)
			{
				sigma = generate_rand_initial();
				async_dynamics(sigma, J);
				updateJ(J, sigma, strenght);
				normalizeJ(NORM_TYPE, J);
			}

			if (D > 0 && t < 3 && D % (int)(Deps_measure[t] / strenght) == 0)
			{
				printf("Deps/N = %Lg\n", (long double)D * strenght);
				for (int l = 0; l < (int)(0.5 * N); l++)
				{
					sigma_new = generate_initial(csi, initial_pattern, (1 - l * (double)(2 / (double)N)) * 0.5);
					async_dynamics(sigma_new, J);
					over[i][t][l] = overlap(csi, sigma_new, initial_pattern);
				}

				t++;
			}
		}
		printf("end sample \n");
	}

	// Analysis of the measures over the different samples

	for (t = 0; t < 3; t++)
	{
		for (int l = 0; l < (int)(0.5 * N); l++)
		{

			m = 0;
			sigma_m = 0;

			for (i = 0; i < N_samp; i++)
			{

				m = m + over[i][t][l] / (double)N_samp;
				sigma_m = sigma_m + over[i][t][l] * over[i][t][l] / (double)N_samp;
			}
			if (t == 0)
			{
				fprintf(fout1, "%Lg\t%lf\t%lf\n", (long double)(l * 2 / (double)N), m, sqrt((sigma_m - m * m) / (double)N_samp));
				fflush(fout1);
			}
			else if (t == 1)
			{
				fprintf(fout2, "%Lg\t%lf\t%lf\n", (long double)(l * 2 / (double)N), m, sqrt((sigma_m - m * m) / (double)N_samp));
				fflush(fout2);
			}
			else if (t == 2)
			{
				fprintf(fout3, "%Lg\t%lf\t%lf\n", (long double)(l * 2 / (double)N), m, sqrt((sigma_m - m * m) / (double)N_samp));
				fflush(fout3);
			}
		}
	}

	fclose(fout1);
	fclose(fout2);
	fclose(fout3);

	return 0;
}
