#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define thresh 0
#define p 0.5
#define p_in 0.0
#define EQ 20
#define T_MAX 1000000

char *NORM_TYPE;
int N, N_samp;
long double alpha, strenghtN, D_maxstrenght;

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

int *generate_initial(int **csi, int pattern)
{ //generate an initial configuration with an average distance from pattern given by p_in

	int i, j;
	int *sigma;
	sigma = (int *)malloc(N * sizeof(int));
	if (sigma == NULL)
	{
		printf("malloc of sigma array failed.\n");
	}

	for (i = 0; i < N; i++)
	{
		if ((lrand48() / (double)RAND_MAX) < p_in)
		{
			sigma[i] = (-1) * csi[pattern][i];
		}
		else
		{
			sigma[i] = csi[pattern][i];
		}
	}
	return sigma;
}

int *generate_rand_initial()
{ //random shooting generator

	int i, j;
	int *sigma;
	sigma = (int *)malloc(N * sizeof(int));
	if (sigma == NULL)
	{
		printf("malloc of sigma array failed.\n");
	}

	for (i = 0; i < N; i++)
	{
		if ((lrand48() / (double)RAND_MAX) < p)
		{
			sigma[i] = 1;
		}
		else
		{
			sigma[i] = -1;
		}
	}
	return sigma;
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
		if (time=99999){printf("ABORT async_dynamics did not converge/n"); exit (-9);}
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
	if (argc != 8)
	{
		printf("Please use 6 input parameters \n");
		printf("Usage: hebbian_unlearning_v2.exe N alpha strenghtN D_max*strenght normalization_tipe n_samples delta_seed\n");
		exit(-9);
	}

	N = atoi(argv[1]);
	alpha = atof(argv[2]);
	strenghtN = atof(argv[3]);
	D_maxstrenght = atof(argv[4]);
	NORM_TYPE = argv[5];
	N_samp = atoi(argv[6]);
	int delta_seed = atoi(argv[7]);
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
	double **over;

	double Overlap;
	int N_samp_over = 10;

	int N_samp_over_percept = 200;
	int index_min, index_max, mu_min, mu_max, over_min, over_max;

	long double strenght = strenghtN / N;
	int D, delta_D = (int)(0.005 / strenght), delta_D_histos = (int)(0.05 / strenght), D_max = (int)(D_maxstrenght / strenght);
	int N_steps = (int)((double)D_max / delta_D_histos) + 1;

	// Initialization of the main configuration/interaction arrays

	sigma = (int *)malloc(N * sizeof(int));
	sigma_new = (int *)malloc(N * sizeof(int));
	over = (double **)malloc(N_samp * sizeof(double));
	J_av = (double **)malloc(N_samp * sizeof(double));
	J_sigma = (double **)malloc(N_samp * sizeof(double));

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
	for (i = 0; i < N_samp; i++)
	{
		over[i] = (double *)malloc(((int)(D_max / delta_D) + 1) * sizeof(double));
		J_av[i] = (double *)malloc(((int)(D_max / delta_D) + 1) * sizeof(double));
		J_sigma[i] = (double *)malloc(((int)(D_max / delta_D) + 1) * sizeof(double));
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

	char string[100], string2[100], string3[100], string4[150], string5[150], string6[150], string7[150], string8[150];
	if (strcmp(NORM_TYPE, "NO_NORM") == 0)
	{
		sprintf(string, "unlearningV2NONORM_overlap_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d.dat", N, alpha, strenghtN, N_samp);
		sprintf(string2, "unlearningV2NONORM_Jmoments_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d.dat", N, alpha, strenghtN, N_samp);
		sprintf(string3, "unlearningV2NONORM_stabilities_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d.dat", N, alpha, strenghtN, N_samp);
		sprintf(string4, "unlearning_histoNONORM_SAT_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d.dat", N, alpha, strenghtN, N_samp);
		sprintf(string5, "unlearning_histoNONORM_UNSAT_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d.dat", N, alpha, strenghtN, N_samp);
		sprintf(string6, "unlearning_histoNONORM_ALL_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d.dat", N, alpha, strenghtN, N_samp);
		sprintf(string7, "unlearning_histoNONORM_MIN_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d.dat", N, alpha, strenghtN, N_samp);
		sprintf(string8, "unlearning_histoNONORM_MAX_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d.dat", N, alpha, strenghtN, N_samp);
	}
	else if (strcmp(NORM_TYPE, "ROW_NORM") == 0)
	{
		sprintf(string, "unlearningV2ROWNORM_overlap_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, N_samp);
		sprintf(string2, "unlearningV2ROWNORM_Jmoments_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, N_samp);
		sprintf(string3, "unlearningV2ROWNORM_stabilities_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, N_samp);
		sprintf(string4, "unlearning_histoROWNORM_SAT_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d.dat", N, alpha, strenghtN, N_samp);
		sprintf(string5, "unlearning_histoROWNORM_UNSAT_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d.dat", N, alpha, strenghtN, N_samp);
		sprintf(string6, "unlearning_histoROWNORM_ALL_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d.dat", N, alpha, strenghtN, N_samp);
		sprintf(string7, "unlearning_histoROWNORM_MIN_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d.dat", N, alpha, strenghtN, N_samp);
		sprintf(string8, "unlearning_histoROWNORM_MAX_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d.dat", N, alpha, strenghtN, N_samp);
	}
	else if (strcmp(NORM_TYPE, "TOT_NORM") == 0)
	{
		sprintf(string, "unlearningV2TOT_NORM_overlap_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, N_samp);
		sprintf(string2, "unlearningV2TOT_NORM_Jmoments_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, N_samp);
		sprintf(string3, "unlearningV2TOT_NORM_stabilities_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, N_samp);
		sprintf(string4, "unlearning_histoTOT_NORM_SAT_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d.dat", N, alpha, strenghtN, N_samp);
		sprintf(string5, "unlearning_histoTOT_NORM_UNSAT_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d.dat", N, alpha, strenghtN, N_samp);
		sprintf(string6, "unlearning_histoTOT_NORM_ALL_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d.dat", N, alpha, strenghtN, N_samp);
		sprintf(string7, "unlearning_histoTOT_NORM_MIN_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d.dat", N, alpha, strenghtN, N_samp);
		sprintf(string8, "unlearning_histoTOT_NORM_MAX_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d.dat", N, alpha, strenghtN, N_samp);
	}
	else
	{
		printf("please select a norm type: NO_NORM, ROW_NORM, TOT_NORM ");
		exit(1);
	}

	FILE *fout1;
	fout1 = fopen(string, "w");
	FILE *fout2;
	fout2 = fopen(string2, "w");
	FILE *fout3;
	fout3 = fopen(string3, "w");
	FILE *fout4;
	fout4 = fopen(string4, "w");
	FILE *fout5;
	fout5 = fopen(string5, "w");
	FILE *fout6;
	fout6 = fopen(string6, "w");
	FILE *fout7;
	fout7 = fopen(string7, "w");
	FILE *fout8;
	fout8 = fopen(string8, "w");

	//fprintf(fout3, "Samples %d N %d alpha %Lg  D_maxstrenghtN %Lg D_maxstrenght %Lg norm %s \n", N_samp, N, alpha, strenghtN, D_maxstrenght, NORM_TYPE);

	//marco
	long double ave_stability_sampled, ave_stability2_sampled, min_stability_sampled, min_stability2_sampled, max_stability_sampled, max_stability2_sampled, asymmetry_sampled, norm, n_sat_sampled, n_sat2_sampled;
	double **stability;
	double **ave_stability;
	double **max_stability;
	double **min_stability;
	double **asymm;
	double **n_sat;
	double **over_SAT, **over_UNSAT;

	n_sat = (double **)malloc(N_samp * sizeof(double *));
	for (i = 0; i < N_samp; i++)
	{
		n_sat[i] = (double *)malloc(((int)(D_max / delta_D) + 1) * sizeof(double));
	}

	asymm = (double **)malloc(N_samp * sizeof(double *));
	for (i = 0; i < N_samp; i++)
	{
		asymm[i] = (double *)malloc(((int)(D_max / delta_D) + 1) * sizeof(double));
	}
	stability = (double **)malloc(P * sizeof(double *));
	for (i = 0; i < P; i++)
	{
		stability[i] = (double *)malloc(N * sizeof(double));
	}
	ave_stability = (double **)malloc(N_samp * sizeof(double *));
	for (i = 0; i < N_samp; i++)
	{
		ave_stability[i] = (double *)malloc(((int)(D_max / delta_D) + 1) * sizeof(double));
	}
	max_stability = (double **)malloc(N_samp * sizeof(double *));
	for (i = 0; i < N_samp; i++)
	{
		max_stability[i] = (double *)malloc(((int)(D_max / delta_D) + 1) * sizeof(double));
	}
	min_stability = (double **)malloc(N_samp * sizeof(double *));
	for (i = 0; i < N_samp; i++)
	{
		min_stability[i] = (double *)malloc(((int)(D_max / delta_D) + 1) * sizeof(double));
	}

	if (n_sat == NULL)
	{
		printf("malloc of n_sat failed.\n");
		fprintf(fout3, "malloc of n_sat failed.\n");
		exit(EXIT_FAILURE);
	}

	if (asymm == NULL)
	{
		printf("malloc of asymmetry failed.\n");
		fprintf(fout3, "malloc of asymmetry failed.\n");
		exit(EXIT_FAILURE);
	}

	if (stability == NULL)
	{
		printf("malloc of stability failed.\n");
		fprintf(fout3, "malloc of stability failed.\n");
		exit(EXIT_FAILURE);
	}
	if (ave_stability == NULL)
	{
		printf("malloc of ave_stability failed.\n");
		fprintf(fout3, "malloc of ave_stability failed.\n");
		exit(EXIT_FAILURE);
	}
	if (max_stability == NULL)
	{
		printf("malloc of max_stability failed.\n");
		fprintf(fout3, "malloc of max_stability failed.\n");
		exit(EXIT_FAILURE);
	}
	if (min_stability == NULL)
	{
		printf("malloc of min_stability failed.\n");
		fprintf(fout3, "malloc of min_stability failed.\n");
		exit(EXIT_FAILURE);
	}

	int **overlap_SAT, **overlap_UNSAT, **overlap_ALL, **overlap_MIN, **overlap_MAX;
	int *N_SAT;

	overlap_SAT = (int **)malloc(N_steps * sizeof(int *));
	for (i = 0; i < N_steps; i++)
	{
		overlap_SAT[i] = (int *)malloc((2 * N + 1) * sizeof(int));
	}

	overlap_UNSAT = (int **)malloc(N_steps * sizeof(int *));
	for (i = 0; i < N_steps; i++)
	{
		overlap_UNSAT[i] = (int *)malloc((2 * N + 1) * sizeof(int));
	}
	overlap_ALL = (int **)malloc(N_steps * sizeof(int *));
	for (i = 0; i < N_steps; i++)
	{
		overlap_ALL[i] = (int *)malloc((2 * N + 1) * sizeof(int));
	}
	overlap_MIN = (int **)malloc(N_steps * sizeof(int *));
	for (i = 0; i < N_steps; i++)
	{
		overlap_MIN[i] = (int *)malloc((2 * N + 1) * sizeof(int));
	}
	overlap_MAX = (int **)malloc(N_steps * sizeof(int *));
	for (i = 0; i < N_steps; i++)
	{
		overlap_MAX[i] = (int *)malloc((2 * N + 1) * sizeof(int));
	}

	N_SAT = (int *)malloc(N_steps * sizeof(int *));

	if (overlap_SAT == NULL)
	{
		printf("malloc of SAT failed.\n");
		exit(EXIT_FAILURE);
	}

	if (overlap_UNSAT == NULL)
	{
		printf("malloc of UNSAT failed.\n");
		exit(EXIT_FAILURE);
	}

	if (overlap_ALL == NULL)
	{
		printf("malloc of ALL failed.\n");
		exit(EXIT_FAILURE);
	}
	if (overlap_MIN == NULL)
	{
		printf("malloc of MIN failed.\n");
		exit(EXIT_FAILURE);
	}
	if (overlap_MAX == NULL)
	{
		printf("malloc of MAX failed.\n");
		exit(EXIT_FAILURE);
	}
	if (N_SAT == NULL)
	{
		printf("malloc of n_SAT failed.\n");
		exit(EXIT_FAILURE);
	}

	for (i = 0; i < N_steps; i++)
	{
		N_SAT[i] = 0;
		for (j = 0; j < 2 * N + 1; j++)
		{
			overlap_SAT[i][j] = 0;
			overlap_UNSAT[i][j] = 0;
			overlap_ALL[i][j] = 0;
			overlap_MIN[i][j] = 0;
			overlap_MAX[i][j] = 0;
		}
	}

	for (i = 0; i < N_samp; i++)
	{ //Repetition of the run over N_samp realizations of disorder

		time1 = clock();
		printf("Sample #%d of %d\n", i + 1, N_samp);
		generate_csi(P, csi);
		generate_J(P, csi, J);

		int t = 0;
		int tt = 0;

		for (D = 0; D < D_max; D++)
		{ //Dreaming..

			if (D > 0)
			{
				sigma = generate_rand_initial();
				async_dynamics(sigma, J);
				updateJ(J, sigma, strenght);
				normalizeJ(NORM_TYPE, J);
			}

			if (D % delta_D == 0)
			{

				//printf("Deps/N = %Lg\n", D*strenght);
				J_av[i][t] = 0;
				J_sigma[i][t] = 0;
				for (l = 0; l < N; l++)
				{
					for (j = l + 1; j < N; j++)
					{
						J_av[i][t] = J_av[i][t] + J[l][j] / (double)(0.5 * N * (N - 1));
						J_sigma[i][t] = J_sigma[i][t] + J[l][j] * J[l][j] / (double)(0.5 * N * (N - 1));
					}
				}
				J_sigma[i][t] = J_sigma[i][t] - J_av[i][t] * J_av[i][t];

				//marco: computing stabilities and n_sat
				for (int mu = 0; mu < P; mu++)
				{
					for (int i = 0; i < N; i++)
					{
						stability[mu][i] = 0;
					}
				}

				for (int mu = 0; mu < P; mu++)
				{
					for (int i = 0; i < N; i++)
					{
						for (int j = 0; j < N; j++)
						{ //J[i][i]=0
							stability[mu][i] += J[i][j] * csi[mu][i] * csi[mu][j];
						}
						norm = get_norm(J, i);
						stability[mu][i] /= norm;
					}
				}

				ave_stability[i][t] = 0;
				max_stability[i][t] = stability[0][0];
				min_stability[i][t] = stability[0][0];
				n_sat[i][t] = 0;

				for (int mu = 0; mu < P; mu++)
				{
					for (int j = 0; j < N; j++)
					{
						ave_stability[i][t] += stability[mu][j];
						if (max_stability[i][t] < stability[mu][j])
						{
							max_stability[i][t] = stability[mu][j];
							index_max = j;
							mu_max = mu;
						}
						if (min_stability[i][t] > stability[mu][j])
						{
							min_stability[i][t] = stability[mu][j];
							index_min = j;
							mu_min = mu;
						}
						if (stability[mu][j] > 0)
						{
							n_sat[i][t]++;
						}
					}
				}

				ave_stability[i][t] /= (N * P);
				n_sat[i][t] /= (N * P);

				//marco computing asymmetry
				asymm[i][t] = 0;
				norm = 0;
				int k;
				for (k = 0; k < N; k++)
				{
					for (j = k + 1; j < N; j++)
					{
						norm += J[k][j] * J[k][j];
						asymm[i][t] += J[k][j] * J[j][k];
					}
				}
				asymm[i][t] /= norm;

				Overlap = 0;
				for (int l = 0; l < N_samp_over; l++)
				{
					sigma_new = generate_initial(csi, initial_pattern);
					async_dynamics(sigma_new, J);
					Overlap += overlap(csi, sigma_new, initial_pattern);
				}
				over[i][t] = Overlap / (double)N_samp_over;

				t++;
			}

			if (D % delta_D_histos == 0)
			{
				for (int l = 0; l < N_samp_over_percept; l++)
				{
					sigma = generate_rand_initial();
					async_dynamics(sigma, J);

					double stab_max = 0;
					double stab_min = 10000;

					for (int i = 0; i < N; i++)
					{
						for (int j = 0; j < P; j++)
						{

							double stab = 0;
							int over = 0;

							for (int k = 0; k < N; k++)
							{
								stab += J[i][k] * csi[j][k] * csi[j][i];
								over += csi[j][k] * csi[j][i] * sigma[k] * sigma[i];
							}
							if (stab <= 0)
							{
								overlap_UNSAT[tt][over + N]++;
							}
							else if (stab > 0)
							{
								overlap_SAT[tt][over + N]++;
								N_SAT[tt]++;
							}
							overlap_ALL[tt][over + N]++;
							if (stab >= stab_max)
							{
								stab_max = stab;
								over_max = over;
							}
							if (stab <= stab_min)
							{
								stab_min = stab;
								over_min = over;
							}
						}
					}
					overlap_MIN[tt][over_min + N]++;
					overlap_MAX[tt][over_max + N]++;
				}

				tt++;
			}
		}
		printf("end sample \n");
		time2 = clock();
		printf("time = %lf\n", (double)(time2 - time1) / ((double)CLOCKS_PER_SEC));
		time_tot += (time2 - time1) / ((double)CLOCKS_PER_SEC * N_samp);
	}
	printf("Total time needed = %lf\n", time_tot);

	// Analysis of the measures over the different samples

	for (t = 0; t < (int)(D_max / delta_D); t++)
	{
		ave_stability_sampled = 0;
		max_stability_sampled = 0;
		min_stability_sampled = 0;
		ave_stability2_sampled = 0;
		max_stability2_sampled = 0;
		min_stability2_sampled = 0;
		asymmetry_sampled = 0;
		n_sat_sampled = 0;
		n_sat2_sampled = 0;

		for (i = 0; i < N_samp; i++)
		{
			ave_stability_sampled += ave_stability[i][t];
			max_stability_sampled += max_stability[i][t];
			min_stability_sampled += min_stability[i][t];
			ave_stability2_sampled += ave_stability[i][t] * ave_stability[i][t];
			max_stability2_sampled += max_stability[i][t] * max_stability[i][t];
			min_stability2_sampled += min_stability[i][t] * min_stability[i][t];
			asymmetry_sampled += asymm[i][t];
			n_sat_sampled += n_sat[i][t];
			n_sat2_sampled += n_sat[i][t] * n_sat[i][t];
		}

		ave_stability_sampled /= N_samp;
		max_stability_sampled /= N_samp;
		min_stability_sampled /= N_samp;
		ave_stability2_sampled /= N_samp;
		max_stability2_sampled /= N_samp;
		min_stability2_sampled /= N_samp;
		asymmetry_sampled /= N_samp;
		n_sat_sampled /= N_samp;
		n_sat2_sampled /= N_samp;

		fprintf(fout3, "%d\t%d\t%Lg\t%Lg\t%Lg\t%d\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\n", N_samp, N, alpha, strenghtN, D_maxstrenght, t * delta_D, (long double)((t * delta_D) * (long double)(strenght)), ave_stability_sampled, (long double)(sqrt((ave_stability2_sampled - ave_stability_sampled * ave_stability_sampled) / (long double)N_samp)), max_stability_sampled, (long double)(sqrt((max_stability2_sampled - max_stability_sampled * max_stability_sampled) / (long double)N_samp)), min_stability_sampled, (long double)(sqrt((min_stability2_sampled - min_stability_sampled * min_stability_sampled) / (long double)N_samp)), asymmetry_sampled, n_sat_sampled, (long double)(sqrt((n_sat2_sampled - n_sat_sampled * n_sat_sampled) / (long double)N_samp)));
		fflush(fout3);
	}

	for (t = 0; t < (int)(D_max / delta_D); t++)
	{
		m = 0;
		J_Av = 0;
		sigma_m = 0;
		J_Sigma = 0;
		for (i = 0; i < N_samp; i++)
		{
			m = m + over[i][t] / (double)N_samp;
			J_Av = J_Av + J_av[i][t] / (double)N_samp;
			J_Sigma = J_Sigma + J_sigma[i][t] / (double)N_samp;
			sigma_m = sigma_m + over[i][t] * over[i][t] / (double)N_samp;
		}
		fprintf(fout1, "%Lg\t%lf\t%lf\n", (long double)(t * delta_D) * strenght, m, sqrt((sigma_m - m * m) / (double)N_samp));
		fprintf(fout2, "%Lg\t%lf\t%lf\n", (long double)(t * delta_D) * strenght, J_Av, J_Sigma);
		fflush(fout1);
		fflush(fout2);
	}

	for (int k = 0; k < N_steps; k++)
	{
		fprintf(fout4, "%d\t", N_SAT[k]);
		fflush(fout4);
		fprintf(fout5, "%d\t", P * N * N_samp * N_samp_over_percept - N_SAT[k]);
		fflush(fout5);
		fprintf(fout6, "%d\t", P * N * N_samp * N_samp_over_percept);
		fflush(fout6);
		fprintf(fout7, "%d\t", N_samp * N_samp_over_percept);
		fflush(fout7);
		fprintf(fout8, "%d\t", N_samp * N_samp_over_percept);
		fflush(fout8);
		for (i = 0; i < 2 * N + 1; i++)
		{
			fprintf(fout4, "%d\t", overlap_SAT[k][i]);
			fflush(fout4);
			fprintf(fout5, "%d\t", overlap_UNSAT[k][i]);
			fflush(fout5);
			fprintf(fout6, "%d\t", overlap_ALL[k][i]);
			fflush(fout6);
			fprintf(fout7, "%d\t", overlap_MIN[k][i]);
			fflush(fout7);
			fprintf(fout8, "%d\t", overlap_MAX[k][i]);
			fflush(fout8);
		}
		fprintf(fout4, "\n");
		fprintf(fout5, "\n");
		fprintf(fout6, "\n");
		fprintf(fout7, "\n");
		fprintf(fout8, "\n");
	}

	fclose(fout1);
	fclose(fout2);
	fclose(fout3);
	fclose(fout4);
	fclose(fout5);
	fclose(fout6);
	fclose(fout7);
	fclose(fout8);

	return 0;
}
