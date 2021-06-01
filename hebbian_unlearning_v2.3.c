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

int myoverlap;

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

int async_dynamics(int *sigma, double **J)
{ //runs asyncronous hopfield dynamics until convergence

	int i, j, k, flag = 0, time = 0, count;
	double field;
	int *sigma_new;
	sigma_new = (int *)malloc(N * sizeof(int));
	if (sigma_new == NULL)
	{
		printf("malloc of sigma array failed.\n");
	}

	while (flag == 0)
	{
		//printf("time = %d\n", time);

		do
		{
			i = (int)((lrand48() / (double)RAND_MAX) * (double)N);
		} while (i == N);
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
			//printf("ok\n");
		}

		count = 0;
		for (k = 0; k < N; k++)
		{
			count = count + sigma[k] * sigma_new[k];
		}
		if (count == N)
		{
			flag++;
		}

		sigma[i] = sigma_new[i];

		time++;
		//printf("time = %d\n", time);
	}
	return time;
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

int overlap_patterns(int **csi, int *sigma, int i, int mu)
{ //marco: this now returns unnormaized overlaps

	int m = 0;

	for (int j = 0; j < N; j++)
	{
		m += csi[mu][i] * csi[mu][j] * sigma[i] * sigma[j];
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

int main(int argc, char *argv[])
{
	if (argc != 7)
	{
		printf("Please use 6 input parameters\n");
		printf("Usage: hebbian_unlearning_v2.exe N alpha strenghtN D_max*strenght normalization_tipe n_samples \n");
		exit(-9);
	}

	N = atoi(argv[1]);
	alpha = atof(argv[2]);
	strenghtN = atof(argv[3]);
	D_maxstrenght = atof(argv[4]);
	NORM_TYPE = argv[5];
	N_samp = atoi(argv[6]);

	int seed = time(0);
	srand48(seed);

	int i, j, l, t, on;
	int *sigma, *sigma_new, **csi;
	double **J, **J_av, **J_sigma, J_Av, J_Sigma;

	int P = (int)(alpha * (double)N);
	double m, sigma_m;
	double **over;

	double Overlap;
	int N_samp_over = 10;

	int N_samp_over_percept = 10;
	int index_min, index_max, mu_min, mu_max;
	int flag_up, flag_down, histo_times[4];
	histo_times[0] = 0;
	histo_times[1] = 0;
	histo_times[2] = 0;

	long double strenght = strenghtN / N; //marco
	int D, delta_D = (int)(0.01 / strenght);
	int D_max = (int)(D_maxstrenght / strenght);

	int time_supp; //supporting time variable for async routine

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

	char string[100], string2[100], string3[100], string4[150], string5[150];
	if (strcmp(NORM_TYPE, "NO_NORM") == 0)
	{
		sprintf(string, "unlearningV2NONORM_overlap_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, seed);
		sprintf(string2, "unlearningV2NONORM_Jmoments_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, seed);
		sprintf(string3, "unlearningV2NONORM_stabilities_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, seed);
		sprintf(string4, "unlearningV2NONORM_perceptron_overlap_histo_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, seed);
	}
	else if (strcmp(NORM_TYPE, "ROW_NORM") == 0)
	{
		sprintf(string, "unlearningV2ROWNORM_overlap_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, seed);
		sprintf(string2, "unlearningV2ROWNORM_Jmoments_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, seed);
		sprintf(string3, "unlearningV2ROWNORM_stabilities_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, seed);
		sprintf(string4, "unlearningV2ROWNORM_perceptron_overlap_histo_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, seed);
	}
	else if (strcmp(NORM_TYPE, "TOT_NORM") == 0)
	{
		sprintf(string, "unlearningV2TOT_NORM_overlap_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, seed);
		sprintf(string2, "unlearningV2TOT_NORM_Jmoments_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, seed);
		sprintf(string3, "unlearningV2TOT_NORM_stabilities_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, seed);
		sprintf(string4, "unlearningV2TOT_NORM_perceptron_overlap_histo_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, seed);
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
	fout3 = fopen(string3, "w"); //marco
	FILE *fout4;
	fout4 = fopen(string4, "w");

	fprintf(fout3, "Samples %d N %d alpha %Lg  D_maxstrenghtN %Lg D_maxstrenght %Lg norm %s \n", N_samp, N, alpha, strenghtN, D_maxstrenght, NORM_TYPE);

	//marco
	long double ave_stability_sampled, ave_stability2_sampled, min_stability_sampled, min_stability2_sampled, max_stability_sampled, max_stability2_sampled, asymmetry_sampled, norm, n_sat_sampled, n_sat2_sampled, steps_sampled, steps2_sampled, phys_time_sampled, phys_time2_sampled;
	double **stability;
	double **ave_stability;
	double **max_stability;
	double **min_stability;
	double **asymm;
	double **n_sat;
	int **over_SAT, **over_UNSAT;
	int **Steps;
	double **Phys_Time;

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
	over_SAT = (int **)malloc(3 * sizeof(int *));
	for (i = 0; i < 3; i++)
	{
		over_SAT[i] = (int *)malloc((N + 1) * sizeof(int));
	}
	over_UNSAT = (int **)malloc(3 * sizeof(int *));
	for (i = 0; i < 3; i++)
	{
		over_UNSAT[i] = (int *)malloc((N + 1) * sizeof(int));
	}
	Steps = (int **)malloc(N_samp * sizeof(int *));
	for (i = 0; i < N_samp; i++)
	{
		Steps[i] = (int *)malloc(((int)(D_max / delta_D) + 1) * sizeof(int));
	}
	Phys_Time = (double **)malloc(N_samp * sizeof(double *));
	for (i = 0; i < N_samp; i++)
	{
		Phys_Time[i] = (double *)malloc(((int)(D_max / delta_D) + 1) * sizeof(double));
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
	if (over_SAT == NULL)
	{
		printf("malloc of over_SAT failed.\n");
		fprintf(fout4, "malloc of over_SAT failed.\n");
		exit(EXIT_FAILURE);
	}
	if (over_UNSAT == NULL)
	{
		printf("malloc of over_UNSAT failed.\n");
		fprintf(fout4, "malloc of over_UNSAT failed.\n");
		exit(EXIT_FAILURE);
	}
	if (Steps == NULL)
	{
		printf("malloc of Steps failed.\n");
		fprintf(fout4, "malloc of Steps failed.\n");
		
		exit(EXIT_FAILURE);
	}
	if (Phys_Time == NULL)
	{
		printf("malloc of Phys_Time failed.\n");
		fprintf(fout4, "malloc of Phys_Time failed.\n");
	}
		

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < (N + 1); j++)
		{
			over_SAT[i][j] = 0;
			over_UNSAT[i][j] = 0;
		}
	}

	for (i = 0; i < N_samp; i++)
	{ //Repetition of the run over N_samp realizations of disorder

		printf("Sample #%d of %d\n", i + 1, N_samp);
		generate_csi(P, csi);
		generate_J(P, csi, J);

		int t = 0;
		int r = 0;
		flag_up = 0;
		flag_down = 0;

		for (D = 0; D < D_max; D++)
		{ //Dreaming..

			if (D > 0)
			{
				
				if(D % delta_D == 0){
					sigma = generate_rand_initial();
					clock_t start = clock();
					Steps[i][t] = async_dynamics(sigma, J);
					clock_t end = clock();
					updateJ(J, sigma, strenght);
					Phys_Time[i][t] = (double)(end - start) / CLOCKS_PER_SEC;
				}else{
					sigma = generate_rand_initial();
					time_supp = async_dynamics(sigma, J);
					updateJ(J, sigma, strenght);
				}
				
				normalizeJ(NORM_TYPE, J);
			}

			if (D % delta_D == 0)
			{ //Check and measure
				// printf("D_max %d D_maxstrenght %Lg strenght %Lg deltaD %d step%d \n", D_max, D_maxstrenght, strenght, delta_D, t);
				// printf("%d D_max %d D_maxstrenght %Lg strenght %Lg deltaD %d totsteps %d step%d \n", (int)(1/0.001), D_max, D_maxstrenght, strenght, delta_D, (int)(D_max / delta_D), t);

				//printf("Deps/N = %Lg\n", D*strenght); //test
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

				//filling histos: CARE average times are non normalized here
				if (t == 0)
				{
					for (int l = 0; l < N_samp_over_percept; l++)
					{
						sigma_new = generate_rand_initial();
						time_supp = async_dynamics(sigma_new, J);
						over_SAT[0][(overlap_patterns(csi, sigma_new, index_max, mu_max) + N) / 2]++;
						over_UNSAT[0][(overlap_patterns(csi, sigma_new, index_min, mu_min) + N) / 2]++;
					}
				}

				if (t != 0 && min_stability[i][t] >= 0)
				{
					if (flag_up < 1)
					{
						histo_times[1] += t;
						//printf("time_1 = %Lg\n", D * strenght); //test
						for (int l = 0; l < N_samp_over_percept; l++)
						{
							sigma_new = generate_rand_initial();
							time_supp = async_dynamics(sigma_new, J);
							over_SAT[1][(overlap_patterns(csi, sigma_new, index_max, mu_max) + N) / 2]++;
							over_UNSAT[1][(overlap_patterns(csi, sigma_new, index_min, mu_min) + N) / 2]++;
						}
					}
					flag_up++;
				}

				if (flag_up > 5 && min_stability[i][t] < 0)
				{
					if (flag_down < 1)
					{
						histo_times[2] += t;
						//printf("time_2 = %Lg %d \n", D * strenght, flag_down); //test
						for (int l = 0; l < N_samp_over_percept; l++)
						{
							sigma_new = generate_rand_initial();
							time_supp = async_dynamics(sigma_new, J);
							over_SAT[2][(overlap_patterns(csi, sigma_new, index_max, mu_max) + N) / 2]++;
							over_UNSAT[2][(overlap_patterns(csi, sigma_new, index_min, mu_min) + N) / 2]++;
						}
					}
					flag_down++;
				}

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
					time_supp = async_dynamics(sigma_new, J);
					Overlap += overlap(csi, sigma_new, initial_pattern);
				}
				over[i][t] = Overlap / (double)N_samp_over;
				
				t++;
			}
		}
		printf("end sample \n");
	}

	// Analysis of the measures over the different samples

	//marco: analize stabilities and n_sat and asymmetry
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
		steps_sampled = 0;
		steps2_sampled = 0;
		phys_time_sampled = 0;
		phys_time2_sampled = 0;

		for (i = 0; i < N_samp; i++)
		{
			ave_stability_sampled += ave_stability[i][t];
			max_stability_sampled += max_stability[i][t];
			min_stability_sampled += min_stability[i][t];
			ave_stability2_sampled += ave_stability[i][t]*ave_stability[i][t];
			max_stability2_sampled += max_stability[i][t]*max_stability[i][t];
			min_stability2_sampled += min_stability[i][t]*min_stability[i][t];
            asymmetry_sampled += asymm[i][t];
			n_sat_sampled += n_sat[i][t];
			n_sat2_sampled += n_sat[i][t]*n_sat[i][t];
			steps_sampled += Steps[i][t];
			steps2_sampled += Steps[i][t]*Steps[i][t];
			phys_time_sampled += Phys_Time[i][t];
			phys_time2_sampled += Phys_Time[i][t]*Phys_Time[i][t];

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
		steps_sampled /= N_samp;
		steps2_sampled /= N_samp;
		phys_time_sampled /= N_samp;
		phys_time2_sampled /= N_samp;

		fprintf(fout3, "Samples %d N %d alpha %Lg strenghtN %Lg D_maxstrenght %Lg dream %d ave_stability %Lg sigma_ave_stability %lf max_stability_sampled %Lg sigma_max_stability %lf min_stability_sampled %Lg sigma_min_stability %lf asymmetry_sampled %Lg n_sat_sampled %Lg sigma_n_sat %lf\tdynamics_steps %Lg\tsigma_dynamics_steps %lf\tphysical time %Lg\tsigma_physical_time %lf\n", N_samp, N, alpha, strenghtN, D_maxstrenght, t * delta_D, ave_stability_sampled, sqrt((ave_stability2_sampled - ave_stability_sampled*ave_stability_sampled )/(double)N_samp), max_stability_sampled, sqrt((max_stability2_sampled - max_stability_sampled*max_stability_sampled )/(double)N_samp),min_stability_sampled, sqrt((min_stability2_sampled - min_stability_sampled*min_stability_sampled )/(double)N_samp),asymmetry_sampled, n_sat_sampled, sqrt((n_sat2_sampled-n_sat_sampled*n_sat_sampled)/(double)N_samp), steps_sampled, sqrt((steps2_sampled-steps_sampled*steps_sampled)/(double)N_samp), phys_time_sampled, sqrt((phys_time2_sampled-phys_time_sampled*phys_time_sampled)/(double)N_samp));
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
		fprintf(fout1, "Dream %d m %lf sigma_m %lf\n", t * delta_D, m, sqrt((sigma_m - m * m) / (double)N_samp));
		fprintf(fout2, "Dream %d mean_J %lf sigma_J %lf\n", t * delta_D, J_Av, J_Sigma);
		fflush(fout1);
		fflush(fout2);
	}

	//marco:printing histograms
	for(int myoverlap=0; myoverlap<(N+1); myoverlap++){
	if(over_SAT[0][myoverlap]!=0){fprintf(fout4, "N %d At_time 0 over %d samples perceptron_overlap_SAT %d counter %d\n",N, N_samp*N_samp_over_percept,myoverlap*2-N,over_SAT[0][myoverlap]);}
}
for(int myoverlap=0; myoverlap<(N+1); myoverlap++){
	if(over_UNSAT[0][myoverlap]!=0){fprintf(fout4, "N %d delta_D %d strenght %Lg At_time 0 over %d samples perceptron_overlap_UNSAT %d counter %d\n",N, delta_D, strenght, N_samp*N_samp_over_percept,myoverlap*2-N,over_UNSAT[0][myoverlap]);}
}

for(int myoverlap=0; myoverlap<(N+1); myoverlap++){
	if(over_SAT[1][myoverlap]!=0){fprintf(fout4, "N %d delta_D %d strenght %Lg ave_t_min %Lg over %d samples perceptron_overlap_SAT %d counter %d\n",N, delta_D, strenght,(long double) histo_times[1]*delta_D*strenght/N_samp, N_samp*N_samp_over_percept,myoverlap*2-N,over_SAT[1][myoverlap]);}
}
for(int myoverlap=0; myoverlap<(N+1); myoverlap++){
	if(over_UNSAT[1][myoverlap]!=0){fprintf(fout4, "N %d delta_D %d strenght %Lg ave_t_min %Lg over %d samples perceptron_overlap_UNSAT %d counter %d\n",N, delta_D, strenght,(long double) histo_times[1]*delta_D*strenght/N_samp, N_samp*N_samp_over_percept,myoverlap*2-N,over_UNSAT[1][myoverlap]);}
}

for(int myoverlap=0; myoverlap<(N+1); myoverlap++){
	if(over_SAT[2][myoverlap]!=0){fprintf(fout4, "N %d delta_D %d strenght %Lg ave_t_max %Lg over %d samples perceptron_overlap_SAT %d counter %d\n",N, delta_D, strenght,(long double) histo_times[2]*delta_D*strenght/N_samp, N_samp*N_samp_over_percept,myoverlap*2-N,over_SAT[2][myoverlap]);}
}
for(int myoverlap=0; myoverlap<(N+1); myoverlap++){
	if(over_UNSAT[2][myoverlap]!=0){fprintf(fout4, "N %d delta_D %d strenght %Lg ave_t_max %Lg over %d samples perceptron_overlap_UNSAT %d counter %d\n",N,delta_D, strenght,(long double) histo_times[2]*delta_D*strenght/N_samp, N_samp*N_samp_over_percept,myoverlap*2-N,over_UNSAT[2][myoverlap]);}
}
	fclose(fout1);
	fclose(fout2);
	fclose(fout3);
	fclose(fout4);

	return 0;
}
