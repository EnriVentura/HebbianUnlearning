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

char * NORM_TYPE;
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

		//time++;
		//printf("time = %d\n", time);
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
void normalizeJ(char* type_of_norm, double **J)
{ //updates couplings according to Hebbian Unlearning
	if (strcmp(type_of_norm , "NO_NORM") == 0)
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

	long double strenght = strenghtN / N;											   //marco
	int D, delta_D = (int)(0.01 / strenght);
    int D_max = (int) (D_maxstrenght / strenght); 
 


    //printf("D_max %d \n", D_max); POOP
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
		over[i] = (double *)malloc(((int)(D_max / delta_D)+1) * sizeof(double));
		J_av[i] = (double *)malloc(((int)(D_max / delta_D)+1) * sizeof(double));
		J_sigma[i] = (double *)malloc(((int)(D_max / delta_D)+1) * sizeof(double));
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
if(strcmp(NORM_TYPE, "NO_NORM")  == 0 ){
	sprintf(string, "unlearningV2NONORM_overlap_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN , seed);
	sprintf(string2, "unlearningV2NONORM_Jmoments_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN , seed);
	sprintf(string3, "unlearningV2NONORM_stabilities_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN , seed);
}
else if(strcmp(NORM_TYPE, "ROW_NORM") == 0 ) {
	sprintf(string, "unlearningV2ROWNORM_overlap_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, seed);
	sprintf(string2, "unlearningV2ROWNORM_Jmoments_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, seed);
	sprintf(string3, "unlearningV2ROWNORM_stabilities_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, seed);
}
else if(strcmp(NORM_TYPE, "TOT_NORM") == 0 ) {
	sprintf(string, "unlearningV2TOT_NORM_overlap_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, seed);
	sprintf(string2, "unlearningV2TOT_NORM_Jmoments_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, seed);
	sprintf(string3, "unlearningV2TOT_NORM_stabilities_N%d_alpha%Lg_strenghtN%Lg_seed%d.dat", N, alpha, strenghtN, seed);
}
else{printf("please select a norm type: NO_NORM, ROW_NORM, TOT_NORM "); exit (1);}

	FILE *fout1;
	fout1 = fopen(string, "w");
	FILE *fout2;
	fout2 = fopen(string2, "w");
	FILE *fout3;
	fout3 = fopen(string3, "w"); //marco

fprintf(fout3, "Samples %d N %d alpha %Lg  D_maxstrenghtN %Lg D_maxstrenght %Lg norm %s \n", N_samp, N, alpha, strenghtN, D_maxstrenght, NORM_TYPE);

	//marco
	long double ave_stability_sampled, min_stability_sampled, max_stability_sampled;
	double **stability;
	double **ave_stability;
	double **max_stability;
	double **min_stability;
	stability = (double **)malloc(P * sizeof(double *));
	for (i = 0; i < P; i++)
	{
		stability[i] = (double *)malloc(N * sizeof(double));
	}
	ave_stability = (double **)malloc(N_samp * sizeof(double *));
	for (i = 0; i < N_samp; i++)
	{
		ave_stability[i] = (double *)malloc(((int)(D_max / delta_D)+1) * sizeof(double));
	}
	max_stability = (double **)malloc(N_samp * sizeof(double *));
	for (i = 0; i < N_samp; i++)
	{
		max_stability[i] = (double *)malloc(((int)(D_max / delta_D)+1) * sizeof(double));
	}
	min_stability = (double **)malloc(N_samp * sizeof(double *));
	for (i = 0; i < N_samp; i++)
	{
		min_stability[i] = (double *)malloc(((int)(D_max / delta_D)+1) * sizeof(double));
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

	for (i = 0; i < N_samp; i++)
	{ //Repetition of the run over N_samp realizations of disorder

		printf("Sample #%d of %d\n", i + 1, N_samp);
		generate_csi(P, csi);
		generate_J(P, csi, J);

		int t = 0;
		for (D = 0; D < D_max; D++)
		{ //Dreaming..

			if (D > 0)
			{
				sigma = generate_rand_initial();
				async_dynamics(sigma, J);
				updateJ(J, sigma, strenght);
				normalizeJ(NORM_TYPE,J);
			}

			if (D % delta_D == 0) //POOP
			{ //Check and measure
               // printf("D_max %d D_maxstrenght %Lg strenght %Lg deltaD %d step%d \n", D_max, D_maxstrenght, strenght, delta_D, t);
               // printf("%d D_max %d D_maxstrenght %Lg strenght %Lg deltaD %d totsteps %d step%d \n", (int)(1/0.001), D_max, D_maxstrenght, strenght, delta_D, (int)(D_max / delta_D), t);
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

				//marco: computing stabilities
				for (int mu = 0; mu < P; mu++)
				{
					for (int i = 0; i < N; i++)
					{
						for (int j = 0; j < N; j++)
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

				for (int mu = 0; mu < P; mu++)
				{
					for (int j = 0; j < N; j++)
					{
						ave_stability[i][t] += stability[mu][j];
						if (max_stability[i][t] < stability[mu][j])
						{
							max_stability[i][t] = stability[mu][j];
						}
						if (min_stability[i][t] > stability[mu][j])
						{
							min_stability[i][t] = stability[mu][j];
						}
					}
				}

				ave_stability[i][t] /= (N * P);

				sigma_new = generate_initial(csi, initial_pattern);
				async_dynamics(sigma_new, J);
				over[i][t] = overlap(csi, sigma_new, initial_pattern);
				t++;
			}
		}
		printf("end sample \n");
	}

	// Analysis of the measures over the different samples

	//marco: analize stabilities
	for (t = 0; t < (int)(D_max / delta_D); t++)
	{
		ave_stability_sampled = 0;
		max_stability_sampled = 0;
		min_stability_sampled = 0;

		for (int i = 0; i < N_samp; i++)
		{
			ave_stability_sampled += ave_stability[i][t];
			max_stability_sampled += max_stability[i][t];
			min_stability_sampled += min_stability[i][t];
		}

		ave_stability_sampled /= N_samp;
		max_stability_sampled /= N_samp;
		min_stability_sampled /= N_samp;

		fprintf(fout3, "Samples %d N %d alpha %Lg strenghtN %Lg D_maxstrenght %Lg dream %d ave_stability %Lg max_stability_sampled %Lg min_stability_sampled %Lg \n", N_samp, N, alpha, strenghtN, D_maxstrenght, t * delta_D, ave_stability_sampled, max_stability_sampled, min_stability_sampled);
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
		fprintf(fout1, "%d\t%lf\t%lf\n", t * delta_D, m, sqrt((sigma_m - m * m) / (double)N_samp));
		fprintf(fout2, "%d\t%lf\t%lf\n", t * delta_D, J_Av, J_Sigma);
		fflush(fout1);
		fflush(fout2);
	}

	fclose(fout1);
	fclose(fout2);
	fclose(fout3);

	return 0;
}
