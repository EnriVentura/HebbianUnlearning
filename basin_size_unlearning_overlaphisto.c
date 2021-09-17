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
int N, N_samp, delta_seed;
long double alpha, strenghtN, D_maxstrenght;
int *sigma1;
int *sigma2;
int *sigma_new;
double *field;
int delta_seed;

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

	int i, j, k, l, flag = 0, flip = 0, time = 0, count = 0;
	for (l = 0; l < N; l++)
	{
		field[l] = 0;
		for (j = 0; j < N; j++)
		{
			field[l] += J[l][j] * sigma[j];
		}
		if (sigma[l] * field[l] > 0)
		{
			count++;
		}
	}
	if (count < N)
	{
		while (flag == 0 && time < 100000)
		{
			do
			{
				i = (int)((lrand48() / (double)RAND_MAX) * (double)N);
			} while (i == N);

			if (field[i] > 0)
			{
				sigma_new[i] = 1;
			}
			else
			{
				sigma_new[i] = -1;
			}
			flip = sigma[i] * sigma_new[i];
			count = 0;

			if (flip == -1)
			{
				for (j = 0; j < N; j++)
				{
					field[j] = field[j] - 2 * J[j][i] * sigma[i];
					if (sigma[j] * field[j] <= 0 && j != i)
					{
						count++;
					}
				}
				if (count == 0)
				{
					flag++;
				}
			}
			else
			{
				for (j = 0; j < N; j++)
				{
					if (sigma[j] * field[j] <= 0 && j != i)
					{
						count++;
					}
				}
				if (count == 0)
				{
					flag++;
				}
			}
			time++;
			sigma[i] = sigma_new[i];
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
    delta_seed= atoi(argv[6]);
	int seed = time(0)+delta_seed;
	srand48(seed);

	char string1[150];
	sprintf(string1, "unlearning_overlaphisto_N%d_alpha%Lg_strenghtN%Lg_Nsamp%d_seed%d.dat", N, alpha, strenghtN, N_samp, seed);

	FILE *fout1;
	fout1 = fopen(string1, "w");

		int i, j, l, t, on;
		int *sigma, **csi;
		double **J;


		double time1, time2, time_tot = 0, time_tott = 0;
		

		int P = (int)(alpha * (double)N);
		double m, sigma_m;
		double **over;
		int **freq;

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
		D_maxstrenght = Deps_measure[0] + 0.05;
		double strenght = strenghtN / (long double)N;
		int D, D_max = (int)(D_maxstrenght / strenght), D_in = (int)(Deps_measure[0]/strenght);
		// Initialization of the main configuration/interaction arrays

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
		field = (double *)malloc(N * sizeof(double));
		if (field == NULL)
		{
			printf("malloc of field array failed.\n");
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

		freq = (int **)malloc(((int)(0.5 * N) + 1) * sizeof(int *));
		over = (double **)malloc(N_samp * sizeof(double *));
		for (int i = 0; i < N_samp; i++)
		{
			over[i] = (double *)malloc(((int)(0.5 * N) + 1) * sizeof(double));
		}
		for(int i = 0; i < (int)(0.5 * N) + 1; i++){
			freq[i] = (int *)malloc((2 * N + 1) * sizeof(int));
		}

		int initial_pattern = 0; //index of the test pattern with respect to the which we compute the overlap

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
					normalizeJ(NORM_TYPE,J);
				}

				if (D % D_in == 0 && D_in != 0)
				{ 
					// Initialization of the output files

					for (int l = 0; l < (int)(0.5 * N) + 1; l++) //this was l < (int)(0.5*N)+1
					{
						sigma_new = generate_initial(csi, initial_pattern, (1 - l * (double)(2 / (double)N)) * 0.5);
						async_dynamics(sigma_new, J);
						over[samp][l] = overlap(csi, sigma_new, initial_pattern);
					}	

					t++;
				}
			
			
			}
		
		}

		for (int l = 0; l < (int)(0.5 * N) + 1; l++)
		{
			for(int i = 0; i < (int)(2 * N) + 1; i++){
				freq[l][i] = 0;
			}
			double m = 0;
			double sigma_m = 0;
			for (int i = 0; i < N_samp; i++)
			{
				m = m + over[i][l] / (double)N_samp;
				sigma_m = sigma_m + over[i][l] * over[i][l] / (double)N_samp;
				freq[l][(int)((over[i][l] + 1)*(double)N)]++;
				//printf("%d\n", (int)(over[i][l] + 1)*N);
			}
	
		}
			
		for(int j = 0; j < (int)(0.5 * N) + 1; j++){
			fprintf(fout1, "%Lg\t", (long double)(j * 2 / (double)N));
			fflush(fout1);
			for(int i = 0; i < 2*N + 1; i++){
				fprintf(fout1, "%d\t", freq[j][i]);
				fflush(fout1);
			}
			fprintf(fout1, "\n");
			fflush(fout1);
		}

	fclose(fout1);
	return 0;
}
