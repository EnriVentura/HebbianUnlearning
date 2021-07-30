#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

int N, P, N_samp, max_iter;
long double lambda, alpha, c, p = 0.5, p_in = 0;

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
		if (time == 99999){printf("ABORT async_dynamics did not converge\n"); exit (-9);}
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


int main(int argc, char *argv[]){

	if (argc != 7)
	{
		printf("Please use 6 input parameters \n");
		printf("Usage: ./get_J_sym_perceptron.exe N alpha lambda max_stability max_number_of_iterations delta_seed\n");
		exit(-9);
	}

	N = atoi(argv[1]);
	alpha = atof(argv[2]);
	lambda = atof(argv[3]);
	c = atof(argv[4]);
	max_iter = atoi(argv[5]);
	int delta_seed = atoi(argv[6]);
	int seed = time(0) + delta_seed;
	srand48(seed);

	FILE *mat;
	FILE *patt;
	char string_mat[150], string_patt[150];
	sprintf(string_mat, "J_matrix_N%d_alpha%Lg.dat", N, alpha);
	sprintf(string_patt, "J_patts_N%d_alpha%Lg.dat", N, alpha);
	mat = fopen(string_mat, "r");
	patt = fopen(string_patt, "r");
	

	P = (int)(alpha*(double)N);

	int **csi, *sigma, *sigma_new, **mask;
	double **J;

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

	long double ave_stability_sampled, ave_stability2_sampled, min_stability_sampled, min_stability2_sampled, max_stability_sampled, max_stability2_sampled, asymmetry_sampled, norm, n_sat_sampled, n_sat2_sampled;
	double **ave_stability;
	double **max_stability;
	double **min_stability;
	double **asymm;
	double **n_sat;
	double **over;
	double **J_av, **J_sigma;
	int *i_neg;

	i_neg = (int *)malloc(N * sizeof(int));

	for(int i = 0; i < N; i++){
		i_neg[i] = -1;
	}

	int N_samp_over = 10;
	over = (double **)malloc(N_samp * sizeof(double));

	J_av = (double **)malloc(N_samp * sizeof(double));
	J_sigma = (double **)malloc(N_samp * sizeof(double));

	n_sat = (double **)malloc(N_samp * sizeof(double *));
	for (int i = 0; i < N_samp; i++)
	{
		n_sat[i] = (double *)malloc(max_iter * sizeof(double));
		J_av[i] = (double *)malloc(max_iter * sizeof(double));
		J_sigma[i] = (double *)malloc(max_iter * sizeof(double));
		over[i] = (double *)malloc(max_iter * sizeof(double));
	}

	asymm = (double **)malloc(N_samp * sizeof(double *));
	for (int i = 0; i < N_samp; i++)
	{
		asymm[i] = (double *)malloc(max_iter * sizeof(double));
	}
	
	ave_stability = (double **)malloc(N_samp * sizeof(double *));
	for (int i = 0; i < N_samp; i++)
	{
		ave_stability[i] = (double *)malloc(max_iter * sizeof(double));
	}
	max_stability = (double **)malloc(N_samp * sizeof(double *));
	for (int i = 0; i < N_samp; i++)
	{
		max_stability[i] = (double *)malloc(max_iter * sizeof(double));
	}
	min_stability = (double **)malloc(N_samp * sizeof(double *));
	for (int i = 0; i < N_samp; i++)
	{
		min_stability[i] = (double *)malloc(max_iter * sizeof(double));
	}

	if (n_sat == NULL)
	{
		printf("malloc of n_sat failed.\n");
		exit(EXIT_FAILURE);
	}

	if (asymm == NULL)
	{
		printf("malloc of asymmetry failed.\n");
		exit(EXIT_FAILURE);
	}

	if (ave_stability == NULL)
	{
		printf("malloc of ave_stability failed.\n");
		exit(EXIT_FAILURE);
	}
	if (max_stability == NULL)
	{
		printf("malloc of max_stability failed.\n");
		exit(EXIT_FAILURE);
	}
	if (min_stability == NULL)
	{
		printf("malloc of min_stability failed.\n");
		exit(EXIT_FAILURE);
	}

	int initial_pattern = 0;
	double max_stab;
	double stab;

	N_samp = 1;
	for(int samp = 0; samp < N_samp; samp++){
		
		printf("Sample #%d of %d\n", samp + 1, N_samp);
		double coup;
		for(int ii = 0; ii < N; ii++){
			for(int jj = 0; jj < N; jj++){
				fscanf(mat, "%lf\t", &coup);
				J[ii][jj] = coup;
			}
			fscanf(mat, "\n");
		}
		int pattt;
		for(int ii = 0; ii < P; ii++){
			for(int jj = 0; jj < N; jj++){
				fscanf(patt, "%d\t", &pattt);
				csi[ii][jj] = pattt;
			}
			fscanf(patt, "\n");
		}

		double k;
		int t = 0;

		do{


			for (int j = 0; j < P; j++)
			{
				int ii = 0;
				for (int i = 0; i < N; i++)
				{
					mask[j][i] = 0;
					stab = 0;
					norm = get_norm(J, i);
					for (int kk = 0; kk < N; kk++)
					{
						stab += J[i][kk] * (double)csi[j][kk] * csi[j][i] / norm;
					}


					if(stab >= c){
						mask[j][i] = 0;
						
					}else{
						max_stab = stab;
						mask[j][i] = 1;
					}
					
		
				}
				
			}

			for(int i = 0; i < N; i++){
				for(int j = i+1; j < N; j++){
					for(int k = 0; k < P; k++){
						J[i][j] += lambda * (mask[k][i] + mask[k][j]) * csi[k][i] * csi[k][j];
					}
					J[j][i] = J[i][j];
				}
				J[i][i] = 0;
				
			}

			if(t == max_iter-1){
				printf("max stab = %lf\n", max_stab);
				char string[150];
	
				sprintf(string, "sym_perceptron_finalJ_N%d_alpha%Lg_lambda%Lg_maxstab%Lg.dat", N, alpha, lambda, c);

				FILE *fout1;
				fout1 = fopen(string, "w");

				for(int i = 0; i < N; i++){
					for(int j = i + 1; j < N; j++){

						fprintf(fout1, "%lf\n", J[i][j]);

					}
				}
				fclose(fout1);
			}

			t++;

		}while(t < max_iter);


	}



}