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

	char string[150], string2[150], string3[150];
	
	sprintf(string, "perceptron_overlap_N%d_alpha%Lg_lambda%Lg_maxstab%Lg_Nsamp%d.dat", N, alpha, lambda, c, N_samp);
	sprintf(string2, "perceptron_Jmoments_N%d_alpha%Lg_lambda%Lg_maxstab%Lg_Nsamp%d.dat", N, alpha, lambda, c, N_samp);
	sprintf(string3, "perceptron_stabilities_N%d_alpha%Lg_lambda%Lg_maxstab%Lg_Nsamp%d.dat", N, alpha, lambda, c, N_samp);

	FILE *fout1;
	fout1 = fopen(string, "w");
	FILE *fout2;
	fout2 = fopen(string2, "w");
	FILE *fout3;
	fout3 = fopen(string3, "w");
	

	P = (int)(alpha*(double)N);

	int **csi, *sigma, *sigma_new;
	double **J;

	sigma = (int *)malloc(N * sizeof(int));
	sigma_new = (int *)malloc(N * sizeof(int));
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
	int *mu_neg;

	mu_neg = (int *)malloc(P * sizeof(int));

	for(int i = 0; i < P; i++){
		mu_neg[i] = -1;
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

	for(int samp = 0; samp < N_samp; samp++){
		
		printf("Sample #%d of %d\n", samp + 1, N_samp);
		generate_csi(P, csi);
		generate_J(P, csi, J);

		double k;
		int t = 0;

		do{

			J_av[samp][t] = 0;
			J_sigma[samp][t] = 0;
			for (int l = 0; l < N; l++)
			{
				for (int j = l + 1; j < N; j++)
				{
					J_av[samp][t] = J_av[samp][t] + J[l][j] / (double)(0.5 * N * (N - 1));
					J_sigma[samp][t] = J_sigma[samp][t] + J[l][j] * J[l][j] / (double)(0.5 * N * (N - 1));
				}
			}
			J_sigma[samp][t] = J_sigma[samp][t] - J_av[samp][t] * J_av[samp][t];

			double stab_max = 0;
			double k = 10000;

			int i_min;
			int mu_min;
			for (int i = 0; i < N; i++)
			{
				int mmu = 0;
				for (int j = 0; j < P; j++)
				{

					double stab = 0;
					double norm = get_norm(J, i);
					for (int kk = 0; kk < N; kk++)
					{
						stab += J[i][kk] * csi[j][kk] * csi[j][i] / norm;
					}

					ave_stability[samp][t] += stab/(double)(N*P);

					if(stab >= c){
						n_sat[samp][t] += 1/((double)(N*P));
					}else{
						mu_neg[mmu] = j;
						mmu++;
					}
					
					if (stab >= stab_max)
					{
						stab_max = stab;
						
					}
					if (stab <= k)
					{
						k = stab;

					}
				}
				mmu = 0;
				while(mu_neg[mmu] != -1 && mmu < N){
					for(int kk = 0; kk < N; kk++){
						J[i][kk] += lambda * csi[mu_neg[mmu]][i] * csi[mu_neg[mmu]][kk];
					}
					mu_neg[mmu] = -1;
					mmu++;
				}
				J[i][i] = 0;
				
			}

			max_stability[samp][t] = stab_max;
			min_stability[samp][t] = k;

			asymm[samp][t] = 0;
			double norm = 0;
			for (int kk = 0; kk < N; kk++)
			{
				for (int j = kk + 1; j < N; j++)
				{
					norm += J[kk][j] * J[kk][j];
					asymm[samp][t] += J[kk][j] * J[j][kk];
				}
			}
			asymm[samp][t] /= norm;

			double Overlap = 0;
			for (int l = 0; l < N_samp_over; l++)
			{
				sigma_new = generate_initial(csi, initial_pattern);
				async_dynamics(sigma_new, J);
				Overlap += overlap(csi, sigma_new, initial_pattern);
			}
			over[samp][t] = Overlap / (double)N_samp_over;

			t++;

		}while(t < max_iter);


	}

	for (int t = 0; t < max_iter; t++)
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

		for (int i = 0; i < N_samp; i++)
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

		fprintf(fout3, "%d\t%d\t%Lg\t%Lg\t%d\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\t%Lg\n", N_samp, N, alpha, lambda, max_iter, (long double)(t*lambda), ave_stability_sampled, (long double)(sqrt((ave_stability2_sampled - ave_stability_sampled * ave_stability_sampled) / (long double)N_samp)), max_stability_sampled, (long double)(sqrt((max_stability2_sampled - max_stability_sampled * max_stability_sampled) / (long double)N_samp)), min_stability_sampled, (long double)(sqrt((min_stability2_sampled - min_stability_sampled * min_stability_sampled) / (long double)N_samp)), asymmetry_sampled, n_sat_sampled, (long double)(sqrt((n_sat2_sampled - n_sat_sampled * n_sat_sampled) / (long double)N_samp)));
		fflush(fout3);
	}

	for (int t = 0; t < max_iter; t++)
	{
		double m = 0;
		double J_Av = 0;
		double sigma_m = 0;
		double J_Sigma = 0;
		for (int i = 0; i < N_samp; i++)
		{
			m = m + over[i][t] / (double)N_samp;
			J_Av = J_Av + J_av[i][t] / (double)N_samp;
			J_Sigma = J_Sigma + J_sigma[i][t] / (double)N_samp;
			sigma_m = sigma_m + over[i][t] * over[i][t] / (double)N_samp;
		}
		fprintf(fout1, "%Lg\t%lf\t%lf\n", (long double)(t*lambda), m, sqrt((sigma_m - m * m) / (double)N_samp));
		fprintf(fout2, "%Lg\t%lf\t%lf\n", (long double)(t*lambda), J_Av, J_Sigma);
		fflush(fout1);
		fflush(fout2);
	}



}