#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

void generate_csi(int P, int **csi);
void generate_J(int P, int **csi, double **J);
double H(double **J, int *sigma);
void async_dynamics(int *sigma, double **J);
void updateJ(double **J, int *sigma, double strenght);
void normalizeJ(char *type_of_norm, double **J);
long double get_norm(double **J, int i);
void generate_rand_initial(int *sigma);
void my_init();

int D_max, max_unlearning_samp, max_iter, seed;
int N, P, N_samp, max_iter, seed, delta_seed;
long double alpha, lambda, c, strenght;
int *sigma, **csi, **mask;
double **J;
char string[200];
FILE *fout1;

int main(int argc, char *argv[])
{
    if (argc != 10)
    {
        printf("Please use 9 input parameters \n");
        printf("Usage: ./perceptronVSunlearning_trajectory.exe N alpha lambda max_stability D_max n_samples unlearning_strenghtN n_unlearning_samples delta_seed\n");
        exit(-9);
    }

    N = atoi(argv[1]);
    alpha = atof(argv[2]);
    lambda = atof(argv[3]);
    c = atof(argv[4]);
    D_max = atoi(argv[5]);
    N_samp = atoi(argv[6]);
    strenght = atof(argv[7]);
    max_unlearning_samp = atoi(argv[8]);
    int delta_seed = atoi(argv[9]);
    seed = time(0) + delta_seed;
    srand48(seed);
    my_init();

    for (int samp = 0; samp < N_samp; samp++)
    {
        //perceptron part
        generate_csi(P, csi);
        generate_J(P, csi, J);

        //printing initial J matrix
        fprintf(fout1, "sym_prc samp %d correspondance 0 step 0 J ", samp); fflush(fout1);
        for (int i = 0; i < N; i++)
        {
            for (int j = i + 1; j < N; j++)
            {

                fprintf(fout1, "%lf ", J[i][j]); fflush(fout1);
            }
        }
        fprintf(fout1, "\n"); fflush(fout1);


        int endcounter, t;
        int correspondance=1;
        for (t = 1; t < max_iter; t++)
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
            //printing J matrix
            if (t % 10 == 0 && endcounter != N * P)
            {
                fprintf(fout1, "sym_prc samp %d correspondance %d step %d J ", samp, correspondance, t); fflush(fout1);
                for (int i = 0; i < N; i++)
                {
                    for (int j = i + 1; j < N; j++)
                    {

                        fprintf(fout1, "%lf ", J[i][j]); fflush(fout1);
                    }
                }
                fprintf(fout1, "\n"); fflush(fout1);
                correspondance++;
            }

            if (endcounter == N * P)
            {
                //printing J matrix
                fprintf(fout1, "sym_prc samp %d correspondance %d step %d J ", samp, correspondance, t); fflush(fout1);
                for (int i = 0; i < N; i++)
                {
                    for (int j = i + 1; j < N; j++)
                    {

                        fprintf(fout1, " %lf ", J[i][j]); fflush(fout1);
                    }
                }
                fprintf(fout1, "\n"); fflush(fout1);
                fprintf(fout1, "sym_prc samp %d max_correspondance %d tot_steps %d \n", samp, correspondance, t); fflush(fout1);
                break;
            }
            if (t == max_iter - 1)
            {
                fprintf(fout1, "sym_perc algorithm did not converge\n");
                break;
            }
        }

        //unlearning part
        for (int unlearning_samp = 0; unlearning_samp < max_unlearning_samp; unlearning_samp++)
        {
            generate_J(P, csi, J);
            //printing initial J matrix
            fprintf(fout1, "unlearning samp %d unlearning_samp %d correspondance 0 step 0 of %d J ", samp, unlearning_samp, D_max); fflush(fout1);
            for (int i = 0; i < N; i++)
            {
                for (int j = i + 1; j < N; j++)
                {

                    fprintf(fout1, " %lf ", J[i][j]); fflush(fout1);
                }
            }
            fprintf(fout1, "\n"); fflush(fout1);

            correspondance=1;
            for (int D = 0; D < D_max; D++)
            { //Dreaming..
                generate_rand_initial(sigma);
                async_dynamics(sigma, J);
                updateJ(J, sigma, strenght);
                //normalizeJ(NORM_TYPE, J); we use NO_NORM anyways

                if (D % (int)(10 * D_max / t) == 0 && D != 0 && D != D_max)
                {
                    fprintf(fout1, "unlearning samp %d unlearning_samp %d correspondance %d step %d of %d J ", samp, unlearning_samp, correspondance, D, D_max); fflush(fout1);

                    for (int i = 0; i < N; i++)
                    {
                        for (int j = i + 1; j < N; j++)
                        {
                            fprintf(fout1, "%lf ", J[i][j]); fflush(fout1);
                        }
                    }
                    fprintf(fout1, "\n"); fflush(fout1);
                    correspondance++;
                }
                if (D == D_max-1)
                {
                    fprintf(fout1, "unlearning samp %d unlearning_samp %d correspondance %d step %d of %d J ", samp, unlearning_samp, correspondance, D, D_max-1); fflush(fout1);
                    for (int i = 0; i < N; i++)
                    {
                        for (int j = i + 1; j < N; j++)
                        {
                            fprintf(fout1, "%lf ", J[i][j]); fflush(fout1);
                        }
                    }
                    fprintf(fout1, "\n"); fflush(fout1);
                }
            }
        }
    }
    fclose(fout1);
    return (0);
}

/******************************************/
void my_init()
{
    max_iter = 1000;
    P = (int)(alpha * (double)N);
    sigma = (int *)malloc(N * sizeof(int));
    csi = (int **)malloc(P * sizeof(int *));
    mask = (int **)malloc(P * sizeof(int *));
    J = (double **)malloc(N * sizeof(double *));
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

    sprintf(string, "perceptronVSunlearning_trajectory_N%d_alpha%Lg_lambda%Lg_prcstab%Lg_Nsamp%d_unlearningSamp%d_seed%d.dat", N, alpha, lambda, c, N_samp, max_unlearning_samp, seed);
    fout1 = fopen(string, "w");
    if (fout1 == NULL)
    {
        printf("couldn't open fout1\n");
        exit(9);
    }
}

/******************************************/
void generate_csi(int P, int **csi)
{ //generate memories as a bernoulli process with probability p
    int i, j;
    for (i = 0; i < P; i++)
    {
        for (j = 0; j < N; j++)
        {
            if ((lrand48() / (double)RAND_MAX) < 0.5)
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

/******************************************/
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

/******************************************/
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

/******************************************/
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
        }
    }
}

/******************************************/
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

/******************************************/
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

/******************************************/
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

/******************************************/
void generate_rand_initial(int *sigma)
{ //random shooting generator
    int i;
    for (i = 0; i < N; i++)
    {
        if ((lrand48() / (double)RAND_MAX) < 0.5)
        {
            sigma[i] = 1;
        }
        else
        {
            sigma[i] = -1;
        }
    }
}
