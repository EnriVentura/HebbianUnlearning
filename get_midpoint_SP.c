#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define epsilon 0.01
#define N_samp 100

int main(){

	int N[] = {100, 150, 200, 300, 400, 500};
	long double alpha[] = {0.3, 0.4, 0.5};
	long double c[] = {1.294, 0.987, 0.764};
	long double lambda = 1;

	long double mi;
	double mf, sigma_mf;

	char string2[150];
	char string1[150];

	for(int a = 0; a < 3; a++){

		FILE *fout;
		sprintf(string1, "sym_c%Lg_alpha%Lg_scaling.dat", c[a], alpha[a]);

		fout = fopen(string1, "w");
		
		for(int i = 0; i < 6; i++){
			FILE *fin;
			sprintf(string2, "sym_perceptron_basin_N%d_alpha%Lg_lambda%Lg_maxstab%Lg_Nsamp%d.dat", N[i], alpha[a], lambda, c[a], N_samp);
			fin = fopen(string2, "r");
			double sum = 0, sigma_sum = 0;
			int num = 0;
			for(int l = 0; l < (int)(0.5*N[i]); l++){
				fscanf(fin, "%Lg\t%lf\t%lf\n", &mi, &mf, &sigma_mf);
				
				if(mf < 0.54 && mf >= 0.46){
					sum += mi;
					sigma_sum += mi*mi;
					num++;
				}
			}
			double ave = sum/(double)num;
			fprintf(fout, "%lf\t%lf\t%lf\n", 1/(double)N[i], ave, sqrt((sigma_sum/(double)num - ave*ave)/(double)num));
		}
		fclose(fout);
	}




}