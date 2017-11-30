/*
 * SRSC.c
 *
 *  Created on: Oct 26, 2017
 *      Author: Dr. Lee Burchett
 *
 *      This program is meant to explore a serial
 *      implementation of the banded matrix
 *      calculation of refractive index spectra
 *      from the HITRAN database
 *
 */

#include "SRSC.h"


#define N_WAVELENGTH 15

static void calc_index_spectrum(HITRAN_DATA *H,
								WEATHER_DATA *W,
								double *nu,
								double complex *spectrum,
								uint N);
int main(int argc, char *argv[]){
	uint i;
	uint M=10; /* The number of moleculular absorption lines */
	uint N=N_WAVELENGTH; /* The number of wavelengths */
	uint W=1; /* The number of points in weather space */
	HITRAN_DATA H;
	WEATHER_DATA wthr;
	double nu[N_WAVELENGTH];
	double complex spectrum[N_WAVELENGTH];
	open_HITRAN_data(M,&H);

	for(i=0;i<N;i++) nu[i] = (double)i;

	calc_index_spectrum(&H,&wthr,nu,spectrum,N);

	free_HITRAN(&H);
	return 0;
}

void calc_index_spectrum(HITRAN_DATA *H,
								WEATHER_DATA *wthr,
								double *nu,
								double complex *spectrum,
								uint N){
	uint i,j;

	for(i=0;i<N;i++){
		spectrum[i] = continuum(wthr,nu[0]);
		for(j=0;j<H->n_lines;j++){
			if(((j-i)<3) && ((i-j)<3)){
				spectrum[i] += calc_line_at(H,wthr,nu[j],j,LORENTZIAN);
			}
		}
	}
}
void open_HITRAN_data(uint M, HITRAN_DATA *h){
	h->n_lines = M;
	h->Sij = malloc(M*sizeof(double));
	h->nuij = malloc(M*sizeof(double));
	h->Aij = malloc(M*sizeof(double));
	h->gamma_air = malloc(M*sizeof(double));
	h->gamma_self = malloc(M*sizeof(double));
	h->Eij = malloc(M*sizeof(double));
	h->n_air = malloc(M*sizeof(double));
	h->delta_air = malloc(M*sizeof(double));
	h->g_p = malloc(M*sizeof(double));
	h->g_pp = malloc(M*sizeof(double));
	h->Q_ratio = malloc(M*sizeof(double));
	h->molec = malloc(M*sizeof(uint));
	h->isotopologue = malloc(M*sizeof(uint));
}
void free_HITRAN(HITRAN_DATA *h){
	if(h->n_lines==0) return;
	h->n_lines = 0;
	free(h->Sij);
	free(h->nuij);
	free(h->Aij);
	free(h->isotopologue);
	free(h->molec);
	free(h->gamma_self);
	free(h->gamma_air);
	free(h->Eij);
	free(h->n_air);
	free(h->delta_air);
	free(h->g_p);
	free(h->g_pp);
	free(h->Q_ratio);

}
double complex calc_line_at(HITRAN_DATA *H,
							WEATHER_DATA *W,
							double nu,
							uint j,
							LINE_SHAPE_FUNCTION fcn)
{
	return (double complex)j + 0*I;
}
double continuum(WEATHER_DATA *w, double nu){
	return 0.01*I;
}
