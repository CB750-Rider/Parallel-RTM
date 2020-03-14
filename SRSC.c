/* Serial Radiative Spectrum Calculator
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
    Copyright (C) 2017 Dr. Lee R. Burchett

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

For more information contact Dr. Lee R. Burchett at 
lee.r.burchett@gmail.com

 */

#include "SRSC.h"
#include "weatherConversion.h"
#include <time.h>
#include <float.h>

#define N_WAVELENGTH 1000
#define P_WEATHER_POINTS 200
#define LINE_WINDOW_WIDTH 25.0

static void calc_index_spectrum(HITRAN_DATA *H,
								WEATHER_DATA *W,
								ftype *nu,
								ctype *spectrum,
								ulong N);
static void calc_index_spectrum_R2(HITRAN_DATA *H,
								WEATHER_DATA *W,
								ftype *nu,
								ctype *spectrum,
								ulong N);
static void calc_index_spectrum_R3(HITRAN_DATA *H,
								WEATHER_DATA *W,
								ftype *nu,
								ctype *spectrum,
								ulong N);
static WEATHER_CONVERSION_ERROR set_input_data( WEATHER_DATA *W, ftype *nu, ulong N);
static void print_spectrum(ctype *spectrum,ftype *nu, ulong N_nu,WEATHER_DATA *wthr);
static void tic(time_t *t);
static void toc(time_t *t);
void write_header(FILE *out, WEATHER_DATA *wthr, ftype *nu, ulong N_nu);
static void write_line(FILE *out, ftype *x, ulong Nx,const char *fmt);
static void read_data_into_vector(FILE *fp, WEATHER_CONVERSION_VECTOR *WV, ftype *Z);

static time_t now[2];

int main(int argc, char *argv[]){
	ulong N=N_WAVELENGTH; /* The number of wavelengths */
	ulong P=P_WEATHER_POINTS; /* The number of points in weather space */
	HITRAN_DATA H;
	WEATHER_DATA wthr;
	ftype *nu = malloc(N_WAVELENGTH*sizeof(ftype));
	ctype *spectrum = malloc(N_WAVELENGTH*P_WEATHER_POINTS*sizeof(ctype));
	HITRAN_PARSE_ERROR rv;
	char function_name[] = "SCSC.c:main()";

	/* Open the data */
	check(LoadHITRANparData(&H));

	/* Open the weather data */
	open_weather_data(P,&wthr);

	set_input_data(&wthr,nu,N);

/*
	printf("Running Version 1.\n");
	tic(now);
	calc_index_spectrum(&H,&wthr,nu,spectrum,N);
	toc(now);printf("\n");

	printf("Running Version 2.\n");
	tic(now);
	calc_index_spectrum_R2(&H,&wthr,nu,spectrum,N);
	toc(now);printf("\n");*/

	printf("Running Version 3.\n");
	tic(now);
	calc_index_spectrum_R3(&H,&wthr,nu,spectrum,N);
	toc(now);printf("\n");

	printf("Done\n");

	print_spectrum(spectrum,nu,N_WAVELENGTH,&wthr);
	free_HITRAN(&H);
	free_weather(&wthr);
	free(spectrum);
	free(nu);
	return 0;
}
WEATHER_CONVERSION_ERROR set_input_data(WEATHER_DATA *W, ftype *nu, ulong N){
	ulong i;
	FILE *fp = fopen("WxProfile.csv","r");
	ftype dP = 200.0/(ftype)W->n_points;
	ftype dT = 100.0/(ftype)W->n_points;
	WEATHER_CONVERSION_VECTOR WV;
	ftype *Z;

	if(fp==NULL){
		for(i=0;i<W->n_points;i++){
			/* Fake some data */
			W->P[i] = 800.0 + dP*(ftype)i;
			W->T[i] = 250.0 + dT*(ftype)i;
			W->mmr[i] = 0.001;
			W->rho[i] = 2.0;
			W->w_rho[i] = 0.2;
			W->ev[i] = 2.5;
		}
	}
	else{
		_wcCheck(openWeatherConversionVector(&WV,(unsigned int)W->n_points));
		Z = malloc(sizeof(ftype)*W->n_points);
		read_data_into_vector(fp,&WV,Z);
		_wcCheck(WV.f.run_conversion(&WV));
		/* TODO Add better integration code here */
		for(i=0;i<W->n_points;i++){
			W->P[i] = WV.val[_PRESSURE][i]/WV.standardPressure;
			W->T[i] = WV.val[_TEMPERATURE_K][i];
			W->ev[i] = WV.val[_VAPOR_PRESSURE][i]/WV.standardPressure;
			W->mmr[i] = WV.val[_MOLE_MIXING_RATIO][i];
			if(i==W->n_points-1)
				W->w_rho[i] = 0.0;
			else
				W->w_rho[i] = (Z[i+1]-Z[i])*WV.val[_ABSOLUTE_HUMIDITY][i]/18.01528;
			if(W->w_rho[i]==0)
				W->rho[i]=0.0;
			else if((log2(W->w_rho[i])-log2(WV.val[_MOLE_MIXING_RATIO][i]))<DBL_MANT_DIG-2)
				W->rho[i] = W->w_rho[i]/WV.val[_MOLE_MIXING_RATIO][i];
			else
				W->rho[i]=0.0;
		}
	}

	ftype d_nu = 950.0/((ftype)N_WAVELENGTH);
	ftype nu_0 = 25.0;
	for(i=0;i<N;i++)
		nu[i] = nu_0 + d_nu*(ftype)i;

	_wcCheck(WV.f.free(&WV));
	free(Z);
	return WEATHER_CONVERSION_SUCCESS;
}
void read_data_into_vector(FILE *fp, WEATHER_CONVERSION_VECTOR *WV, ftype *Z){
	char *line=NULL;
	size_t n=0;
	ulong i=0;
	ftype *P,*T,*vp,*ah;

	P = malloc(sizeof(ftype)*WV->N);
	T = malloc(sizeof(ftype)*WV->N);
	vp = malloc(sizeof(ftype)*WV->N);
	ah = malloc(sizeof(ftype)*WV->N);

	rewind(fp);
	getline(&line,&n,fp);/* Skip the header */
	while(!feof(fp)){
		if(getline(&line,&n,fp)>0){
			if(sscanf(line,"%lg, %lg, %lg, %lg, %lg\n",Z+i,P+i,T+i,vp+i,ah+i)!=5){
				fprintf(stderr,"Error reading in Wx data.\n");
				exit(-1);
			}
			i++;
		}
	}
	WV->f.set_field(WV,P,_PRESSURE);
	WV->f.set_field(WV,T,_TEMPERATURE_C);
	WV->f.set_field(WV,vp,_VAPOR_PRESSURE);
	WV->f.run_conversion(WV);
	if(line!=NULL)free(line);
	free(P);
	free(T);
	free(vp);
	free(ah);
}
void calc_index_spectrum(HITRAN_DATA *H,
								WEATHER_DATA *wthr,
								ftype *nu,
								ctype *spectrum,
								ulong N){
	ulong i,j,k,out_idx;

	for(k=0;k<wthr->n_points;k++){
		for(i=0;i<N;i++){
			out_idx = i*wthr->n_points + k;
			spectrum[out_idx] = continuum(wthr,nu[i],k);
			for(j=0;j<H->n_lines;j++){
				if(fabs(nu[i] - H->nuij[j])<25.0){
					spectrum[out_idx] += calc_line_at(H,wthr,nu[i],j,k,LORENTZIAN);
				} /* If we are close enough */
			} /* For each line in the database */
		} /* For each wavelength (row) in the output */
	} /* For each weather data (column) in the output */
}
void calc_index_spectrum_R2(HITRAN_DATA *H,
								WEATHER_DATA *wthr,
								ftype *nu,
								ctype *spectrum,
								ulong N){
	/* Here we try to speed things up by removing the repeated if statement.
	 * We will assume that all wavelengths and line positions are monotonic
	 * increasing. This way we really only need to keep track of where the
	 * last set started and stopped, then we can reduce the number of "if"
	 * clauses. */
	ulong i,j,k,out_idx,start=0,stop=0;

	for(k=0;k<wthr->n_points;k++){
		for(i=0;i<N;i++){
			out_idx = i*wthr->n_points + k;
			spectrum[out_idx] = continuum(wthr,nu[i],k);
			/* Starting where the last wavelength left off, keep going until
			 * we are in bounds */
			for(j=start;
					(j<H->n_lines) && (nu[i] - H->nuij[j]>25.0);
					j++);
			/* Run all the way until the last place we stopped */
			for(;j<stop;j++)
					spectrum[out_idx] += calc_line_at(H,wthr,nu[i],j,k,LORENTZIAN);
			for(;(j<H->n_lines) && (H->nuij[j] - nu[i]<25.0);j++){
				spectrum[out_idx] += calc_line_at(H,wthr,nu[i],j,k,LORENTZIAN);
				stop=j;
			} /* For each line in the database */
		} /* For each wavelength (row) in the output */
	} /* For each weather data (column) in the output */
}
void calc_index_spectrum_R3(HITRAN_DATA *H,
								WEATHER_DATA *wthr,
								ftype *nu,
								ctype *spectrum,
								ulong N){
	/* Here we try to speed things up by removing the repeated if statement,
	 * and by running the N wavelengths in the innermost loop.
	 * We will assume that all wavelengths and line positions are monotonic
	 * increasing. */
	ulong i,j,k,out_idx,LI[2],WL[2];
	LINE_SHAPE_PRE_CALC lspc;

	/* Find the first line number we need */
	LI[START] = find_index_after(nu[0]-LINE_WINDOW_WIDTH,H->nuij,H->n_lines);
	if(LI[START] >= H->n_lines)LI[START]--;
	LI[STOP] = find_index_after(nu[N-1]+LINE_WINDOW_WIDTH,H->nuij,H->n_lines)-1;

	/*
	printf("The lines used are %lu to %lu with wavenumbers %lf to %lf\n",
			LI[START],LI[STOP],H->nuij[LI[START]],H->nuij[LI[STOP]]);
	printf("The wavenumbers used go from %lf to %lf",nu[0],nu[N-1]);*/

	for(k=0;k<wthr->n_points;k++){
		for(i=0;i<N;i++){
			out_idx = i*wthr->n_points + k;
			spectrum[out_idx] = continuum(wthr,nu[i],k);
		}
		for(j=LI[START];j<LI[STOP];j++){
			pre_calc(H,wthr,&lspc,j,k,LORENTZIAN);
			WL[START] = find_index_after(H->nuij[j]-LINE_WINDOW_WIDTH,nu,N);
			WL[START] = WL[START] > 0 ? WL[START]-1:0;
			WL[STOP] = find_index_after(H->nuij[j]+LINE_WINDOW_WIDTH,nu,N);
			WL[STOP] = WL[STOP] < N ? WL[STOP] : N-1;
			for(i=WL[START];i<WL[STOP];i++){
				out_idx = i*wthr->n_points + k;
				spectrum[out_idx] += calc_line_final(H,wthr,nu[i],j,k,&lspc,LORENTZIAN);
			}
		}
	}

}

void print_spectrum(ctype *spectrum,ftype *nu, ulong N_nu,WEATHER_DATA *wthr){
	FILE *out = fopen("OutputSpectrum.csv","w");
	ulong i,j,P=wthr->n_points;

	write_header(out,wthr,nu,N_nu);

	for(i=0;i<N_nu;i++){
		fprintf(out,"%g+i%g",creal(spectrum[i*P]),cimag(spectrum[i*P]));
		for(j=1;j<P;j++)
			fprintf(out,",%g+i%g",creal(spectrum[i*P + j]),cimag(spectrum[i*P + j]));
		fprintf(out,"\n");
	}
	fclose(out);
}
void write_header(FILE *out, WEATHER_DATA *wthr, ftype *nu, ulong N_nu){
	fprintf(out,"#Temperature(K),Pressure(atm),mole mixing ratio(ppm),vapor pressure(atm),moist air rho(N/m^2),water vapor rho(N/m^2),nu(cm^-1),data(%lu X %lu)\n",N_nu,wthr->n_points);

	write_line(out,wthr->T,wthr->n_points,"%lf");
	write_line(out,wthr->P,wthr->n_points,"%lf");
	write_line(out,wthr->mmr,wthr->n_points,"%g");
	write_line(out,wthr->ev,wthr->n_points,"%g");
	write_line(out,wthr->rho,wthr->n_points,"%g");
	write_line(out,wthr->w_rho,wthr->n_points,"%g");
	write_line(out,nu,N_nu,"%lf");

}
void write_line(FILE *out, ftype *x, ulong Nx,const char *fmt){
	ulong i;
	fprintf(out,fmt,x[0]);
	for(i=1;i<Nx;i++){
		fprintf(out,",");
		fprintf(out,fmt,x[i]);
	}
	fprintf(out,"\n");


}
static void tic(time_t *t){
	time(t);
}
static void toc(time_t *t){
	int days, hours, minutes;
	ftype dt;
	time(t+1);
	dt = difftime(t[1],t[0]);
	days = (int)floor(dt/86400.0);
	dt = dt-86400.0*(ftype)days;
	hours  = (int)floor(dt/3600.0);
	dt = dt-3600.0*(ftype)hours;
	minutes = (int)floor(dt/60.0);
	dt = dt - 60.0*(ftype)minutes;
	printf("%d:%2d:%2d:%lf",days,hours,minutes,dt);
}

