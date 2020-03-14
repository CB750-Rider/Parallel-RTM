/*  Serial Radiative Spectrum Calculator
 * SRSC.h
 *
 *  Created on: Oct 26, 2017
 *      Author: Dr. Lee Burchett
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

#ifndef SRSC_H_
#define SRSC_H_
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <float.h>
#include "HITRAN_Parser.h"
#include "weatherConversion.h"

#define N_HITRAN_PAR_FILES 1

typedef enum{
	LORENTZIAN,
	GAUSSIAN,
	VOIGT,
	N_LINE_SHAPES
} LINE_SHAPE_FUNCTION;

/* Define as needed
typedef float ftype;
typedef float complex ctype;
#define _ftype_max FLT_MAX
 * */
typedef double ftype;
typedef double complex ctype;
#define _ftype_max DBL_MAX

/* TODO:
 * 1. Define constants
 * 2. Create a HITRAN Parser to HITRAN_DATA converter
 * 		a. Should sort the HITRAN data by wavenumber
 * 		b. Should open HITRAN_DATA and copy over
 * 3. Create a function to calculate all the Q_ratios
 * 4. Write the various line shape functions */

typedef struct {
	ulong n_lines;
	ftype *Sij; /* Spectral Line Intensity */
	ftype *nuij; /* Spectral Transition Wavenumber */
	ftype *Aij; /* Einstein A coefficient of the transition */
	ftype *gamma_air; /* The air broadened half width */
	ftype *gamma_self; /* The self-broadened half width */
	ftype *Eij; /* The lower state energy of the transition */
	ftype *n_air; /* the coefficient of the temperature dependence of the air broadened half width */
	ftype *delta_air; /* Pressure shift at T ref = 296 K and P ref = 1 atm*/
	ftype *g_p; /* lower and upper statistical weights */
	ftype *g_pp;
	ftype *Q_ratio; /* Q(T_ref)/Q(T) Partition sum ratio */
	unsigned short *molec;
	unsigned short *isotopologue;
	ftype *mmr; /* Standard mole mixing ratio */
	HITRAN_FILE *data; /* Pointer to the HITRAN par data */
} HITRAN_DATA;

typedef struct {
	ulong n_points;
	ftype *T; /* Temperature in Kelvin */
	ftype *P; /* Pressure in atmospheres */
	ftype *mmr; /* Mole mixing ratio in ppm */
	ftype *ev; /* Partial pressure of water vapor in atmospheres */
	ftype *rho; /* atmosphere number density in number per meter squared */
	ftype *w_rho; /* water vapor number density in number per meter squared */
} WEATHER_DATA;

typedef struct {
	ftype S; /* Line intensity */
	ftype gamma; /* Line Width */
	ftype gamma2;
	ftype nu0; /* Line Center */
	ftype nu02;
	ftype rhoS_pi; /* Number Density*S/pi */
} LINE_SHAPE_PRE_CALC;

typedef enum {
	START,
	STOP
} VECTOR_END_NAMES;

/* DataManagement.c */
void open_HITRAN_data(ulong M, HITRAN_DATA *h);
void free_HITRAN(HITRAN_DATA *h);
void open_weather_data(ulong P, WEATHER_DATA *w);
void free_weather(WEATHER_DATA *w);
ulong find_index_after(ftype x0, ftype *x, ulong N);

/* FileIO.c */
HITRAN_PARSE_ERROR LoadHITRANparData(HITRAN_DATA *H);

/* SRSC.c */
ctype calc_line_at(HITRAN_DATA *H,
					WEATHER_DATA *W,
					ftype nu,
					ulong line_idx,
					ulong wx_idx,
					LINE_SHAPE_FUNCTION fcn);
void pre_calc(HITRAN_DATA *H,
		      WEATHER_DATA *W,
			  LINE_SHAPE_PRE_CALC *out,
			  ulong line_idx,
			  ulong wx_idx,
			  LINE_SHAPE_FUNCTION fcn);
ctype calc_line_final(HITRAN_DATA *H,
		WEATHER_DATA *W,
		ftype nu,
		ulong line_idx,
		ulong wx_idx,
		LINE_SHAPE_PRE_CALC *p,
		LINE_SHAPE_FUNCTION *fcn);
ctype continuum(WEATHER_DATA *w, ftype nu,ulong i);
ftype find_Q_ratio(HITRAN_FILE *H,
		           unsigned short molec,
				   unsigned short isotope_num,
				   ftype T);
#endif /* SRSC_H_ */
