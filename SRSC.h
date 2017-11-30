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
#include "HITRAN_Parser.h"

typedef enum{
	LORENTZIAN,
	GAUSSIAN,
	VOIGT,
	N_LINE_SHAPES
} LINE_SHAPE_FUNCTION;

typedef unsigned int uint;

/* TODO:
 * 1. Define constants
 * 2. Create a HITRAN Parser to HITRAN_DATA converter
 * 		a. Should sort the HITRAN data by wavenumber
 * 		b. Should open HITRAN_DATA and copy over
 * 3. Create a function to calculate all the Q_ratios
 * 4. Write the various line shape functions */

typedef struct {
	uint n_lines;
	double *Sij; /* Spectral Line Intensity */
	double *nuij; /* Spectral Transition Wavenumber */
	double *Aij; /* Einstein A coefficient of the transition */
	double *gamma_air; /* The air broadened half width */
	double *gamma_self; /* The self-broadened half width */
	double *Eij; /* The lower state energy of the transition */
	double *n_air; /* the coefficient of the temperature dependence of the air broadened half width */
	double *delta_air; /* Pressure shift at T ref = 296 K and P ref = 1 atm*/
	double *g_p; /* lower and upper statistical weights */
	double *g_pp;
	double *Q_ratio; /* Q(T_ref)/Q(T) Partition sum ratio */
	unsigned short *molec;
	unsigned short *isotopologue;
} HITRAN_DATA;

typedef struct {
	double T; /* Temperature in Kelvin */
	double P; /* Pressure in atmospheres */
	double mmr; /* Mole mixing ratio in ppm */
	double rho; /* number density in TODO define this*/

} WEATHER_DATA;

void open_HITRAN_data(uint M, HITRAN_DATA *h);
void free_HITRAN(HITRAN_DATA *h);
double complex calc_line_at(HITRAN_DATA *H,
							WEATHER_DATA *W,
							double nu,
							uint j,
							LINE_SHAPE_FUNCTION fcn);
double continuum(WEATHER_DATA *w, double nu);

#endif /* SRSC_H_ */
