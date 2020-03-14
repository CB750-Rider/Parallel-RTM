/*
 * DataManagement.c
 *
 *  Created on: Feb 24, 2018
 *      Author: Dr. Lee Burchett
 */

#include "SRSC.h"

void open_HITRAN_data(ulong M, HITRAN_DATA *h){
	h->n_lines = M;
	h->Sij = calloc(M,sizeof(ftype));
	h->nuij = calloc(M,sizeof(ftype));
	h->Aij = calloc(M,sizeof(ftype));
	h->gamma_air = calloc(M,sizeof(ftype));
	h->gamma_self = calloc(M,sizeof(ftype));
	h->Eij = calloc(M,sizeof(ftype));
	h->n_air = calloc(M,sizeof(ftype));
	h->delta_air = calloc(M,sizeof(ftype));
	h->g_p = calloc(M,sizeof(ftype));
	h->g_pp = calloc(M,sizeof(ftype));
	h->Q_ratio = calloc(M,sizeof(ftype));
	h->molec = calloc(M,sizeof(ulong));
	h->isotopologue = calloc(M,sizeof(ulong));
	h->mmr = calloc(M,sizeof(ulong));
}
void free_HITRAN(HITRAN_DATA *h){
	int i;
	for(i=0;i<N_HITRAN_PAR_FILES;i++){
		free(h->data[i].line);
		free_molecules(h->data+i);
	}
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
	free(h->mmr);
}
void open_weather_data(ulong P, WEATHER_DATA *w){
	w->n_points = P;
	w->P = malloc(P*sizeof(ftype));
	w->T = malloc(P*sizeof(ftype));
	w->mmr = malloc(P*sizeof(ftype));
	w->rho = malloc(P*sizeof(ftype));
	w->ev = malloc(P*sizeof(ftype));
	w->w_rho = malloc(P*sizeof(ftype));
}
void free_weather(WEATHER_DATA *w){
	w->n_points=0;
	free(w->P);
	free(w->T);
	free(w->mmr);
	free(w->rho);
	free(w->ev);
	free(w->w_rho);
}

ulong find_index_after(ftype x0, ftype *x, ulong N){
	/* Use a recursive binary search to find the index of the entry in x that
	comes before x0. Assume x is monotonic increasing. Be careful! This can
	return an index that does not exist! */
	ulong idx=0,rem=N/2;
	switch(N){
	    case 0:
	        return 0;
	        break;
	    case 1:
	        return x0 >= x[0];
	        break;
	    case 2:
	        return (x0 >= x[0]) + (x0 >= x[1]);
	        break;
	    case 3:
	        return (x0 >= x[0]) + (x0 >= x[1]) + (x0 >= x[2]);
	        break;
	    case 4:
	        return (x0 >= x[0]) + (x0 >= x[1]) + (x0 >= x[2]) + (x0 >= x[3]);
	        break;
	    default:
	        if (x0 >= x[N/2]){
	            idx = rem;
	            rem = N - idx;
	        }
	        return idx + find_index_after(x0,x + idx, rem );
	        break;
	}
}
