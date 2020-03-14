/*
 * FileIO.c
 *
 *  Created on: Feb 23, 2018
 *      Author: Dr. Lee Burchett
 */

#include "SRSC.h"


/*
static char *input_files[N_HITRAN_PAR_FILES] = {
		"05_HITEMP2010new.par",
		"08_HITEMP2010.par",
		"13_HITEMP2010.par",
		"CO2/02_02000-02125_HITEMP2010.par",
		"CO2/02_2125-2250_HITEMP2010.par",
		"H2O/01_2000-2250_HITEMP2010.par"

};
static char *hitran_source_folder="/mnt/lee/DataStore/HITRAN/HITEMP/";*/
static char *input_files[N_HITRAN_PAR_FILES] = {
		"ATMOS_plus_NUDET_MM_Wave.par"
/*		"test.par"*/
};
static HITRAN_FILE data[N_HITRAN_PAR_FILES];
static char *hitran_source_folder="/mnt/lee/DataStore/HITRAN/";
static ftype mixingRatio[N_HITRAN_MOLECULES] = {
		 0.0, 		/* Unknown (0)*/
		 0.1, 		/* H2O (1)*/
		 400.0e-6, 	/* CO2 (2)*/
		 1.0e-6,	/* O3  (3)*/
	     320.0e-9,  /* N2O (4)*/
		 0.15e-6,	/* CO  (5)*/
		 1.794e-6,	/* CH4 (6)*/
		 0.21, 		/* O2  (7)*/
		 2.99e-10,  /* NO  (8)*/
		 2.93e-10,  /* SO2 (9)*/
		 328.0e-9,  /* NO2 (10)*/
		 5.03e-11,  /* NH3 (11)*/
		 5.3e-11,	/* HNO3(12)*/
	     4.3987e-14,/* OH  (13)*/
		 0.0,/*	    HF, (14)*/
		 0.0,/*	   HCl, (15)*/
		 0.0,/*HBr, (16)*/
		 0.0,/*	    HI, (17)*/
		 0.0,/*	   ClO, (18)*/
		 0.0,/*	   OCS, (19)*/
		 0.0,/*	  H2CO, (20)*/
		 0.0,/*	  HOCl, (21)*/
		 0.0,/*	    N2, (22)*/
		 1.0e-6,/*	   HCN, (23)*/
		 0.0,/*	 CH3Cl, (24)*/
		 0.0,/*	  H2O2, (25)*/
		 0.0,/*	  C2H2, (26)*/
		 0.0,/*	  C2H6, (27)*/
		 0.0,/*	   PH3, (28)*/
		 0.0,/*	  COF2, (29)*/
		 0.0,/*	   SF6, (30)*/
		 0.0,/*	   H2S, (31)*/
		 0.0,/*	 HCOOH, (32)*/
		 0.0,/*	   HO2, (33)*/
		 0.0,/*	     O, (34)*/
		 0.0,/*	ClONO2, (35)*/
		 0.0,/*	   NOPlus, (36)*/
		 0.0,/*	  HOBr, (37)*/
		 0.0,/*	  C2H4, (38)*/
		 0.0,/*	 CH3OH, (39)*/
		 0.0,/*	 CH3Br, (40)*/
		 1.0e-6,/*	 CH3CN, (41)*/
		 0.0,/*	   CF4, (42)*/
		 0.0,/*	  C4H2, (43)*/
		 0.0,/*	  HC3N, (44)*/
		 0.0,/*	    H2, (45)*/
		 0.0,/*	    CS, (46)*/
		 0.0,/*	   SO3, (47)*/
};

static void sort_and_copy(HITRAN_FILE data[N_HITRAN_PAR_FILES],HITRAN_DATA *H);
static int all_data_are_read(ulong index[N_HITRAN_PAR_FILES],HITRAN_FILE data[N_HITRAN_PAR_FILES]);
static int find_next(ulong index[N_HITRAN_PAR_FILES],HITRAN_FILE data[N_HITRAN_PAR_FILES]);
static int check_loaded_data(HITRAN_DATA *H);
static void copy_line(HITRAN_DATA *H, ulong idx, HITRAN_LINE *line);

HITRAN_PARSE_ERROR LoadHITRANparData(HITRAN_DATA *H){
	/* Load up and sort all the data from several files. */
	int fnum;
	ulong M=0;
	HITRAN_PARSE_ERROR rv;
	char function_name[] = "FileIO.c:LoadHITRANparData()";
	check(read_molecular_parameters(data,hitran_source_folder));

	H->data = data;

	for(fnum=0;fnum<N_HITRAN_PAR_FILES;fnum++){
		check(copy_molecular_parameters(data,data+fnum));
		printf("Importing %s. ",input_files[fnum]);
		if (import_hitran_data(input_files[fnum], data + fnum,hitran_source_folder)){
			printf("Error reading %s%s.\n",hitran_source_folder,input_files[fnum]);
			return FILE_IO_ERROR;
		}
		printf("%lu lines loaded.\n",data[fnum].n_line);
		M += data[fnum].n_line;
	}

	open_HITRAN_data(M,H);

	sort_and_copy(data,H);

	if(check_loaded_data(H)==1)
		printf("Checked %lu points. Data loaded just fine.\n",H->n_lines);

	printf("Finished Loading and Sorting data.\n");
	/* TODO Add this functionality save_as_HITRAN_binary(H);*/

	return SUCCESS;
}

int check_loaded_data(HITRAN_DATA *H){
	ulong i;
	ftype nu=H->nuij[0];

	for(i=0;i<H->n_lines;i++){
		if(H->nuij[i]<nu){
			printf("ERROR found in HITRAN data! Nu decreased at index %lu.\n",i);
			return -1;
		}
		if(H->nuij[i]==0.0){
			printf("ERROR found in HITRAN data! 0.0 Value for Nu at %lu.\n",i);
			return -1;
		}
		if(H->Sij[i]==0.0){
			printf("ERROR found in HITRAN data! 0.0 Value for S at %lu.\n",i);
			return -1;
		}
		nu = H->nuij[i];
	}
	return 1;
}

void sort_and_copy(HITRAN_FILE data[N_HITRAN_PAR_FILES],HITRAN_DATA *H){
	ulong index[N_HITRAN_PAR_FILES] = {0},i;
	int next;

	for(i=0;i<H->n_lines;i++){
		if((next = find_next(index,data))<0) break;
		if(data[next].n_line <= index[next])continue;
		copy_line(H,i,data[next].line + index[next]);
		index[next]++;
		if(all_data_are_read(index,data))break;
	}
}

void copy_line(HITRAN_DATA *H, ulong idx, HITRAN_LINE *line){
	if(H==NULL)return;
	if(idx>H->n_lines)return;
	if(line==NULL)return;
	H->Sij[idx] = (ftype)line->S;
	H->nuij[idx] = (ftype)line->nu;
	H->Aij[idx] = (ftype)line->A;
	H->gamma_air[idx] = (ftype)line->gamma_air;
	H->gamma_self[idx] = (ftype)line->gamma_self;
	H->Eij[idx] = (ftype)line->E_prime_prime;
	H->n_air[idx] = (ftype)line->n_air;
	H->delta_air[idx] = (ftype)line->delta_air;
	H->g_p[idx] = (ftype)line->g_prime;
	H->g_pp[idx] = (ftype)line->g_prime_prime;
	H->molec[idx] = (unsigned short)line->molecule_id;
	H->isotopologue[idx] = (unsigned short)line->isotopologue_id;
	H->mmr[idx] = mixingRatio[H->molec[idx]];
}

int find_next(ulong index[N_HITRAN_PAR_FILES],HITRAN_FILE data[N_HITRAN_PAR_FILES]){
	int next=-1;
	ftype min=_ftype_max ,nu_i;
	int i;

	/* Sort of like mergesort. We want to know which input file has the
	 * lowest value of nu. So we have to go over all files.*/
	for(i=0;i<N_HITRAN_PAR_FILES;i++){
		if(index[i] < data[i].n_line){ /* Make sure we are safe to access this index */
			nu_i = data[i].line[index[i]].nu;
			if(nu_i<min){
				min=nu_i;
				next=i;
			}
		}
	}

	return next;
}

int all_data_are_read(ulong index[N_HITRAN_PAR_FILES],HITRAN_FILE data[N_HITRAN_PAR_FILES]){
	int i,is_read=0;
	for(i=0;i<N_HITRAN_PAR_FILES;i++){
		is_read += (index[i]== (data[i].n_line));
	}
	return (is_read >= N_HITRAN_PAR_FILES);
}
