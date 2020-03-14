#include "SRSC.h"
#include "HITRAN_Parser.h"

static ctype Lorentzian(ftype Sij,
				 ftype Qr,
				 ftype Epp,
				 ftype nuij,
				 ftype n_air,
				 ftype gamma_air,
				 ftype gamma_self,
				 ftype delta_air,
				 ftype P,
				 ftype T,
				 ftype partial_p,
				 ftype nu
);
static ftype Lorentzian_Gamma(
				ftype T,
				ftype n,
				ftype g_a,
				ftype g_s,
				ftype P_a,
				ftype P_s);

ctype calc_line_at(HITRAN_DATA *H,
							WEATHER_DATA *W,
							ftype nu,
							ulong line_idx,
							ulong wx_idx,
							LINE_SHAPE_FUNCTION fcn)
{
	switch(fcn){
	case LORENTZIAN:
		if(H->molec[line_idx]==H2O)
			return Lorentzian(H->Sij[line_idx],
							  H->Q_ratio[line_idx],
							  H->Eij[line_idx],
							  H->nuij[line_idx],
							  H->n_air[line_idx],
							  H->gamma_air[line_idx],
							  H->gamma_self[line_idx],
							  H->delta_air[line_idx],
							  W->P[wx_idx],
							  W->T[wx_idx],
							  W->ev[wx_idx],
							  nu);
		else
			return Lorentzian(H->Sij[line_idx],
							  H->Q_ratio[line_idx],
							  H->Eij[line_idx],
							  H->nuij[line_idx],
							  H->n_air[line_idx],
							  H->gamma_air[line_idx],
							  H->gamma_self[line_idx],
							  H->delta_air[line_idx],
							  W->P[wx_idx],
							  W->T[wx_idx],
							  W->P[wx_idx]*H->mmr[line_idx],
							  nu);
		break;
	case GAUSSIAN:
		if(H->molec[line_idx]==H2O)
			return Lorentzian(H->Sij[line_idx],
							  H->Q_ratio[line_idx],
							  H->Eij[line_idx],
							  H->nuij[line_idx],
							  H->n_air[line_idx],
							  H->gamma_air[line_idx],
							  H->gamma_self[line_idx],
							  H->delta_air[line_idx],
							  W->P[wx_idx],
							  W->T[wx_idx],
							  W->ev[wx_idx],
							  nu);
		else
			return Lorentzian(H->Sij[line_idx],
							  H->Q_ratio[line_idx],
							  H->Eij[line_idx],
							  H->nuij[line_idx],
							  H->n_air[line_idx],
							  H->gamma_air[line_idx],
							  H->gamma_self[line_idx],
							  H->delta_air[line_idx],
							  W->P[wx_idx],
							  W->T[wx_idx],
							  W->P[wx_idx]*H->mmr[line_idx],
							  nu);
		break;
	case VOIGT:
		if(H->molec[line_idx]==H2O)
			return Lorentzian(H->Sij[line_idx],
							  H->Q_ratio[line_idx],
							  H->Eij[line_idx],
							  H->nuij[line_idx],
							  H->n_air[line_idx],
							  H->gamma_air[line_idx],
							  H->gamma_self[line_idx],
							  H->delta_air[line_idx],
							  W->P[wx_idx],
							  W->T[wx_idx],
							  W->ev[wx_idx],
							  nu);
		else
			return Lorentzian(H->Sij[line_idx],
							  H->Q_ratio[line_idx],
							  H->Eij[line_idx],
							  H->nuij[line_idx],
							  H->n_air[line_idx],
							  H->gamma_air[line_idx],
							  H->gamma_self[line_idx],
							  H->delta_air[line_idx],
							  W->P[wx_idx],
							  W->T[wx_idx],
							  W->P[wx_idx]*H->mmr[line_idx],
							  nu);
		break;
	default:
		if(H->molec[line_idx]==H2O)
			return Lorentzian(H->Sij[line_idx],
							  H->Q_ratio[line_idx],
							  H->Eij[line_idx],
							  H->nuij[line_idx],
							  H->n_air[line_idx],
							  H->gamma_air[line_idx],
							  H->gamma_self[line_idx],
							  H->delta_air[line_idx],
							  W->P[wx_idx],
							  W->T[wx_idx],
							  W->ev[wx_idx],
							  nu);
		else
			return Lorentzian(H->Sij[line_idx],
							  H->Q_ratio[line_idx],
							  H->Eij[line_idx],
							  H->nuij[line_idx],
							  H->n_air[line_idx],
							  H->gamma_air[line_idx],
							  H->gamma_self[line_idx],
							  H->delta_air[line_idx],
							  W->P[wx_idx],
							  W->T[wx_idx],
							  W->P[wx_idx]*H->mmr[line_idx],
							  nu);

		break;
	}
	return 1.1*nu + 0.0*I;
}
ctype continuum(WEATHER_DATA *w, ftype nu,ulong i){
	return 0.0 + 0.0*I;
}
ctype Lorentzian(ftype Sij,
				 ftype Qr,
				 ftype Epp,
				 ftype nuij,
				 ftype n_air,
				 ftype gamma_air,
				 ftype gamma_self,
				 ftype delta_air,
				 ftype P,
				 ftype T,
				 ftype partial_p,
				 ftype nu
){
	ftype S = Sij*Qr*
			exp(-_HITRAN_c2*Epp/T)*
			(1.0-exp(-_HITRAN_c2*nuij/T))/
			exp(-_HITRAN_c2*Epp/_HITRAN_Tref)/
			(1.0-exp(-_HITRAN_c2*nuij/_HITRAN_Tref));
	/* TODO gamma is weather and line dependent, but not dependent on nu. It should be
	 * calculated elsewhere. */
	ftype gamma = Lorentzian_Gamma(T,n_air,gamma_air,gamma_self,P-partial_p,partial_p);
	/* TODO nu_0 is weather and line dependent, but not dependent on nu. It should be
	 * calculated elsewhere. */
	ftype nu_0 = nuij + delta_air*P;
	/* TODO Add the complex part to this*/
	ctype f = gamma/(M_PI*(gamma*gamma + (nu - nu_0)*(nu - nu_0))) + 0.0*I;

	return S*f;
}
ftype Lorentzian_Gamma(ftype T, ftype n, ftype g_a, ftype g_s, ftype P_a, ftype P_s){
	/* T_ref for HITRAN is 296.0 K */
	ftype out = pow(_HITRAN_Tref/T,n);
	return out*(g_a*P_a + g_s*P_s);
}
void pre_calc(HITRAN_DATA *H,
		      WEATHER_DATA *W,
			  LINE_SHAPE_PRE_CALC *out,
			  ulong line_idx,
			  ulong wx_idx,
			  LINE_SHAPE_FUNCTION fcn){
	/* TODO: Calculate this here;*/
	ftype Qr = find_Q_ratio(H->data,
			                H->molec[line_idx],
							H->isotopologue[line_idx],
							W->T[wx_idx]);
	/* TODO Eventually, Add in functions for other line shapes */
	if(H->molec[line_idx]==H2O){
		out->gamma = Lorentzian_Gamma(W->T[wx_idx],H->n_air[line_idx],
			H->gamma_air[line_idx],H->gamma_self[line_idx],
			W->P[wx_idx],W->ev[wx_idx]);
		out->rhoS_pi = W->w_rho[wx_idx]/M_PI;
	}
	else{
		out->gamma = Lorentzian_Gamma(W->T[wx_idx],H->n_air[line_idx],
					H->gamma_air[line_idx],H->gamma_self[line_idx],
					W->P[wx_idx],W->P[wx_idx]*H->mmr[line_idx]);
		out->rhoS_pi = W->rho[wx_idx]*H->mmr[line_idx]/M_PI;
	}
	out->nu0 = H->nuij[line_idx] + H->delta_air[line_idx]*W->P[wx_idx];
	out->S = H->Sij[line_idx]*Qr*
			exp(-_HITRAN_c2 * H->Eij[line_idx] / W->T[wx_idx] )*
			(1.0 - exp(-_HITRAN_c2 * H->nuij[line_idx]/W->T[wx_idx]))/
			exp(-_HITRAN_c2*H->Eij[line_idx]/_HITRAN_Tref)/
			(1.0-exp(-_HITRAN_c2*H->nuij[line_idx]/_HITRAN_Tref));
	out->rhoS_pi *= out->S;
	out->gamma2 = out->gamma*out->gamma;
	out->nu02 = out->nu0*out->nu0;

}
ctype calc_line_final(
		HITRAN_DATA *H,
		WEATHER_DATA *W,
		ftype nu,
		ulong j,
		ulong k,
		LINE_SHAPE_PRE_CALC *p,
		LINE_SHAPE_FUNCTION *fcn){
	/* The pre calculation part should have been run */
	ftype real = p->rhoS_pi*p->gamma/(p->gamma +
			nu*nu - 2.0*nu*p->nu0 + p->nu02);
	ftype nu0MinusNu = p->nu02 - nu*nu;
	ftype fourGammaNu = 4.0*p->gamma2*nu*nu;
	ftype term1 = real/M_PI/nu/4.0;
	ftype term2 = p->rhoS_pi*nu0MinusNu/(fourGammaNu + nu0MinusNu*nu0MinusNu);
	return 1.0e5*real + I*5.0e5*(term1*term1 + term2);
}
ftype find_Q_ratio(HITRAN_FILE *H,
				   unsigned short molec,
				   unsigned short isotope_num,
				   ftype T){
	double Q_T;
	double dq_dt;
	double dt;

	if(H->molec[molec].Q_T==NULL){
		printf("Unable to find the Q ratio for molecule #%u\n",molec);
		return 1.0;
	}
	if(isotope_num==0){
		/*printf("Improper Isotope Number, 0, given for molecule #%u.\n",molec);*/
		isotope_num=1;
	}
	isotope_num -=1;
	if(H->molec[molec].n_isotope <= (int)isotope_num){
		printf("Unable to find isotope #%u for molecule #%u\n",isotope_num, molec);
		return 1.0;
	}
	if(H->molec[molec].Q_T[isotope_num].Q_T==NULL){
		printf("Unable to find the Q ratio for molecule #%u, isotope #%u\n",molec, isotope_num);
		return 1.0;
	}
	ulong idx = find_index_after(T,
			H->molec[molec].Q_T[isotope_num].T,
			H->molec[molec].Q_T[isotope_num].n);
	if(idx >= H->molec[molec].Q_T[isotope_num].n){
		 Q_T = H->molec[molec].Q_T[isotope_num].Q_T[
				 H->molec[molec].Q_T[isotope_num].n-1];
	}
	else if(idx==0){
		 Q_T = H->molec[molec].Q_T[isotope_num].Q_T[0];
	}
	else{
		/* Linear interpolation */
		 Q_T = H->molec[molec].Q_T[isotope_num].Q_T[idx-1];
		 dq_dt = (H->molec[molec].Q_T[isotope_num].Q_T[idx] -
				 H->molec[molec].Q_T[isotope_num].Q_T[idx-1])/
				(H->molec[molec].Q_T[isotope_num].T[idx] -
				 H->molec[molec].Q_T[isotope_num].T[idx-1]);
		 dt = T - H->molec[molec].Q_T[isotope_num].T[idx-1];
		 Q_T += dt*dq_dt;
	}
	return (ftype) H->molec[molec].Q_296K[isotope_num] / Q_T;

}
