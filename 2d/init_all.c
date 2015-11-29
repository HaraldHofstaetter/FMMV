#include"_fmmv.h"
#include <math.h>


int init_all(FmmvHandle *FMMV)
{
	init_quad_coeffs(FMMV);
	init_coeffs(FMMV);

	return 0;
}


int finish_all(FmmvHandle *FMMV)
{
	finish_quad_coeffs(FMMV);
	finish_coeffs(FMMV);

	return 0;
}


void init_coeffs(FmmvHandle *FMMV)
{
	_FLOAT_ r, f, sqrt_w;
	int n, m, p;
	int q = FMMV->s_eps;
	_FLOAT_ scale = FMMV->scale;
	_FLOAT_ sqrt_scale = sqrt(scale);

	if (FMMV->pL+1>FMMV->pM) {
		p = FMMV->pL+1;
	}
	else {
		p = FMMV->pM;
	}

	FMMV->vander_over_fact[0] = FMMV_MALLOC(FMMV, q*(p+2)*sizeof(_FLOAT_));
	for (n=0; n<q; n++) {
		sqrt_w = sqrt(FMMV->w[n]);
		FMMV->vander_over_fact[n] = FMMV->vander_over_fact[0] + (p+2)*n;
		FMMV->vander_over_fact[n][0] = 1.0/(FMMV->lambda[n]*scale); 
		FMMV->vander_over_fact[n][1] = 1.0; 
		for (m=1; m<=p; m++) {
			FMMV->vander_over_fact[n][m+1] = FMMV->vander_over_fact[n][m]*(FMMV->lambda[n]*scale)/((_FLOAT_) m);
		}
		f = sqrt_scale*sqrt_w;
		for (m=0; m<=p+1; m++) {
			FMMV->vander_over_fact[n][m] *= f;
		}	
	}	
        FMMV->C -= log(scale);
        

        FMMV->recip = FMMV_MALLOC(FMMV, (p+1)*sizeof(_FLOAT_));
        FMMV->recip[0] = 1.0;
	for (n=1;n<=p;n++) {
	     FMMV->recip[n] = 1.0/(_FLOAT_) n;
	}


        if (FMMV->beta!=0.0) {         
            double *FF = FMMV_MALLOC(FMMV, (p+1)*sizeof(double));
            double *FF_INV = FMMV_MALLOC(FMMV, (p+1)*sizeof(double));
             /* FF[n] = n!! = 2*4*...*(2*n) = 2^n*n!
                FF_INV[n] = 1/FF[n]
             */   
            FF[0] = 1.0;
            FF_INV[0] = 1.0;
            for (n=1;n<=p;n++) {
                FF[n] = FF[n-1]*((_FLOAT_) (2*n));
                FF_INV[n] = 1.0/FF[n];
            }
	    FMMV->FF = FF;
	    FMMV->FF_INV = FF_INV;

	    /* log(2) - gamma = 0.11593151565841244881 */
	    //FMMV->k0_correction = log(FMMV->beta) - 0.11593151565841244881;
	    FMMV->k0_correction = 0.0;
        }

}	

void finish_coeffs(FmmvHandle *FMMV)
{
	int p;
	int q = FMMV->s_eps;

	if (FMMV->pL+1>FMMV->pM) {
		p = FMMV->pL+1;
	}
	else {
		p = FMMV->pM;
	}
	
	FMMV_FREE(FMMV, FMMV->recip, (p+1)*sizeof(_FLOAT_));
	FMMV_FREE(FMMV, FMMV->vander_over_fact[0], q*(p+2)*sizeof(_FLOAT_));
        if (FMMV->beta!=0.0) {         
	     FMMV_FREE(FMMV, FMMV->FF, (p+1)*sizeof(double));
	     FMMV_FREE(FMMV, FMMV->FF_INV, (p+1)*sizeof(double));
	}
}	
			
