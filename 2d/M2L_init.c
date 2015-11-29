#include "_fmmv.h"
#include <math.h>

#ifdef DOUBLE_FOR_INIT
  #define _FLOAT1_ double
#else
  #define _FLOAT1_ _FLOAT_
#endif


typedef struct {
	int label;
	int dir;
        double re;	
	double im;
} DiagData;

enum {
	D_N_1_1, D_N_m1_1, D_N_1_m1, D_N_m1_m1,
	D_S_1_1, D_S_m1_1, D_S_1_m1, D_S_m1_m1,
	D_W_1_1, D_W_m1_1, D_W_1_m1, D_W_m1_m1,
	D_E_1_1, D_E_m1_1, D_E_1_m1, D_E_m1_m1,
	D_N_m1_5, D_N_1_5, D_S_m1_m5, D_S_1_m5,
	D_W_m5_m1, D_W_m5_1, D_E_5_m1, D_E_5_1,
	D_W_m5_m3, D_S_m3_m5, D_S_m5_m5, D_W_m5_3,
	D_N_m3_5, D_W_m5_5, D_N_3_5, D_E_5_3,
	D_N_5_5, D_S_3_m5, D_E_5_m3, D_E_5_m5,
	D_N_m2_4, D_N_0_4, D_N_2_4, D_N_4_4,
	D_E_4_2, D_E_4_0, D_E_4_m2, D_E_4_m4,
	D_S_2_m4, D_S_0_m4, D_S_m2_m4, D_S_m4_m4,
	D_W_m4_m2, D_W_m4_0, D_W_m4_2, D_W_m4_4,
};

#define N_D_X2X 52

DiagData diagData[N_D_X2X] = {
	{D_N_1_1, XN, 0.5, 0.5},
	{D_N_m1_1, XN, -0.5, 0.5},
	{D_N_1_m1, XN, 0.5, -0.5},
	{D_N_m1_m1, XN, -0.5, -0.5},
	{D_S_1_1, XS, 0.5, 0.5},
	{D_S_m1_1, XS, -0.5, 0.5},
	{D_S_1_m1, XS, 0.5, -0.5},
	{D_S_m1_m1, XS, -0.5, -0.5},
	{D_W_1_1, XW, 0.5, 0.5},
	{D_W_m1_1, XW, -0.5, 0.5},
	{D_W_1_m1, XW, 0.5, -0.5},
	{D_W_m1_m1, XW, -0.5, -0.5},
	{D_E_1_1, XE, 0.5, 0.5},
	{D_E_m1_1, XE, -0.5, 0.5},
	{D_E_1_m1, XE, 0.5, -0.5},
	{D_E_m1_m1, XE, -0.5, -0.5},
	{D_N_m1_5, XN, -0.5, 2.5},
	{D_N_1_5, XN, 0.5, 2.5},
	{D_S_m1_m5, XS, -0.5, -2.5},
	{D_S_1_m5, XS, 0.5, -2.5},
	{D_W_m5_m1, XW, -2.5, -0.5},
	{D_W_m5_1, XW, -2.5, 0.5},
	{D_E_5_m1, XE, 2.5, -0.5},
	{D_E_5_1, XE, 2.5, 0.5},
	{D_W_m5_m3, XW, -2.5, -1.5},
	{D_S_m3_m5, XS, -1.5, -2.5},
	{D_S_m5_m5, XS, -2.5, -2.5},
	{D_W_m5_3, XW, -2.5, 1.5},
	{D_N_m3_5, XN, -1.5, 2.5},
	{D_W_m5_5, XW, -2.5, 2.5},
	{D_N_3_5, XN, 1.5, 2.5},
	{D_E_5_3, XE, 2.5, 1.5},
	{D_N_5_5, XN, 2.5, 2.5},
	{D_S_3_m5, XS, 1.5, -2.5},
	{D_E_5_m3, XE, 2.5, -1.5},
	{D_E_5_m5, XE, 2.5, -2.5},
	{D_N_m2_4, XN, -1.0, 2.0},
	{D_N_0_4, XN, 0.0, 2.0},
	{D_N_2_4, XN, 1.0, 2.0},
	{D_N_4_4, XN, 2.0, 2.0},
	{D_E_4_2, XE, 2.0, 1.0},
	{D_E_4_0, XE, 2.0, 0.0},
	{D_E_4_m2, XE, 2.0, -1.0},
	{D_E_4_m4, XE, 2.0, -2.0},
	{D_S_2_m4, XS, 1.0, -2.0},
	{D_S_0_m4, XS, 0.0, -2.0},
	{D_S_m2_m4, XS, -1.0, -2.0},
	{D_S_m4_m4, XS, -2.0, -2.0},
	{D_W_m4_m2, XW, -2.0, -1.0},
	{D_W_m4_0, XW, -2.0, 0.0},
	{D_W_m4_2, XW, -2.0, 1.0},
	{D_W_m4_4, XW, -2.0, 2.0},
};

/* TODO: ws2 and ws2_reduced */
#define N_D_X2X_ws2 0
DiagData diagData_ws2[] = {};

#define N_D_X2X_ws2_reduced 0
DiagData diagData_ws2_reduced[] = {};

/* TODO: X_RRII */


/* static */ void gen_diag_X2X(_FLOAT1_ beta, _FLOAT1_ *lambda, int s_exp, int dir, _FLOAT1_ w_re, _FLOAT1_ w_im, _FLOAT1_ *d)
{
	int j;
        _FLOAT1_ exp_lw;
        if (beta==0.0) {
	    for (j=0; j<s_exp; j++) {	
		switch(dir) {
		case XE: /* k==1 */	
                        exp_lw = exp(-lambda[j]*w_re);
			d[2*j  ] = exp_lw*cos(-lambda[j]*w_im);
			d[2*j+1] = exp_lw*sin(-lambda[j]*w_im);
			break;
		case XN: /* k==2 */	
                        exp_lw = exp(-lambda[j]*w_im);
			d[2*j  ] = exp_lw*cos(+lambda[j]*w_re);
			d[2*j+1] = exp_lw*sin(+lambda[j]*w_re);
			break;
		case XW: /* k==3 */	
                        exp_lw = exp(+lambda[j]*w_re);
			d[2*j  ] = exp_lw*cos(+lambda[j]*w_im);
			d[2*j+1] = exp_lw*sin(+lambda[j]*w_im);
			break;
		case XS: /* k==4 */	
                        exp_lw = exp(+lambda[j]*w_im);
			d[2*j  ] = exp_lw*cos(-lambda[j]*w_re);
			d[2*j+1] = exp_lw*sin(-lambda[j]*w_re);
			break;
		}	
	    }	
        }
        else {
            _FLOAT1_ ll;
	    for (j=0; j<s_exp; j++) {	
                ll = sqrt(lambda[j]*lambda[j]+beta*beta);
		switch(dir) {
		case XE: /* k==1 */	
                        exp_lw = exp(-ll*w_re);
			d[2*j  ] = exp_lw*cos(-lambda[j]*w_im);
			d[2*j+1] = exp_lw*sin(-lambda[j]*w_im);
			break;
		case XN: /* k==2 */	
                        exp_lw = exp(-ll*w_im);
			d[2*j  ] = exp_lw*cos(+lambda[j]*w_re);
			d[2*j+1] = exp_lw*sin(+lambda[j]*w_re);
			break;
		case XW: /* k==3 */	
                        exp_lw = exp(+ll*w_re);
			d[2*j  ] = exp_lw*cos(+lambda[j]*w_im);
			d[2*j+1] = exp_lw*sin(+lambda[j]*w_im);
			break;
		case XS: /* k==4 */	
                        exp_lw = exp(+ll*w_im);
			d[2*j  ] = exp_lw*cos(-lambda[j]*w_re);
			d[2*j+1] = exp_lw*sin(-lambda[j]*w_re);
			break;
		}	
	    }	
        }
}	


static void init_D_X2X(FmmvHandle *FMMV, _FLOAT1_ beta)
{
	int s_exp = FMMV->s_exp;
	int k,n;
	DiagData *dd;


	if (FMMV->ws==1) {
		n = N_D_X2X;
		dd = diagData;
	}
	else { /* ws == 2 */
	   if (FMMV->reducedScheme) {
		n = N_D_X2X_ws2_reduced;
		dd = diagData_ws2_reduced;
	   }
	   else {
		n = N_D_X2X_ws2;
		dd = diagData_ws2;
	   }	
	}

	for (k=0; k<n; k++) {
  	    gen_diag_X2X(beta, FMMV->lambda, s_exp, dd[k].dir, dd[k].re, dd[k].im, 
                         FMMV->D_X2X + 2*k*s_exp);
	}	

}

static void init_M2L_Coulomb(FmmvHandle *FMMV, int level) 
{
     init_D_X2X(FMMV, 0.0);
}

static void init_M2L_Yukawa(FmmvHandle *FMMV, int level) 
{
     //init_D_X2X(FMMV, ldexp(FMMV->beta, -level));
     init_D_X2X(FMMV, ldexp(FMMV->beta/FMMV->scale, -level));  ///!!! SCALE !!!???
}    


void init_M2L(FmmvHandle *FMMV, int level) 
{
    if (level<0) {
        int n;
        int s_exp = FMMV->s_exp;

	if (FMMV->ws==1) {
		n = N_D_X2X;
	}
	else { /* ws == 2 */
	   if (FMMV->reducedScheme) {
		n = N_D_X2X_ws2_reduced;
	   }
	   else {
		n = N_D_X2X_ws2;
	   }	
	}

	FMMV->D_X2X = FMMV_MALLOC(FMMV, 2*n*s_exp*sizeof(_FLOAT1_));
        if (FMMV->beta==0) {
            init_M2L_Coulomb(FMMV, level);
        }
   }
   else {
        if (FMMV->beta!=0) {
            init_M2L_Yukawa(FMMV, level);
        }
   }
}



void finish_M2L(FmmvHandle *FMMV) 
{ 
	int s_exp = FMMV->s_exp;
	int n;

	if (FMMV->ws==1) {
		n = N_D_X2X;
	}
	else { /* ws == 2 */
	   if (FMMV->reducedScheme) {
		n = N_D_X2X_ws2_reduced;
	   }
	   else {
		n = N_D_X2X_ws2;
	   }	
	}

	FMMV_FREE(FMMV, FMMV->D_X2X, 2*n*s_exp*sizeof(_FLOAT_));


}








