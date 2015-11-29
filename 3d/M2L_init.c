/*
 * FMMV - the Fastest Multipole Method of Vienna
 * Copyright (c) 2006 Harald Hofstaetter
 * Institute for Analysis and Scientific Computing
 * Vienna University of Technology
 * 
 * This file is part of FMMV.
 * 
 * FMMV is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * FMMV is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with FMMV; if not, write to the Free Software  Foundation, Inc.,
 * 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 * 
 */

#include "_fmmv.h"
#include<math.h>
#include<stdio.h>

//#ifdef DOUBLE_FOR_INIT
//  #define _FLOAT1_ double
//#else
  #define _FLOAT1_ _FLOAT_
//#endif

extern double A[];
extern double B[];
extern double R[];


void legendreP(int p, double *P, double x)
{

	double	pmm;
	double	pm1;
	double	pm2;
	double	pml;
	double	u;

	int m, l, k, kk;
        
        if (fabs(x)<=1.0) {
	    u = sqrt(1.0 - x*x);
        }
        else {
	     u = sqrt(x*x-1.0);
        }

	pmm = 1.0;
	
	/* m==0: ***************/
	P[0] = pmm;
	pm2 = pmm;
	pml = pmm*x;
	P[1] = pml;
	k = 1;
	for (l=2; l<=p; l++) {
		pm1 = pml;
		pml = R[l]*((2*l-1)*x*pm1 - (l-1)*pm2);
		pm2 = pm1;
		k += l;
		P[k] = pml;
	}

	/* 1 <=m <= p-1: ***************/
	kk = 0;
	for (m=1; m<p; m++) {
		pmm *= (1-2*m)*u;
		kk += m;
		k = kk + m;
		P[k] = pmm;
		pm2 = pmm;
		pml = (2*m+1)*pmm*x;
		k += m+1;
		P[k] = pml;
		for (l=m+2; l<=p; l++) {
			pm1 = pml;
			pml = R[l-m] * ((2*l-1)*x*pm1 - (l+m-1)*pm2);
			pm2 = pm1;
			k += l;
			P[k] = pml;
		}
	}

	/* m==p: ***************/
        pmm *= (1-2*p)*u;
	kk += p;
	k = kk + p;
        P[k] = pmm;
}





/* from M2L.h (so inclusion of M2L.h is not necessary) */
#define N_D_X2X  30 
#define N_D_X2X_ws2 122 
#define N_D_X2X_ws2_reduced 37 

#define PI 3.1415926535897932385
/* static */ void gen_diag_X2X(_FLOAT1_ beta, _FLOAT1_ *lambda, int *M, int s_eps, int s_exp, _FLOAT1_ x, _FLOAT1_ y, _FLOAT1_ z, _FLOAT1_ *d)
{
        int k,j;
        _FLOAT1_ exp_mlkz;
        _FLOAT1_ alpha, arg;
	int s_exp2 = s_exp/2;
        int jj = 0;
        if (beta==0.0) {
	    if ((x==0.0)&&(y==0.0)) {
        	for (k=0; k<s_eps; k++) {
                	exp_mlkz = exp(-z*lambda[k]);
                	for (j=0; j<M[k]/2; j++) {
                        	#ifdef X_RRII
                        	d[jj] = (_FLOAT1_) exp_mlkz;
                        	d[s_exp2+jj] = 0.0; //(_FLOAT1_) exp_mlkz;
                        	#else
                        	d[2*jj] = (_FLOAT1_) exp_mlkz;
                        	d[2*jj+1] = 0.0; //(_FLOAT1_) exp_mlkz;
                        	#endif
                        	jj++;
                	}
        	}
	    }
	    else { 
        	for (k=0; k<s_eps; k++) {
                	exp_mlkz = exp(-z*lambda[k]);
                	for (j=0; j<M[k]/2; j++) {
                        	alpha = (2.0*PI*j)/M[k];
                        	arg = -lambda[k]*(x*cos(alpha) + y*sin(alpha));
                        	#ifdef X_RRII
                        	d[jj] = (_FLOAT1_) exp_mlkz*cos(arg);
                        	d[s_exp2+jj] = (_FLOAT1_) exp_mlkz*sin(arg);
                        	#else
                        	d[2*jj] = (_FLOAT1_) exp_mlkz*cos(arg);
                        	d[2*jj+1] = (_FLOAT1_) exp_mlkz*sin(arg);
                        	#endif
                        	jj++;
                	}
        	}
	    }
  	}
	else {
	    if ((x==0.0)&&(y==0.0)) {
        	for (k=0; k<s_eps; k++) {
                	exp_mlkz = exp(-z*(lambda[k]+beta));
                	for (j=0; j<M[k]/2; j++) {
                        	#ifdef X_RRII
                        	d[jj] = (_FLOAT1_) exp_mlkz;
                        	d[s_exp2+jj] = 0.0; //(_FLOAT1_) exp_mlkz;
                        	#else
                        	d[2*jj] = (_FLOAT1_) exp_mlkz;
                        	d[2*jj+1] = 0.0; //(_FLOAT1_) exp_mlkz;
                        	#endif
                        	jj++;
                	}
        	}
	    }
	    else { 
        	for (k=0; k<s_eps; k++) {
                	exp_mlkz = exp(-z*(lambda[k]+beta));
                	for (j=0; j<M[k]/2; j++) {
                        	alpha = (2.0*PI*j)/M[k];
                        	arg = -sqrt(lambda[k]*(lambda[k]+2.0*beta))*(x*cos(alpha) + y*sin(alpha));
                        	#ifdef X_RRII
                        	d[jj] = (_FLOAT1_) exp_mlkz*cos(arg);
                        	d[s_exp2+jj] = (_FLOAT1_) exp_mlkz*sin(arg);
                        	#else
                        	d[2*jj] = (_FLOAT1_) exp_mlkz*cos(arg);
                        	d[2*jj+1] = (_FLOAT1_) exp_mlkz*sin(arg);
                        	#endif
                        	jj++;
                	}
        	}
	    }
  	}
}


static void init_D_X2X(FmmvHandle *FMMV, _FLOAT1_ beta)
{
	int s_exp = FMMV->s_exp;
	int s_eps = FMMV->s_eps;
	_FLOAT1_ *D_X2X = FMMV->D_X2X;
	_FLOAT1_ *lambda = (_FLOAT1_*) FMMV->lambda;
	int *M = FMMV->M;
	int n;
	_FLOAT1_ x,y,z;

	n = 0;
	y = 0.5;
	for (x=-0.5; x<=+0.5; x++) {
	    for (z=-0.5; z<=+0.5; z++) {
	    	gen_diag_X2X(beta, lambda, M, s_eps, s_exp, x, y, z, D_X2X + n*s_exp); 
		n++;
	    }
	}    
	if (FMMV->ws==1) {
		z = 1.5;
		for (x=-1.5; x<=+1.5; x++) {
		    for (y=0.5; y<=+1.5; y++) {
		    	gen_diag_X2X(beta, lambda, M, s_eps, s_exp, x, y, z, D_X2X + n*s_exp); 
			n++;
		    }
		}    
		z = 2.5;
		for (x=-2.5; x<=+2.5; x++) {
		    for (y=0.5; y<=+2.5; y++) {
		    	gen_diag_X2X(beta, lambda, M, s_eps, s_exp, x, y, z, D_X2X + n*s_exp); 
			n++;
		    }
		}    
	}
	else { /* ws == 2 */
		z = 2.5;
		for (x=-2.5; x<=+2.5; x++) {
		    for (y=0.5; y<=+2.5; y++) {
		    	gen_diag_X2X(beta, lambda, M, s_eps, s_exp, x, y, z, D_X2X + n*s_exp); 
			n++;
		    }
		}    
	        if (FMMV->reducedScheme) {
		    z = 4;
		    for (x=-4; x<=+4; x+=2) {
		        for (y=0; y<=+4; y+=2) {
		        	gen_diag_X2X(beta, lambda, M, s_eps, s_exp, x, y, z, D_X2X + n*s_exp); 
		    	n++;
		        }
		    }    
	       }
	       else {
		    for (x=-4.5; x<=+4.5; x++) {
		      for (y=0.5; y<=+4.5; y++) {
		        for (z=3.5; z<=+4.5; z++) {
		    	    gen_diag_X2X(beta, lambda, M, s_eps, s_exp, x, y, z, D_X2X + n*s_exp); 
			    n++;
		        }
		      }  
		    }    

	       }
	}
}

static void init_M2L_Coulomb(FmmvHandle *FMMV, int level) 
{
        // TODO: DOUBLE_FOR_INIT !!!!
	int pM = FMMV->pM;
	int pL = FMMV->pL;
	int s_eps = FMMV->s_eps;
        int p = (pM>pL?pM:pL);
        int k,m,n, ii;
        _FLOAT1_ V[FMM_S_EPS_MAX][FMM_P_MAX+1];
        _FLOAT1_ h;
        
        init_D_X2X(FMMV, 0.0);

        //generate Vandermonde Matrix
        for (k=0;k<s_eps;k++) {
            V[k][0] = 1.0;
            V[k][1] = FMMV->lambda[k]*FMMV->scale;
            for (n=2;n<=p;n++) {
                V[k][n] = V[k][n-1]*V[k][1];
            }
            /* factor 2 to compensate missing factor
             * in FFT_X2L_finish
             * sqrt, so for M2X resp. X2L these are tranposes of each other...
             */
             h = sqrt(2.0*FMMV->w[k]*FMMV->scale/((_FLOAT1_) FMMV->M[k])); 
             for (n=0;n<=p;n++) {
                 V[k][n] *= h;
             }
        }

        ii = 0;
        for (m=0; m<=pM; m++) {
            for (n=m; n<=pM; n++) {
                for (k=0; k<s_eps; k++) { 
                    FMMV->CMX_coeffs[ii] =  fabs(A[I(n,m)]) * V[k][n];
                    ii++;
                }
            }    
        }

        if (pM!=pL) {
            ii = 0;
            for (m=0; m<=pL; m++) {
                for (n=m; n<=pL; n++) {
                    for (k=0; k<s_eps; k++) { 
                        FMMV->CXL_coeffs[ii] =  fabs(A[I(n,m)]) * V[k][n];
                        ii++;
                    }
                }    
            }
        }

        if ((FMMV->ws==2)&&(FMMV->reducedScheme)) {
            ii = 0;
            for (m=0; m<=pL; m++) {
                for (n=m; n<=pL; n++) {
                    for (k=0; k<s_eps; k++) { 
                        FMMV->CXL_reduced_coeffs[ii] =  ldexp(fabs(A[I(n,m)]) * V[k][n], n); // * 2^n
                        ii++;
                    }
                }    
            }
        }

}        


extern double YUKSCALE_INV[];

static void init_M2L_Yukawa(FmmvHandle *FMMV, int level) 
{
        int pM = FMMV->pM;
        int pL = FMMV->pL;
	int s_eps = FMMV->s_eps;
        //_FLOAT1_ P[(FMM_P_MAX+1)*(FMM_P_MAX+2)/2];
        double P[(FMM_P_MAX+1)*(FMM_P_MAX+2)/2];
        _FLOAT1_ h;
	
        int k, m, n, ii;

        //init_D_X2X(FMMV, ldexp(FMMV->beta, -level));
        init_D_X2X(FMMV, ldexp(FMMV->beta/FMMV->scale, -level));  ///!!! SCALE !!!???

        for (k=0; k<s_eps; k++) {
            h = sqrt(2.0*FMMV->w[k]*FMMV->scale/((_FLOAT1_) FMMV->M[k]))/FMMV->beta;
            legendreP(pM, P,  1.0 + FMMV->lambda[k]*FMMV->scale/ldexp(FMMV->beta, -level));
            for (m=0; m<=pM; m++) {
                h = (m%2?-fabs(h):+fabs(h));
                for (n=m; n<=pM; n++) {
                    ii = (m*(2*pM-m+3)/2  + (n-m))*s_eps + k;
                    FMMV->CMX_coeffs[ii] = ldexp( ((_FLOAT_) (2*n+1)) * h * B[J(n,m)] * P[J(n,m)],  -level*n) * YUKSCALE_INV[n];
                }
            }    
        }
        if (pM!=pL) {
            for (k=0; k<s_eps; k++) {
                h = sqrt(2.0*FMMV->w[k]*FMMV->scale/((_FLOAT1_) FMMV->M[k]))/FMMV->beta;
                legendreP(pL, P,  1.0 + FMMV->lambda[k]*FMMV->scale/ldexp(FMMV->beta, -level));
                for (m=0; m<=pL; m++) {
                h = (m%2?-fabs(h):+fabs(h));
                    for (n=m; n<=pL; n++) {
                        ii = (m*(2*pL-m+3)/2  + (n-m))*s_eps + k;
                        FMMV->CXL_coeffs[ii] = ldexp( ((_FLOAT_) (2*n+1)) * h * B[J(n,m)] * P[J(n,m)],  -level*n)  * YUKSCALE_INV[n] ;
                    }
                }    
    
            }
       } 

       if ((FMMV->ws==2)&&(FMMV->reducedScheme)) {
            for (k=0; k<s_eps; k++) {
                h = sqrt(2.0*FMMV->w[k]*FMMV->scale/((_FLOAT1_) FMMV->M[k]))/FMMV->beta;
                legendreP(pL, P,  1.0 + 2.0*FMMV->lambda[k]*FMMV->scale/ldexp(FMMV->beta, -level)); //note factor 2 !!!
                for (m=0; m<=pL; m++) {
                h = (m%2?-fabs(h):+fabs(h));
                    for (n=m; n<=pL; n++) {
                        ii = (m*(2*pL-m+3)/2  + (n-m))*s_eps + k;
                        FMMV->CXL_reduced_coeffs[ii] = ldexp( ((_FLOAT_) (2*n+1)) * h * B[J(n,m)] * P[J(n,m)],  -level*n)  * YUKSCALE_INV[n] ;
                    }
                }    
    
            }
       }
}

void init_M2L(FmmvHandle *FMMV, int level)
{
	int pM = FMMV->pM;
	int pL = FMMV->pL;
	int s_eps = FMMV->s_eps;
        int s_exp = FMMV->s_exp;
	
        int ppM, ppL, n;

        if (level<0) {
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
            FMMV->D_X2X = FMMV_MALLOC(FMMV, n*s_exp*sizeof(_FLOAT1_));

	    ppM = (pM+1)*(pM+2)*(2*pM+3)/6;
	    FMMV->CMX_coeffs = FMMV_MALLOC(FMMV, s_eps*ppM*sizeof(_FLOAT1_));
    
            if (pM!=pL) {
	        ppL = (pL+1)*(pL+2)*(2*pL+3)/6;
	        FMMV->CXL_coeffs = FMMV_MALLOC(FMMV, s_eps*ppL*sizeof(_FLOAT1_));
            }  
            else {
	        FMMV->CXL_coeffs = FMMV->CMX_coeffs;
            }
            if ((FMMV->ws==2)&&(FMMV->reducedScheme)) {
	        ppL = (pL+1)*(pL+2)*(2*pL+3)/6;
	        FMMV->CXL_reduced_coeffs = FMMV_MALLOC(FMMV, s_eps*ppL*sizeof(_FLOAT1_));
            }
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
	int pM = FMMV->pM;
	int pL = FMMV->pL;
	int s_eps = FMMV->s_eps;
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
	FMMV_FREE(FMMV, FMMV->D_X2X, n*s_exp*sizeof(_FLOAT_));
        
	FMMV_FREE(FMMV, FMMV->CMX_coeffs, (s_eps*(pM+1)*(pM+2)*(2*pM+3)/6)*sizeof(_FLOAT_));
        if (pM!=pL) {
	    FMMV_FREE(FMMV, FMMV->CXL_coeffs, (s_eps*(pL+1)*(pL+2)*(2*pL+3)/6)*sizeof(_FLOAT_));
        }
        if ((FMMV->ws==2)&&(FMMV->reducedScheme)) {
	    FMMV_FREE(FMMV, FMMV->CXL_reduced_coeffs, (s_eps*(pL+1)*(pL+2)*(2*pL+3)/6)*sizeof(_FLOAT_));
        }
}



