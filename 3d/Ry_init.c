/*
 * FMMV - the Fastest Multipole Method of Vienna
 * Copyright (c) 2006-2015 Harald Hofstaetter
 * http://www.harald-hofstaetter.at
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
#include <math.h> /* sin, cos, fabs, sqrt */

/* 
 *   Holger Dachsel: Fast and accurate determination of the Wigner
 *   rotation matrices in the fast multipole method
 */

#ifdef DOUBLE_FOR_INIT
  #define _FLOAT1_ double
#else
  #define _FLOAT1_ _FLOAT_
#endif

static _FLOAT1_ sqrt_tab[2*FMM_P_MAX*(FMM_P_MAX+1)];
static _FLOAT1_ inv_sqrt_tab[2*FMM_P_MAX*(FMM_P_MAX+1)];

static void wigner(int p, _FLOAT1_ theta, _FLOAT1_ **D)
{
	_FLOAT1_ pp[FMM_P_MAX+1][FMM_P_MAX+1];
	_FLOAT1_ g[FMM_P_MAX+1][FMM_P_MAX+1];
	_FLOAT1_ d1[FMM_P_MAX+1][2*FMM_P_MAX+1];
	_FLOAT1_ a1[FMM_P_MAX+1][2*FMM_P_MAX+1];
	_FLOAT1_ d2[FMM_P_MAX+1][2*FMM_P_MAX+1];
	_FLOAT1_ a2[FMM_P_MAX+1][2*FMM_P_MAX+1];
	_FLOAT1_ cos_theta_plus1_n[FMM_P_MAX+1];
	_FLOAT1_ sin_theta_n[FMM_P_MAX+1];
	_FLOAT1_ t1, t2;
	
	_FLOAT1_ cos_theta = cos(theta);
	_FLOAT1_ sin_theta = sin(theta);
	_FLOAT1_ xxx_theta = sin_theta/(1+cos_theta);

	int n, m, l, k;

	cos_theta_plus1_n[0] = 1.0;
	sin_theta_n[0] = 1.0;
	for (n=1; n<=p; n++) {
    		cos_theta_plus1_n[n] = cos_theta_plus1_n[n-1]*(1+cos_theta);
    		sin_theta_n[n] = sin_theta_n[n-1]*sin_theta;
	}	

	pp[0][0] = 1.0;
	for (l=1; l<=p; l++) {
    		pp[l][l] = sqrt_tab[2*l-1]*inv_sqrt_tab[2*l]*sin_theta*pp[l-1][l-1];
	    	pp[l][l-1] = sqrt_tab[2*l-1]*cos_theta*pp[l-1][l-1];
    		if (l>1) {
			for (k=0; k<=l-2; k++) {
            			pp[l][k] = ((2*l-1)*cos_theta*pp[l-1][k]-sqrt_tab[(l-k-1)*(l+k-1)]*pp[l-2][k])*inv_sqrt_tab[(l-k)*(l+k)];
        		}    
    		}    
	}

	g[0][0] = 1.0;
	for (l=1; l<=p; l++) {
    		g[l][0] = sqrt_tab[2*l-1]*inv_sqrt_tab[2*l]*g[l-1][0];
		for (m=1; m<=l; m++) {
        		g[l][m] = sqrt_tab[l-m+1]*inv_sqrt_tab[l+m]*g[l][m-1];
		}	
	}

	for (l=1; l<=p; l++) {
    		t1 = pp[l][0];
    		d1[0][l] = t1;
    		a1[0][l] = fabs(t1);
		
		for (k=1; k<=l; k++) {
        		t1 = pp[l][k];
        		d1[0][l+k] = (k%2?-1:+1) * t1;
        		d1[0][l-k] = t1;
        		a1[0][l+k] = fabs(t1);
        		a1[0][l-k] = fabs(t1);
		}
		
		for (m=0; m<=l-1; m++) {
			for (k=-l+1; k<=l; k++) {
            			t1 = sqrt_tab[l*(l+1)-k*(k-1)]*inv_sqrt_tab[l*(l+1)-m*(m+1)];
            			t2 = (m+k)*inv_sqrt_tab[l*(l+1)-m*(m+1)]*xxx_theta;
            			d1[m+1][l+k] = t1*d1[m][l+k-1] - t2*d1[m][l+k];
            			a1[m+1][l+k] = fabs(t1)*a1[m][l+k-1] + fabs(t2)*a1[m][l+k];
        		}
        		t1 = (l-m)*inv_sqrt_tab[l*(l+1)-m*(m+1)]*xxx_theta;
        		d1[m+1][l-l] = t1*d1[m][l-l];
        		a1[m+1][l-l] = fabs(t1)*a1[m][l-l];
    		}
    
		for (m=0; m<=l; m++) {
         		t1 = ((l+m)%2?-1:+1)*g[l][m]*cos_theta_plus1_n[m]*sin_theta_n[l-m];
         		d2[m][l+l] = t1;
         		a2[m][l+l] = fabs(t1);
		}	
		
		for (k=l; k>=-l+1; k--) {
         		t1 = (l+k)*inv_sqrt_tab[l*(l+1)-k*(k-1)]*xxx_theta;
         		d2[l][l+k-1] = t1*d2[l][l+k];
         		a2[l][l+k-1] = fabs(t1)*a2[l][l+k];
		}	
	
		for (m=l-1; m>=0; m--) {
			for (k=l; k>=-l+1; k--) {
             			t1 = sqrt_tab[l*(l+1)-m*(m+1)]*inv_sqrt_tab[l*(l+1)-k*(k-1)];
             			t2 = (m+k)*inv_sqrt_tab[l*(l+1)-k*(k-1)]*xxx_theta;
             			d2[m][l+k-1] = t1*d2[m+1][l+k] + t2*d2[m][l+k];
             			a2[m][l+k-1] = fabs(t1)*a2[m+1][l+k] + fabs(t2)*a2[m][l+k];
			}	
		}	
     
		for (m=0; m<=l; m++) {
			for (k=-l; k<=l; k++) {
             			if (a1[m][l+k] < a2[m][l+k]) {
                 			D[l][m*(2*l+1)+l+k] = d1[m][l+k];
				}	
             			else {
                 			D[l][m*(2*l+1)+l+k] = d2[m][l+k];
             			}    
         		}
     		}
	}
}	
	


void init_Ry(FmmvHandle *FMMV)
{
	int p=(FMMV->pM>FMMV->pL ? FMMV->pM : FMMV->pL);
	int coeffs_size;
	
	_FLOAT1_ **Ry_theta;
	_FLOAT1_ **Ry_minus_theta;
	_FLOAT1_ **Ry_minus_pi_minus_theta;
	_FLOAT1_ **Ry_pi_minus_theta;
	_FLOAT1_ **Ry_pi2;
	_FLOAT1_ **Ry_minus_pi2;
	_FLOAT1_ *Ry_theta_coeffs;
	_FLOAT1_ *Ry_minus_theta_coeffs;
	_FLOAT1_ *Ry_minus_pi_minus_theta_coeffs;
	_FLOAT1_ *Ry_pi_minus_theta_coeffs;
	_FLOAT1_ *Ry_pi2_coeffs;
	_FLOAT1_ *Ry_minus_pi2_coeffs;
	_FLOAT1_ *dd;
	_FLOAT1_ *D[FMM_P_MAX+1];
	_FLOAT1_ d_real, d_imag, sign;
	
	int m, n, l, ii, jj;
	
	Ry_theta = FMMV_MALLOC(FMMV, 2*(p+1)*sizeof(_FLOAT1_*));
	Ry_minus_theta = FMMV_MALLOC(FMMV, 2*(p+1)*sizeof(_FLOAT1_*));
	Ry_pi_minus_theta = FMMV_MALLOC(FMMV, 2*(p+1)*sizeof(_FLOAT1_*));
	Ry_minus_pi_minus_theta = FMMV_MALLOC(FMMV, 2*(p+1)*sizeof(_FLOAT1_*));
	Ry_pi2 = FMMV_MALLOC(FMMV, 2*(p+1)*sizeof(_FLOAT1_*));
	Ry_minus_pi2 = FMMV_MALLOC(FMMV, 2*(p+1)*sizeof(_FLOAT1_*));
	
	coeffs_size = (p+1)*(2*p*p+4*p+3)/3;
	Ry_theta_coeffs = FMMV_MALLOC(FMMV, coeffs_size*sizeof(_FLOAT1_));
	Ry_minus_theta_coeffs = FMMV_MALLOC(FMMV, coeffs_size*sizeof(_FLOAT1_));
	Ry_pi_minus_theta_coeffs = FMMV_MALLOC(FMMV, coeffs_size*sizeof(_FLOAT1_));
	Ry_minus_pi_minus_theta_coeffs = FMMV_MALLOC(FMMV, coeffs_size*sizeof(_FLOAT1_));
	Ry_pi2_coeffs = FMMV_MALLOC(FMMV, coeffs_size*sizeof(_FLOAT1_));
	Ry_minus_pi2_coeffs = FMMV_MALLOC(FMMV, coeffs_size*sizeof(_FLOAT1_));
	dd = FMMV_MALLOC(FMMV, (((4*p+3)*(p+2)*(p+1))/6)*sizeof(_FLOAT1_));

	Ry_theta[0] = Ry_theta_coeffs;
	Ry_minus_theta[0] = Ry_minus_theta_coeffs;
	Ry_pi_minus_theta[0] = Ry_pi_minus_theta_coeffs;
	Ry_minus_pi_minus_theta[0] = Ry_minus_pi_minus_theta_coeffs;
	Ry_pi2[0] = Ry_pi2_coeffs;
	Ry_minus_pi2[0] = Ry_minus_pi2_coeffs;
	Ry_theta[1] = 0;
	Ry_minus_theta[1] = 0;
	Ry_pi_minus_theta[1] = 0;
	Ry_minus_pi_minus_theta[1] = 0;
	Ry_pi2[1] = 0;
	Ry_minus_pi2[1] = 0;
	ii = 1;
	for (n=1; n<=p; n++) {
		Ry_theta[2*n] = Ry_theta_coeffs + ii;
		Ry_minus_theta[2*n] = Ry_minus_theta_coeffs + ii;
		Ry_pi_minus_theta[2*n] = Ry_pi_minus_theta_coeffs + ii;
		Ry_minus_pi_minus_theta[2*n] = Ry_minus_pi_minus_theta_coeffs + ii;
		Ry_pi2[2*n] = Ry_pi2_coeffs + ii;
		Ry_minus_pi2[2*n] = Ry_minus_pi2_coeffs + ii;
		ii += (n+1)*(n+1); 
		Ry_theta[2*n+1] = Ry_theta_coeffs + ii;
		Ry_minus_theta[2*n+1] = Ry_minus_theta_coeffs + ii;
		Ry_pi_minus_theta[2*n+1] = Ry_pi_minus_theta_coeffs + ii;
		Ry_minus_pi_minus_theta[2*n+1] = Ry_minus_pi_minus_theta_coeffs + ii;
		Ry_pi2[2*n+1] = Ry_pi2_coeffs + ii;
		Ry_minus_pi2[2*n+1] = Ry_minus_pi2_coeffs + ii;
		ii += n*n;
	}	

	D[0] = dd;
	ii = 1;
	for (n=1; n<=p; n++) {
	 	D[n] = dd + ii;
		ii += (n+1)*(2*n+1);
	}	

 	sqrt_tab[0] = 0;
 	inv_sqrt_tab[0] = 0; /* dummy */
	for (n=1; n<2*p*(p+1); n++) {
	 	sqrt_tab[n] = sqrt((_FLOAT1_) n);
		inv_sqrt_tab[n] = 1.0/sqrt_tab[n];
	}	

	wigner(p, 0.95531661812450927816, D);  /* theta = arctan(sqrt(2)) */

	Ry_minus_theta[0][0] = 1.0;
	Ry_theta[0][0] = 1.0;
	Ry_minus_pi_minus_theta[0][0] = 1.0;
	Ry_pi_minus_theta[0][0] = 1.0;
	for (l=1; l<=p; l++) {
		for (m=0; m<=l; m++) {
        		d_real = D[l][m*(2*l+1)+l];
			ii = m*(l+1);
		        Ry_minus_theta[2*l][ii] = d_real;
		        Ry_theta[2*l][ii] = (m%2?-1:+1)*d_real;
		        Ry_minus_pi_minus_theta[2*l][ii] = ((l+m)%2?-1:+1)*d_real;
		        Ry_pi_minus_theta[2*l][ii] = (l%2?-1:+1)*d_real;
			for (n=1; n<=l; n++) {
				sign = (n%2?-1:+1);
        			d_real = D[l][m*(2*l+1)+l+n] + sign*D[l][m*(2*l+1)+l-n]; 
        			d_imag = D[l][m*(2*l+1)+l+n] - sign*D[l][m*(2*l+1)+l-n]; 
				ii = m*(l+1)+n; /* real */
				jj = (m-1)*l+n-1;   /* imag */
		       		Ry_minus_theta[2*l][ii] = d_real;
			      	Ry_minus_theta[2*l+1][jj] = d_imag;
				sign = ((m+n)%2?-1:+1);
			        Ry_theta[2*l][ii] = sign*d_real;
			        Ry_theta[2*l+1][jj] = sign*d_imag;
				sign = ((m+n+l)%2?-1:+1);
			        Ry_minus_pi_minus_theta[2*l][ii] = sign*d_real;
			        Ry_minus_pi_minus_theta[2*l+1][jj] = -sign*d_imag;
				sign = (l%2?-1:+1);
			        Ry_pi_minus_theta[2*l][ii] = sign*d_real;
			        Ry_pi_minus_theta[2*l+1][jj] = -sign*d_imag;
			}	
		}	
        }
		
	wigner(p, 1.5707963267948966192, D);  /* theta = pi/2 */
	Ry_minus_pi2[0][0] = 1.0;
	Ry_pi2[0][0] = 1.0;
	for (l=1; l<=p; l++) {
		for (m=0; m<=l; m++) {
			ii = m*(l+1);
			if ((l+m)%2) {
			        Ry_minus_pi2[2*l][ii] = 0.0;
		        	Ry_pi2[2*l][ii] = 0.0;
			}
			else {
	        		d_real = D[l][m*(2*l+1)+l];
			        Ry_minus_pi2[2*l][ii] = d_real;
		        	Ry_pi2[2*l][ii] = (m%2?-1:+1)*d_real;
			}	
			for (n=1; n<=l; n++) {
				ii = m*(l+1)+n; /* real */
				jj = (m-1)*l+n-1;   /* imag */
				if ((l+m+n)%2) {
					sign = (n%2?-1:+1);
        				d_imag = D[l][m*(2*l+1)+l+n] - sign*D[l][m*(2*l+1)+l-n]; 
				      	Ry_minus_pi2[2*l+1][jj] = d_imag;
					sign = ((m+n)%2?-1:+1);
			        	Ry_pi2[2*l+1][jj] = sign*d_imag;
			       		Ry_minus_pi2[2*l][ii] = 0.0;
				        Ry_pi2[2*l][ii] = 0.0;
				}	
				else {	
					sign = (n%2?-1:+1);
        				d_real = D[l][m*(2*l+1)+l+n] + sign*D[l][m*(2*l+1)+l-n]; 
			       		Ry_minus_pi2[2*l][ii] = d_real;
					sign = ((m+n)%2?-1:+1);
				        Ry_pi2[2*l][ii] = sign*d_real;
				      	Ry_minus_pi2[2*l+1][jj] = 0.0;
			        	Ry_pi2[2*l+1][jj] = 0.0;
				}	
			}	
		}	
        }
		
	FMMV_FREE(FMMV, dd, (((4*p+3)*(p+2)*(p+1))/6)*sizeof(_FLOAT1_));

	#ifdef DOUBLE_FOR_INIT
	FMMV_FREE(FMMV, Ry_theta, 2*(p+1)*sizeof(_FLOAT1_*));
	FMMV_FREE(FMMV, Ry_minus_theta, 2*(p+1)*sizeof(_FLOAT1_*));
	FMMV_FREE(FMMV, Ry_pi_minus_theta, 2*(p+1)*sizeof(_FLOAT1_*));
	FMMV_FREE(FMMV, Ry_minus_pi_minus_theta, 2*(p+1)*sizeof(_FLOAT1_*));
	FMMV_FREE(FMMV, Ry_pi2, 2*(p+1)*sizeof(_FLOAT1_*));
	FMMV_FREE(FMMV, Ry_minus_pi2, 2*(p+1)*sizeof(_FLOAT1_*));

	FMMV->Ry_theta_coeffs = FMMV_MALLOC(FMMV, coeffs_size*sizeof(_FLOAT_));
	FMMV->Ry_minus_theta_coeffs = FMMV_MALLOC(FMMV, coeffs_size*sizeof(_FLOAT_));
	FMMV->Ry_pi_minus_theta_coeffs = FMMV_MALLOC(FMMV, coeffs_size*sizeof(_FLOAT_));
	FMMV->Ry_minus_pi_minus_theta_coeffs = FMMV_MALLOC(FMMV, coeffs_size*sizeof(_FLOAT_));
	FMMV->Ry_pi2_coeffs = FMMV_MALLOC(FMMV, coeffs_size*sizeof(_FLOAT_));
	FMMV->Ry_minus_pi2_coeffs = FMMV_MALLOC(FMMV, coeffs_size*sizeof(_FLOAT_));

	double2single(coeffs_size, Ry_theta_coeffs, FMMV->Ry_theta_coeffs);
	double2single(coeffs_size, Ry_minus_theta_coeffs, FMMV->Ry_minus_theta_coeffs);
	double2single(coeffs_size, Ry_pi_minus_theta_coeffs, FMMV->Ry_pi_minus_theta_coeffs);
	double2single(coeffs_size, Ry_minus_pi_minus_theta_coeffs, FMMV->Ry_minus_pi_minus_theta_coeffs);
	double2single(coeffs_size, Ry_pi2_coeffs, FMMV->Ry_pi2_coeffs);
	double2single(coeffs_size, Ry_minus_pi2_coeffs, FMMV->Ry_minus_pi2_coeffs);

	FMMV->Ry_theta = FMMV_MALLOC(FMMV, 2*(p+1)*sizeof(_FLOAT_*));
	FMMV->Ry_minus_theta = FMMV_MALLOC(FMMV, 2*(p+1)*sizeof(_FLOAT_*));
	FMMV->Ry_pi_minus_theta = FMMV_MALLOC(FMMV, 2*(p+1)*sizeof(_FLOAT_*));
	FMMV->Ry_minus_pi_minus_theta = FMMV_MALLOC(FMMV, 2*(p+1)*sizeof(_FLOAT_*));
	FMMV->Ry_pi2 = FMMV_MALLOC(FMMV, 2*(p+1)*sizeof(_FLOAT_*));
	FMMV->Ry_minus_pi2 = FMMV_MALLOC(FMMV, 2*(p+1)*sizeof(_FLOAT_*));
	
	FMMV->Ry_theta[0] = FMMV->Ry_theta_coeffs;
	FMMV->Ry_minus_theta[0] = FMMV->Ry_minus_theta_coeffs;
	FMMV->Ry_pi_minus_theta[0] = FMMV->Ry_pi_minus_theta_coeffs;
	FMMV->Ry_minus_pi_minus_theta[0] = FMMV->Ry_minus_pi_minus_theta_coeffs;
	FMMV->Ry_pi2[0] = FMMV->Ry_pi2_coeffs;
	FMMV->Ry_minus_pi2[0] = FMMV->Ry_minus_pi2_coeffs;
	FMMV->Ry_theta[1] = 0;
	FMMV->Ry_minus_theta[1] = 0;
	FMMV->Ry_pi_minus_theta[1] = 0;
	FMMV->Ry_minus_pi_minus_theta[1] = 0;
	FMMV->Ry_pi2[1] = 0;
	FMMV->Ry_minus_pi2[1] = 0;
	ii = 1;
	for (n=1; n<=p; n++) {
		FMMV->Ry_theta[2*n] = FMMV->Ry_theta_coeffs + ii;
		FMMV->Ry_minus_theta[2*n] = FMMV->Ry_minus_theta_coeffs + ii;
		FMMV->Ry_pi_minus_theta[2*n] = FMMV->Ry_pi_minus_theta_coeffs + ii;
		FMMV->Ry_minus_pi_minus_theta[2*n] = FMMV->Ry_minus_pi_minus_theta_coeffs + ii;
		FMMV->Ry_pi2[2*n] = FMMV->Ry_pi2_coeffs + ii;
		FMMV->Ry_minus_pi2[2*n] = FMMV->Ry_minus_pi2_coeffs + ii;
		ii += (n+1)*(n+1); 
		FMMV->Ry_theta[2*n+1] = FMMV->Ry_theta_coeffs + ii;
		FMMV->Ry_minus_theta[2*n+1] = FMMV->Ry_minus_theta_coeffs + ii;
		FMMV->Ry_pi_minus_theta[2*n+1] = FMMV->Ry_pi_minus_theta_coeffs + ii;
		FMMV->Ry_minus_pi_minus_theta[2*n+1] = FMMV->Ry_minus_pi_minus_theta_coeffs + ii;
		FMMV->Ry_pi2[2*n+1] = FMMV->Ry_pi2_coeffs + ii;
		FMMV->Ry_minus_pi2[2*n+1] = FMMV->Ry_minus_pi2_coeffs + ii;
		ii += n*n;
	}	
	FMMV_FREE(FMMV, Ry_theta_coeffs, coeffs_size*sizeof(_FLOAT1_));
	FMMV_FREE(FMMV, Ry_minus_theta_coeffs, coeffs_size*sizeof(_FLOAT1_));
	FMMV_FREE(FMMV, Ry_pi_minus_theta_coeffs, coeffs_size*sizeof(_FLOAT1_));
	FMMV_FREE(FMMV, Ry_minus_pi_minus_theta_coeffs, coeffs_size*sizeof(_FLOAT1_));
	FMMV_FREE(FMMV, Ry_pi2_coeffs, coeffs_size*sizeof(_FLOAT1_));
	FMMV_FREE(FMMV, Ry_minus_pi2_coeffs, coeffs_size*sizeof(_FLOAT1_));
	#else
	FMMV->Ry_theta = Ry_theta;
	FMMV->Ry_minus_theta = Ry_minus_theta;
	FMMV->Ry_pi_minus_theta = Ry_pi_minus_theta;
	FMMV->Ry_minus_pi_minus_theta = Ry_minus_pi_minus_theta;
	FMMV->Ry_pi2 = Ry_pi2;
	FMMV->Ry_minus_pi2 = Ry_minus_pi2;

	FMMV->Ry_theta_coeffs = Ry_theta_coeffs;
	FMMV->Ry_minus_theta_coeffs = Ry_minus_theta_coeffs;
	FMMV->Ry_pi_minus_theta_coeffs = Ry_pi_minus_theta_coeffs;
	FMMV->Ry_minus_pi_minus_theta_coeffs = Ry_minus_pi_minus_theta_coeffs;
	FMMV->Ry_pi2_coeffs = Ry_pi2_coeffs;
	FMMV->Ry_minus_pi2_coeffs = Ry_minus_pi2_coeffs;
	#endif
}


void finish_Ry(FmmvHandle *FMMV)
{
	int p=(FMMV->pM>FMMV->pL ? FMMV->pM : FMMV->pL);
	int coeffs_size;
	
	FMMV_FREE(FMMV, FMMV->Ry_theta, 2*(p+1)*sizeof(_FLOAT_*));
	FMMV_FREE(FMMV, FMMV->Ry_minus_theta, 2*(p+1)*sizeof(_FLOAT_*));
	FMMV_FREE(FMMV, FMMV->Ry_pi_minus_theta, 2*(p+1)*sizeof(_FLOAT_*));
	FMMV_FREE(FMMV, FMMV->Ry_minus_pi_minus_theta, 2*(p+1)*sizeof(_FLOAT_*));
	FMMV_FREE(FMMV, FMMV->Ry_pi2, 2*(p+1)*sizeof(_FLOAT_*));
	FMMV_FREE(FMMV, FMMV->Ry_minus_pi2, 2*(p+1)*sizeof(_FLOAT_*));

	coeffs_size = (p+1)*(2*p*p+4*p+3)/3;
	FMMV_FREE(FMMV, FMMV->Ry_theta_coeffs, coeffs_size*sizeof(_FLOAT_));
	FMMV_FREE(FMMV, FMMV->Ry_minus_theta_coeffs, coeffs_size*sizeof(_FLOAT_));
	FMMV_FREE(FMMV, FMMV->Ry_pi_minus_theta_coeffs, coeffs_size*sizeof(_FLOAT_));
	FMMV_FREE(FMMV, FMMV->Ry_minus_pi_minus_theta_coeffs, coeffs_size*sizeof(_FLOAT_));
	FMMV_FREE(FMMV, FMMV->Ry_pi2_coeffs, coeffs_size*sizeof(_FLOAT_));
	FMMV_FREE(FMMV, FMMV->Ry_minus_pi2_coeffs, coeffs_size*sizeof(_FLOAT_));
}
		
