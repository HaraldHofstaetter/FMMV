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
#include <math.h>
#include <stdlib.h> /* abs */

/* global variables and auxilliary routines */

/* Factorials */
//_FLOAT_ F[2*FMM_P_MAX+8];
double F[2*FMM_P_MAX+8];

/* A[I(n,m)] = (-1)**n/sqrt((n-m)!(n+m)!) */
/* n=0..FMM_P_MAX, m=0..n, A[n][-m]=A[n][+m] */
double A[FMM_P_MAX*(FMM_P_MAX+2)+1];

/* B[J(n,m)] = sqrt((n-abs(m))!/(n+abs(m))!) */
double B[(FMM_P_MAX+2)*(FMM_P_MAX+4)+1];

/* R[k] = 1/k (k>=1, dummy for k==0) */
double R[2*FMM_P_MAX+4];

_FLOAT_ _ONE_F_[1] = { 1.0 };
_FLOAT_ _ZERO_F_[1] = { 0.0 };
int _ZERO_[FMM_P_MAX+2*FMM_S_EPS_MAX+1];

_FLOAT_ zeros[10*(FMM_P_MAX+1)*(FMM_P_MAX+2)];
/* zeros hopefully large enough, should be initialized to zero */

double FF[FMM_P_MAX+1];
double FF_INV[FMM_P_MAX+1];
double YUKSCALE[FMM_P_MAX+1];
double YUKSCALE_INV[FMM_P_MAX+1];

#if (FMM_PRECISION==0)
_FLOAT_ R_single[2*FMM_P_MAX+4];
_FLOAT_ B_single[(FMM_P_MAX+2)*(FMM_P_MAX+4)+1];
_FLOAT_ FF_single[FMM_P_MAX+1];
_FLOAT_ FF_INV_single[FMM_P_MAX+1];
#endif


void double2single(int n, double *source, float *target)
{
	int i;
	for (i=0; i<n; i++) {
		target[i] = (float) source[i];
        }
}

void init_coeffs(FmmvHandle *FMMV)
{
	int n,m;
        
	
	F[0] = 1.0;
	for (n=1; n<=2*FMM_P_MAX+4; n++) {
		F[n]=n*F[n-1];
	}

	/* init A coefficients */
	for (n=0; n<=FMM_P_MAX; n++) {
		for (m=-n; m<=+n; m++) {
			//A[I(n,m)] = (n%2 ? -1.0 : +1.0)/sqrt(F[n-m]*F[n+m]);
			A[I(n,m)] = (n%2 ? -1.0 : +1.0)/sqrt(F[n-m])/sqrt(F[n+m]); // cave canem (overflow)
		}
	}
	
	/* init B coefficients */
	for (n=0; n<=FMM_P_MAX+2; n++) {
		for (m=0; m<=n; m++) {
			B[J(n,m)] = sqrt(F[n-abs(m)]/F[n+abs(m)]);
		}
	}	

	/* init R */
	R[0] = 0.0;
	for (n=1; n<2*FMM_P_MAX+4; n++) {
		R[n] = ((double) 1.0)/n;
	}


	for (n=0; n<=FMM_P_MAX+2*FMM_S_EPS_MAX; n++) {
		_ZERO_[n] = n;
	}	

        if (FMMV->beta!=0.0) {         
             /* FF[n] = n!! = 1*3*...*(2*n+1)
                FF_INV[n] = 1/FF[n]
                YUKSCALE[n] = FF[n]/beta^(n+1)
                YUKSCALE_INV[n] = 1/YUKSCALE[n]
             */   
            _FLOAT_ h = 1.0/FMMV->beta;                          
            YUKSCALE[0] = h;
            YUKSCALE_INV[0] = FMMV->beta;
            FF[0] = 1.0;
            FF_INV[0] = 1.0;
            for (n=1;n<=FMM_P_MAX;n++) {
                YUKSCALE[n] = YUKSCALE[n-1]*((_FLOAT_) (2*n+1))*h;
                YUKSCALE_INV[n] = 1.0/YUKSCALE[n];
                FF[n] = FF[n-1]*((_FLOAT_) (2*n+1));
                FF_INV[n] = 1.0/FF[n];
            }
        }

#if (FMM_PRECISION==0)
        double2single((FMM_P_MAX+2)*(FMM_P_MAX+4)+1, B, B_single);
        double2single(2*FMM_P_MAX+4, R, R_single);
        double2single(FMM_P_MAX+1, FF, FF_single);
        double2single(FMM_P_MAX+1, FF_INV, FF_INV_single);
#endif
}


