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
#include<math.h>



extern double A[];

//#ifdef DOUBLE_FOR_INIT
//  #define _FLOAT1_ double
//#else
  #define _FLOAT1_ _FLOAT_
//#endif




static void init_M2M_Coulomb(FmmvHandle *FMMV, int level)
{
	 int pM = FMMV->pM;
	 _FLOAT1_ **Tz_M2M = FMMV->Tz_M2M;
	 _FLOAT1_ r_n[FMM_P_MAX+1];
         _FLOAT1_ h_n[FMM_P_MAX+1];
	 _FLOAT1_ r;
	 int ii,i,j,k;
	    
	 r = 2*0.4330127018922193233818615/FMMV->scale; /* sqrt(3)/4 */
	 // r = 2.0/FMMV->scale; /* sqrt(3)/4 */
	 r_n[0] = 1.0;
	 h_n[0] = 1.0;
	 for (i=1;i<=pM; i++) {
		r_n[i] = r_n[i-1] * r;
		h_n[i] = h_n[i-1] * 0.5;
         }	

	 for (k=0; k<=pM; k++) {
		ii=0;
		for (i=0; i<=pM-k; i++) {
			for (j=0; j<=i; j++) {
				Tz_M2M[k][ii] = A[I(i-j,0)]*A[I(k+j,k)]*r_n[i]*h_n[k+i]/(A[I(k+i,k)] * r_n[j]);
				ii++;
			}
		}
	 }
}

static void init_L2L_Coulomb(FmmvHandle *FMMV, int level)
{
 	 int pL = FMMV->pL;
	 _FLOAT1_ **Tz_L2L = FMMV->Tz_L2L;
	 _FLOAT1_ r_n[FMM_P_MAX+1];
         _FLOAT1_ h_n[FMM_P_MAX+1];
	 _FLOAT1_ r;
	 int ii,i,j,k;
	    
	 r = 2*0.4330127018922193233818615/FMMV->scale; /* sqrt(3)/4 */
	 // r = 2.0/FMMV->scale; /* sqrt(3)/4 */
	 r_n[0] = 1.0;
	 h_n[0] = 1.0;
	 for (i=1;i<=pL; i++) {
		r_n[i] = r_n[i-1] * r;
		h_n[i] = h_n[i-1] * 0.5;
         }	

	 for (k=0; k<=pL; k++) {
		ii=0;
		for (i=0; i<=pL-k; i++) {
			for (j=i; j<=pL-k; j++) {
				Tz_L2L[k][ii] = ((i+j)%2?-1:1)*A[I(j-i,0)]*A[I(k+i,k)]*r_n[j]*h_n[k+j]/(A[I(k+j,k)] * r_n[i]);
				ii++;
			}
		}
	 }
}



#if (FMM_PRECISION==0)
#define SINH(x) sinh(x)
void bessel_i_scaled_double(double x, double i[], int n, int n1)
/* computes i_n(x)/x^(n+1) */
{
     int l;
     double s, x2;

     x2 = x*x;
     i[n1-1] = 0.0;
     i[n1-2] = 1.0;

     for (l=n1-3; l>=0; l--) {
         i[l] = x2*i[l+2] + (2*l+3)*i[l+1];
     }

     s = SINH(x)/(x2*i[0]);
     for (l=0;l<n; l++) {
         i[l] *= s;
     }    
}
#endif



void bessel_i(_FLOAT_ x, _FLOAT_ i[], int n, int n1);
void bessel_i_scaled(_FLOAT_ x, _FLOAT_ i[], int n, int n1);
extern double F[];
extern double YUKSCALE[];
extern double YUKSCALE_INV[];

static void init_M2M_Yukawa(FmmvHandle *FMMV, int level)
{
      int pM = FMMV->pM;
      _FLOAT_ **Tz_M2M = FMMV->Tz_M2M;

#if (FMM_PRECISION==0)
      double II[4*FMM_P_MAX+4];
#else
      _FLOAT_ II[4*FMM_P_MAX+4];
#endif
      _FLOAT_ RR[2*FMM_P_MAX+3];
      int m,n,nn,k, ii;
      _FLOAT_ r, h, hh;

      r = ldexp(0.4330127018922193233818615/FMMV->scale, -level);  /* sqrt(3)/4 */
#if 0
      RR[0] = 1.0;
      RR[1] = 1.0/(r*FMMV->beta);
      bessel_i(r*FMMV->beta, II, 2*pM+2, 4*pM+4);
      for (k=2;k<=pM;k++) {
          RR[k] = RR[k-1]*RR[1];
      }

      for (m=0; m<=pM; m++) {
          ii=0;
	  for (n=m; n<=pM; n++) {
	      for (nn=m; nn<=pM; nn++) {
	         h = 0.0;
	         for (k=m; k<=(n<nn?n:nn); k++) {
		     hh = ldexp( (F[2*k]*II[nn+n-k]*RR[k]) /
	                         (F[k+m]*F[k-m]*F[nn-k]*F[n-k]*F[k]),
   	                         -k + (level)*(n-nn) -nn );
		     h += hh;
	         }
	         Tz_M2M[m][ii] =  ((_FLOAT_) (2*nn+1))*sqrt(F[nn-m]
                                  * F[n+m]*F[nn+m]*F[n-m])*h 
			          * YUKSCALE_INV[nn]*YUKSCALE[n];
	 	 ii++;
		 
	      }
	  }
      }	  
#endif

#if (FMM_PRECISION==0)
      bessel_i_scaled_double(r*FMMV->beta, II, 2*pM+2, 4*pM+4);
#else
      bessel_i_scaled(r*FMMV->beta, II, 2*pM+2, 4*pM+4);
#endif
      RR[0] = 1.0;
      RR[1] = r*FMMV->beta;
      for (k=2;k<=2*pM+2;k++) {
          RR[k] = RR[k-1]*RR[1];
      }

      for (m=0; m<=pM; m++) {
          ii=0;
	  for (n=m; n<=pM; n++) {
	      for (nn=m; nn<=pM; nn++) {
	         h = 0.0;
	         for (k=m; k<=(n<nn?n:nn); k++) {
		     hh = ldexp( (F[2*k]*II[nn+n-k]*RR[nn+n-2*k+1]) /
	                         (F[k+m]*F[k-m]*F[nn-k]*F[n-k]*F[k]),
   	                         -k + (level)*(n-nn) -nn );
		     h += hh;
	         }
	         Tz_M2M[m][ii] =  ((_FLOAT_) (2*nn+1))*sqrt(F[nn-m]
                                  * F[n+m]*F[nn+m]*F[n-m])*h 
			          * YUKSCALE_INV[nn]*YUKSCALE[n];
	 	 ii++;
	      }
	  }
      }	  

}

static void init_L2L_Yukawa(FmmvHandle *FMMV, int level)
{
      int pL = FMMV->pL;
      _FLOAT_ **Tz_L2L = FMMV->Tz_L2L;

#if (FMM_PRECISION==0)
      double II[4*FMM_P_MAX+4];
#else
      _FLOAT_ II[4*FMM_P_MAX+4];
#endif
      _FLOAT_ RR[2*FMM_P_MAX+3];
      int m,n,nn,k, ii;
      _FLOAT_ r, h, hh;
      
      r = ldexp(0.4330127018922193233818615/FMMV->scale, -level);  /* sqrt(3)/4 */
#if 0
      RR[0] = 1.0;
      RR[1] = 1.0/(r*FMMV->beta);
      bessel_i(r*FMMV->beta, II, 2*pL+2, 4*pL+4);
      for (k=2;k<=pL;k++) {
          RR[k] = RR[k-1]*RR[1];
      }

      for (m=0; m<=pL; m++) {
          ii=0;
	  for (nn=m; nn<=pL; nn++) {
	      for (n=m; n<=pL; n++) {
	         h = 0.0;
	         for (k=m; k<=(n<nn?n:nn); k++) {
		     hh = ldexp( (F[2*k]*II[nn+n-k]*RR[k]) /
	                         (F[k+m]*F[k-m]*F[nn-k]*F[n-k]*F[k]),
   	                         -k + level*(n-nn) -nn );
		     h += hh;
	         }
	         Tz_L2L[m][ii] = ((n+nn)%2?-1:1)*((_FLOAT_) (2*nn+1))
                                 * sqrt(F[nn-m]*F[n+m]*F[nn+m]*F[n-m])*h
			         * YUKSCALE_INV[nn]*YUKSCALE[n];
	 	 ii++;
		 
	      }
	  }
      }	 
#endif

#if (FMM_PRECISION==0)
      bessel_i_scaled_double(r*FMMV->beta, II, 2*pL+2, 4*pL+4);
#else
      bessel_i_scaled(r*FMMV->beta, II, 2*pL+2, 4*pL+4);
#endif
      RR[0] = 1.0;
      RR[1] = r*FMMV->beta;
      for (k=2;k<=2*pL+1;k++) {
          RR[k] = RR[k-1]*RR[1];
      }

      for (m=0; m<=pL; m++) {
          ii=0;
	  for (nn=m; nn<=pL; nn++) {
	      for (n=m; n<=pL; n++) {
	         h = 0.0;
	         for (k=m; k<=(n<nn?n:nn); k++) {
		     hh = ldexp( (F[2*k]*II[nn+n-k]*RR[nn+n-2*k+1]) /
	                         (F[k+m]*F[k-m]*F[nn-k]*F[n-k]*F[k]),
   	                         -k + level*(n-nn) -nn );
		     h += hh;
	         }
	         Tz_L2L[m][ii] = ((n+nn)%2?-1:1)*((_FLOAT_) (2*nn+1))
                                 * sqrt(F[nn-m]*F[n+m]*F[nn+m]*F[n-m])*h
			         * YUKSCALE_INV[nn]*YUKSCALE[n];
	 	 ii++;
		 
	      }
	  }
      }	  



}

void init_M2M(FmmvHandle *FMMV, int level) 
{
    int pM = FMMV->pM;
    int ppM, bb, k;
   
    if (level<0) {
        if (FMMV->beta==0.0) {
           bb = 1;
        }
        else {
           bb = 2;
        }

    	ppM = (pM+1)*(pM+2)*(bb*pM+3)/6;

	FMMV->Tz_M2M = FMMV_MALLOC(FMMV, (pM+1)*sizeof(_FLOAT1_*));
	FMMV->Tz_M2M_coeffs = FMMV_MALLOC(FMMV, ppM*sizeof(_FLOAT1_));

	for (k=0; k<=pM; k++) {
   	    FMMV->Tz_M2M[k] = FMMV->Tz_M2M_coeffs + ppM - (pM-k+1)*(pM-k+2)*(bb*(pM-k)+3)/6;
	}	
        if (FMMV->beta==0) {
            init_M2M_Coulomb(FMMV, level);
        }
    }
    else {
        if (FMMV->beta!=0) {
            init_M2M_Yukawa(FMMV, level);
        }
    }
}

void init_L2L(FmmvHandle *FMMV, int level) 
{
    int pL = FMMV->pL;
    int ppL, bb, k;

    if (level<0) {
        if (FMMV->beta==0.0) {
           bb = 1;
        }
        else {
           bb = 2;
        }
    	ppL = (pL+1)*(pL+2)*(bb*pL+3)/6;

	FMMV->Tz_L2L = FMMV_MALLOC(FMMV, (pL+1)*sizeof(_FLOAT1_*));
	FMMV->Tz_L2L_coeffs = FMMV_MALLOC(FMMV, ppL*sizeof(_FLOAT1_));

	for (k=0; k<=pL; k++) {
   	    FMMV->Tz_L2L[k] = FMMV->Tz_L2L_coeffs + ppL - (pL-k+1)*(pL-k+2)*(bb*(pL-k)+3)/6;
	}	
        if (FMMV->beta==0) {
            init_L2L_Coulomb(FMMV, level);
        }
    }
    else {
        if (FMMV->beta!=0) {
            init_L2L_Yukawa(FMMV, level);
        }
    }
}


void finish_M2M(FmmvHandle *FMMV)
{
	int pM = FMMV->pM;
	int bb;
	if (FMMV->beta==0.0) {
	    bb = 1;
	}
	else {
	    bb = 2;
	}    
	FMMV_FREE(FMMV, FMMV->Tz_M2M, (pM+1)*sizeof(_FLOAT_*));
	FMMV_FREE(FMMV, FMMV->Tz_M2M_coeffs, ((pM+1)*(pM+2)*(bb*pM+3)/6)*sizeof(_FLOAT_));
}

void finish_L2L(FmmvHandle *FMMV)
{
	int pL = FMMV->pL;
	int bb;
	if (FMMV->beta==0.0) {
	    bb = 1;
	}
	else {
	    bb = 2;
	}    
	FMMV_FREE(FMMV, FMMV->Tz_L2L, (pL+1)*sizeof(_FLOAT_*));
	FMMV_FREE(FMMV, FMMV->Tz_L2L_coeffs, ((pL+1)*(pL+2)*(bb*pL+3)/6)*sizeof(_FLOAT_));
}


