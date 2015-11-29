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

#define reJ(k,m) ((2*(m)-1)*s_eps + (k))
	
#define imJ(k,m) (2*(m)*s_eps + (k))

void init_perm(FmmvHandle *FMMV)
{
	int pM =  FMMV->pM;
	int pL =  FMMV->pL;
	int s_exp =  FMMV->s_exp;
	int s_eps =  FMMV->s_eps;
	int *P_MRT;
        int *P_MRT_plus;
        int *P_MRT_minus;
        int *P_Mriri2rrii;
	int *P_LRT;
        int *P_LRT_plus;
        int *P_LRT_minus;
        int *P_Lriri2rrii;
        int *P_X_riri2rrii;
	int *P_VF;
	int *neg_F;
	int *P_FV;
	int *neg_V;
	

	int ii, ii_plus, ii_minus, n, m, jj, j, M;

	P_MRT = FMMV_MALLOC(FMMV, ((pM+1)*(pM+1))*sizeof(int));
	P_MRT_plus = FMMV_MALLOC(FMMV, ((pM+1)*(pM+1)-((pM+1)*(pM+1))/2)*sizeof(int));
	P_MRT_minus = FMMV_MALLOC(FMMV, (((pM+1)*(pM+1))/2)*sizeof(int));
	P_Mriri2rrii = FMMV_MALLOC(FMMV, ((pM+1)*(pM+1))*sizeof(int));
	P_LRT = FMMV_MALLOC(FMMV, ((pL+1)*(pL+1))*sizeof(int));
	P_LRT_plus = FMMV_MALLOC(FMMV, ((pL+1)*(pL+1)-((pL+1)*(pL+1))/2)*sizeof(int));
	P_LRT_minus = FMMV_MALLOC(FMMV, (((pL+1)*(pL+1))/2)*sizeof(int));
	P_Lriri2rrii = FMMV_MALLOC(FMMV, ((pL+1)*(pL+1))*sizeof(int));
	P_X_riri2rrii = FMMV_MALLOC(FMMV, s_exp*sizeof(int));
	P_VF = FMMV_MALLOC(FMMV, 2*FMMV->len_F*sizeof(int));
	P_FV = FMMV_MALLOC(FMMV, (2*pL+1)*s_eps*sizeof(int));
	neg_F = FMMV_MALLOC(FMMV, 2*(2*pM+1)*s_eps*sizeof(int));
	neg_V = FMMV_MALLOC(FMMV, (2*pL+1)*s_eps*sizeof(int));
	FMMV->P_MRT = P_MRT;
	FMMV->P_MRT_plus = P_MRT_plus;
	FMMV->P_MRT_minus = P_MRT_minus;
	FMMV->P_Mriri2rrii = P_Mriri2rrii;
	FMMV->P_LRT = P_LRT;
	FMMV->P_LRT_plus = P_LRT_plus;
	FMMV->P_LRT_minus = P_LRT_minus;
	FMMV->P_Lriri2rrii = P_Lriri2rrii;
	FMMV->P_X_riri2rrii = P_X_riri2rrii;
	FMMV->P_VF = P_VF;
	FMMV->P_FV = P_FV;
	FMMV->neg_F = neg_F;
	FMMV->neg_V = neg_V;

        ii = 0;
	ii_plus = 0;
	ii_minus = 0;
        for (n=0; n<=pM; n++) { 
                P_MRT[ii] = n*n;
		if (n%2==0) {
			P_MRT_plus[ii_plus] = P_MRT[ii];
			ii_plus++;
		}
		else {
			P_MRT_minus[ii_minus] = P_MRT[ii];
			ii_minus++;
		}
                ii++;
	}
        for (m=1; m<=pM; m++) { 
        	for (n=m; n<=pM; n++) { 
                        P_MRT[ii] = n*n + m;
			if (n%2==0) {
				P_MRT_plus[ii_plus] = P_MRT[ii];
				ii_plus++;
			}
			else {
				P_MRT_minus[ii_minus] = P_MRT[ii];
				ii_minus++;
			}
                        ii++;
		}
        	for (n=m; n<=pM; n++) { 
                        P_MRT[ii] = n*n + n + m;
			if (n%2==1) {
				P_MRT_plus[ii_plus] = P_MRT[ii];
				ii_plus++;
			}
			else {
				P_MRT_minus[ii_minus] = P_MRT[ii];
				ii_minus++;
			}
                        ii++;
		}
	}

        ii = 0;
	ii_plus = 0;
	ii_minus = 0;
        for (n=0; n<=pL; n++) { 
                P_LRT[ii] = n*n;
		if (n%2==0) {
			P_LRT_plus[ii_plus] = P_LRT[ii];
			ii_plus++;
		}
		else {
			P_LRT_minus[ii_minus] = P_LRT[ii];
			ii_minus++;
		}
                ii++;
	}
        for (m=1; m<=pL; m++) { 
        	for (n=m; n<=pL; n++) { 
                        P_LRT[ii] = n*n + m;
			if (n%2==0) {
				P_LRT_plus[ii_plus] = P_LRT[ii];
				ii_plus++;
			}
			else {
				P_LRT_minus[ii_minus] = P_LRT[ii];
				ii_minus++;
			}
                        ii++;
		}
        	for (n=m; n<=pL; n++) { 
                        P_LRT[ii] = n*n + n + m;
			if (n%2==1) {
				P_LRT_plus[ii_plus] = P_LRT[ii];
				ii_plus++;
			}
			else {
				P_LRT_minus[ii_minus] = P_LRT[ii];
				ii_minus++;
			}
                        ii++;
		}
	}

        ii = 0;
        for (n=0; n<=pM; n++) { 
                /* re: */
        	for (m=0; m<=n; m++) { 
                        P_Mriri2rrii[ii] = n*(n+1) + 2*m;
                        ii += 1;
		}	
                /* im: */
        	for (m=1; m<=n; m++) { 
                        P_Mriri2rrii[ii] = n*(n+1) + 2*m + 1;
                        ii += 1;
		}
	}

        ii = 0;
        for (n=0; n<=pL; n++) { 
                /* re: */
        	for (m=0; m<=n; m++) { 
                        P_Lriri2rrii[ii] = n*(n+1) + 2*m;
                        ii += 1;
		}	
                /* im: */
        	for (m=1; m<=n; m++) { 
                        P_Lriri2rrii[ii] = n*(n+1) + 2*m + 1;
                        ii += 1;
		}
	}

	ii = 0;
	for (n=0; n<s_exp/2; n++) {
        	P_X_riri2rrii[ii] = 2*n;
		ii++;
	}	
	for (n=0; n<s_exp/2; n++) {
        	P_X_riri2rrii[ii] = 2*n+1;
		ii++;
	}	

	ii = 0;
	jj = 0;
	for (n=0; n<s_eps; n++) {
		M = FMMV->M[n];
		P_VF[jj] = n;
		P_VF[jj+1] = (2*pM+1)*s_eps; 
		jj += 2;
		for (m=1; m<=pM; m++) {
			j = m;
			if (j%M<=M/2) {
				switch (m%4) {
				case 0:	
					P_VF[jj] = reJ(n,m);
					P_VF[jj+1] = imJ(n,m);
					break;
				case 1:	
					P_VF[jj] = imJ(n,m);
					P_VF[jj+1] = reJ(n,m);
					neg_F[ii] = jj+1; 
					ii++;
					break;
				case 2:	
					P_VF[jj] = reJ(n,m);
					P_VF[jj+1] = imJ(n,m);
					neg_F[ii] = jj; 
					neg_F[ii+1] = jj+1; 
					ii += 2;
					break;
				case 3:	
					P_VF[jj] = imJ(n,m);
					P_VF[jj+1] = reJ(n,m);
					neg_F[ii] = jj; 
					ii++;
					break;
				}
				jj += 2;	
			}	
		}
		for (j=pM+1; j<M*FMMV->FFT_rep[n]-pM; j++) {
			if (j%M<=M/2) {
				P_VF[jj] = (2*pM+1)*s_eps; 
				P_VF[jj+1] = (2*pM+1)*s_eps; 
				jj +=2;
			}
		}	
		for (m=-pM; m<0; m++) {
			j = m+M*FMMV->FFT_rep[n];
			if (j%M<=M/2) {
				switch (m%4) {
				case 0:
					P_VF[jj] = reJ(n,-m);
					P_VF[jj+1] = imJ(n,-m);
					neg_F[ii] = jj+1; 
					ii++;
					break;
				case -3:
					P_VF[jj] = imJ(n,-m);
					P_VF[jj+1] = reJ(n,-m);
					break;
				case -2:
					P_VF[jj] = reJ(n,-m);
					P_VF[jj+1] = imJ(n,-m);
					neg_F[ii] = jj; 
					ii++;
					break;
				case -1:
					P_VF[jj] = imJ(n,-m);
					P_VF[jj+1] = reJ(n,-m);
					neg_F[ii] = jj; 
					neg_F[ii+1] = jj+1; 
					ii += 2;
					break;
				}	
				jj += 2;
			}
		}	
	}
	FMMV->len_neg_F = ii;	


	
	ii = 0;
	jj = 0;
	for (n=0; n<s_eps; n++) {
		M = FMMV->M[n];
		P_FV[n] = 2*jj;
		for (m=1; m<=pL; m++) {
			if (m%M<=M/2) {
				j = 2*(jj+m%M);
				switch (m%4) {
				case 0:	
					P_FV[reJ(n,m)] = j;
					P_FV[imJ(n,m)] = j+1;
					break;
				case 1:	
					P_FV[imJ(n,m)] = j;
					P_FV[reJ(n,m)] = j+1;
					neg_V[ii] = imJ(n,m);
					ii++;
					break;
				case 2:
					P_FV[reJ(n,m)] = j;
					P_FV[imJ(n,m)] = j+1;
					neg_V[ii] = reJ(n,m);
					neg_V[ii+1] = imJ(n,m);
					ii += 2;
					break;
				case 3:	
					P_FV[imJ(n,m)] = j;
					P_FV[reJ(n,m)] = j+1;
					neg_V[ii] = reJ(n,m);
					ii++;
					break;
				}	
			}
			else {
		  	    j = 2*(jj+M-m%M);
			    if ((m%M)%2) { /* re -> -re */
				switch (m%4) {  
				case 0:	
					P_FV[reJ(n,m)] = j;
					P_FV[imJ(n,m)] = j+1;
					neg_V[ii] = reJ(n,m);
					ii++;
					break;
				case 1:	  
					P_FV[imJ(n,m)] = j;
					P_FV[reJ(n,m)] = j+1;
					break;
				case 2:  
					P_FV[reJ(n,m)] = j;
					P_FV[imJ(n,m)] = j+1;
					neg_V[ii] = imJ(n,m);
					ii ++;
					break;
				case 3:	 
					P_FV[imJ(n,m)] = j;
					P_FV[reJ(n,m)] = j+1;
					neg_V[ii] = imJ(n,m);
					neg_V[ii+1] = reJ(n,m);
					ii+=2;
					break;
				}	
			    }
			    else {
				switch (m%4) { /* im -> -im */ 
				case 0:	
					P_FV[reJ(n,m)] = j;
					P_FV[imJ(n,m)] = j+1;
					neg_V[ii] = imJ(n,m);
					ii++;
					break;
				case 1:	  
					P_FV[imJ(n,m)] = j;
					P_FV[reJ(n,m)] = j+1;
					neg_V[ii] = reJ(n,m);
					neg_V[ii+1] = imJ(n,m);
					ii += 2;
					break;
				case 2:  
					P_FV[reJ(n,m)] = j;
					P_FV[imJ(n,m)] = j+1;
					neg_V[ii] = reJ(n,m);
					ii ++;
					break;
				case 3:	 
					P_FV[imJ(n,m)] = j;
					P_FV[reJ(n,m)] = j+1;
					break;
				}	
			    }
			}	
		}
		jj += M/2 + 1;	
	}
	FMMV->len_neg_V = ii;	

}	

void finish_perm(FmmvHandle *FMMV)
{
	int pM =  FMMV->pM;
	int pL =  FMMV->pL;
	int s_exp =  FMMV->s_exp;
	int s_eps =  FMMV->s_eps;

	FMMV_FREE(FMMV, FMMV->P_MRT, ((pM+1)*(pM+1))*sizeof(int));
	FMMV_FREE(FMMV, FMMV->P_MRT_plus, ((pM+1)*(pM+1)-((pM+1)*(pM+1))/2)*sizeof(int));
	FMMV_FREE(FMMV, FMMV->P_MRT_minus, (((pM+1)*(pM+1))/2)*sizeof(int));
	FMMV_FREE(FMMV, FMMV->P_Mriri2rrii, ((pM+1)*(pM+1))*sizeof(int));
	FMMV_FREE(FMMV, FMMV->P_LRT, ((pL+1)*(pL+1))*sizeof(int));
	FMMV_FREE(FMMV, FMMV->P_LRT_plus, ((pL+1)*(pL+1)-((pL+1)*(pL+1))/2)*sizeof(int));
	FMMV_FREE(FMMV, FMMV->P_LRT_minus, (((pL+1)*(pL+1))/2)*sizeof(int));
	FMMV_FREE(FMMV, FMMV->P_Lriri2rrii, ((pL+1)*(pL+1))*sizeof(int));
	FMMV_FREE(FMMV, FMMV->P_X_riri2rrii, s_exp*sizeof(int));
	FMMV_FREE(FMMV, FMMV->P_VF, 2*FMMV->len_F*sizeof(int));
	FMMV_FREE(FMMV, FMMV->P_FV, (2*pL+1)*s_eps*sizeof(int));
	FMMV_FREE(FMMV, FMMV->neg_F, 2*(2*pM+1)*s_eps*sizeof(int));
	FMMV_FREE(FMMV, FMMV->neg_V, (2*pL+1)*s_eps*sizeof(int));
}
