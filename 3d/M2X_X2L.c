/*
 * FMMV - the Fastest Multipole Method of Vienna
 * Copyright (c) 2006-2010 Harald Hofstaetter
 * Institute of Mathematics
 * University of Vienna, Austria
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

static void CMX(FmmvHandle *FMMV, _FLOAT_ *x, _FLOAT_ *y)
{
	int p = FMMV->pM;
	int s_eps = FMMV->s_eps;
	_FLOAT_ *V = FMMV->CMX_coeffs;
	_FLOAT_ *xp = x;
	_FLOAT_ *yp = y;
	int m, k;
	
	k = p+1;
	gemv_trans(k, s_eps, V, s_eps, xp, yp);
	xp += k;
	yp += s_eps;
        V += s_eps*k;
	
	for (m=1; m<=p; m++) {
		k--;
		gemv_trans(k, s_eps, V, s_eps, xp, yp);
		xp += k;
		yp += s_eps;
		gemv_trans(k, s_eps, V, s_eps, xp, yp);
		xp += k;
		yp += s_eps;
                V += s_eps*k;
	}	
}

static void CXL(FmmvHandle *FMMV, _FLOAT_ *x, _FLOAT_ *y, int reduced)
{
	int p = FMMV->pL;
	int s_eps = FMMV->s_eps;
	_FLOAT_ *V;
	_FLOAT_ *xp = x;
	_FLOAT_ *yp = y;
	int m, k;

        if (reduced) {
	    V = FMMV->CXL_reduced_coeffs;
        }
        else {
	    V = FMMV->CXL_coeffs;
        }

	k = p+1;
	gemv(k, s_eps, V, s_eps, xp, yp);
	xp += s_eps;
	yp += k;
        V += s_eps*k;

	for (m=1; m<=p; m++) {
		k--;
		gemv(k, s_eps, V, s_eps, xp, yp);
		xp += s_eps;
		yp += k;
		gemv(k, s_eps, V, s_eps, xp, yp);
		xp += s_eps;
		yp += k;
                V += s_eps*k;
	}	
}	

static void CXL_minus(FmmvHandle *FMMV,  _FLOAT_ *x, _FLOAT_ *y, int reduced)
{
	int p = FMMV->pL;
	int s_eps = FMMV->s_eps;
	_FLOAT_*V;
	_FLOAT_ *xp = x;
	_FLOAT_ *yp = y;
	int k, kk;

        if (reduced) {
	    V = FMMV->CXL_reduced_coeffs;
        }
        else {
	    V = FMMV->CXL_coeffs;
        }

        kk = p+1;
	if (p%2) {
		k = (p+1)/2;
		gemv(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
		xp += s_eps; yp += k;
                V += s_eps*kk; kk--; /****/
		for (k=p/2; k>=1; k--) {
			gemv(k+1, s_eps, V, 2*s_eps, xp, yp);
			xp += s_eps; yp += k+1;
			gemv(k , s_eps, V+s_eps, 2*s_eps, xp, yp);
			xp += s_eps; yp += k;
                        V += s_eps*kk; kk--; /****/
			gemv(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			xp += s_eps; yp += k;
			gemv(k, s_eps, V, 2*s_eps, xp, yp);
			xp += s_eps; yp += k;
                        V += s_eps*kk; kk--; /****/
		}
		
	}
	else {
		k = p/2;
		gemv(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
		xp += s_eps; yp += k;
                V += s_eps*kk; kk--; /****/
		gemv(k, s_eps, V, 2*s_eps, xp, yp);
		xp += s_eps; yp += k;
		gemv(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
		xp += s_eps; yp += k;
                V += s_eps*kk; kk--; /****/
		for (k=(p-1)/2; k>=1; k--) {
			gemv(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			xp += s_eps; yp += k;
			gemv(k+1, s_eps, V, 2*s_eps, xp, yp);
			xp += s_eps; yp += k+1;
                        V += s_eps*kk; kk--; /****/
			gemv(k, s_eps, V, 2*s_eps, xp, yp);
			xp += s_eps; yp += k;
			gemv(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			xp += s_eps; yp += k;
                        V += s_eps*kk; kk--; /****/
		}
	}	
	gemv(1 , s_eps, V, 2*s_eps, xp, yp);

}



static void CXL_plus(FmmvHandle *FMMV, _FLOAT_ *x, _FLOAT_ *y, int reduced)
{
	int p = FMMV->pL;
	int s_eps = FMMV->s_eps;
	_FLOAT_ *V;
	_FLOAT_ *xp = x;
	_FLOAT_ *yp = y;
	int k, kk;

        if (reduced) {
	    V = FMMV->CXL_reduced_coeffs;
        }
        else {
	    V = FMMV->CXL_coeffs;
        }
	
        kk = p+1;
	if (p%2) {
		k = (p+1)/2;
		gemv(k, s_eps, V, 2*s_eps, xp, yp);
		xp += s_eps; yp += k;
                V += s_eps*kk; kk--; /****/
		for (k=p/2; k>=1; k--) {
			gemv(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			xp += s_eps; yp += k;
			gemv(k+1 , s_eps, V, 2*s_eps, xp, yp);
			xp += s_eps; yp += k+1;
                        V += s_eps*kk; kk--; /****/
			gemv(k, s_eps, V, 2*s_eps, xp, yp);
			xp += s_eps; yp += k;
			gemv(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			xp += s_eps; yp += k;
                        V += s_eps*kk; kk--; /****/
		}
	}
	else {
		k = p/2;
		gemv(k+1, s_eps, V, 2*s_eps, xp, yp);
		xp += s_eps; yp += k+1;
                V += s_eps*kk; kk--; /****/
		gemv(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
		xp += s_eps; yp += k;
		gemv(k, s_eps, V, 2*s_eps, xp, yp);
		xp += s_eps; yp += k;
                V += s_eps*kk; kk--; /****/
		for (k=(p-1)/2; k>=1; k--) {
			gemv(k+1 , s_eps, V, 2*s_eps, xp, yp);
			xp += s_eps; yp += k+1;
			gemv(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			xp += s_eps; yp += k;
                        V += s_eps*kk; kk--; /****/
			gemv(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			xp += s_eps; yp += k;
			gemv(k, s_eps, V, 2*s_eps, xp, yp);
			xp += s_eps; yp += k;
                        V += s_eps*kk; kk--; /****/
		}
	}	
	gemv(1 , s_eps, V, 2*s_eps, xp, yp);
}
	
static void CMX_minus(FmmvHandle *FMMV,  _FLOAT_ *x, _FLOAT_ *y)
{
	int p = FMMV->pM;
	int s_eps = FMMV->s_eps;
	_FLOAT_ *V = FMMV->CMX_coeffs;
	_FLOAT_ *xp = x;
	_FLOAT_ *yp = y;
	int k, kk;

        kk = p+1;
	if (p%2) {
		k = (p+1)/2;
		gemv_trans(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
		yp += s_eps; xp += k;
                V += s_eps*kk; kk--; /****/
		for (k=p/2; k>=1; k--) {
			gemv_trans(k+1, s_eps, V, 2*s_eps, xp, yp);
			yp += s_eps; xp += k+1;
			gemv_trans(k , s_eps, V+s_eps, 2*s_eps, xp, yp);
			yp += s_eps; xp += k;
                        V += s_eps*kk; kk--; /****/
			gemv_trans(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			yp += s_eps; xp += k;
			gemv_trans(k, s_eps, V, 2*s_eps, xp, yp);
			yp += s_eps; xp += k;
                        V += s_eps*kk; kk--; /****/
		}
		gemv_trans(1, s_eps, V, 2*s_eps, xp, yp);
		memset(yp+s_eps, 0, s_eps*sizeof(_FLOAT_));
		
	}
	else {
		k = p/2;
		gemv_trans(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
		yp += s_eps; xp += k;
                V += s_eps*kk; kk--; /****/
		gemv_trans(k, s_eps, V, 2*s_eps, xp, yp);
		yp += s_eps; xp += k;
		gemv_trans(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
		yp += s_eps; xp += k;
                V += s_eps*kk; kk--; /****/
		for (k=(p-1)/2; k>=1; k--) {
			gemv_trans(k , s_eps, V+s_eps, 2*s_eps, xp, yp);
			yp += s_eps; xp += k;
			gemv_trans(k+1, s_eps, V, 2*s_eps, xp, yp);
			yp += s_eps; xp += k+1;
                        V += s_eps*kk; kk--; /****/
			gemv_trans(k, s_eps, V, 2*s_eps, xp, yp);
			yp += s_eps; xp += k;
			gemv_trans(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			yp += s_eps; xp += k;
                        V += s_eps*kk; kk--; /****/
		}
		memset(yp, 0, s_eps*sizeof(_FLOAT_));
		gemv_trans(1, s_eps, V, 2*s_eps, xp, yp+s_eps);
	}	

}


	
static void CMX_plus(FmmvHandle *FMMV, _FLOAT_ *x, _FLOAT_ *y)
{
	int p = FMMV->pM;
	int s_eps = FMMV->s_eps;
	_FLOAT_ *V = FMMV->CMX_coeffs;
	_FLOAT_ *xp = x;
	_FLOAT_ *yp = y;
	int k, kk;
	
        kk = p+1;
	if (p%2) {
		k = (p+1)/2;
		gemv_trans(k, s_eps, V, 2*s_eps, xp, yp);
		yp += s_eps; xp += k;
                V += s_eps*kk; kk--; /****/
		for (k=p/2; k>=1; k--) {
			gemv_trans(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			yp += s_eps; xp += k;
			gemv_trans(k+1 , s_eps, V, 2*s_eps, xp, yp);
			yp += s_eps; xp += k+1;
                        V += s_eps*kk; kk--; /****/
			gemv_trans(k, s_eps, V, 2*s_eps, xp, yp);
			yp += s_eps; xp += k;
			gemv_trans(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			yp += s_eps; xp += k;
                        V += s_eps*kk; kk--; /****/
		}
		memset(yp, 0, s_eps*sizeof(_FLOAT_));
		gemv_trans(1, s_eps, V, 2*s_eps, xp, yp+s_eps);
	}
	else {
		k = p/2;
		gemv_trans(k+1, s_eps, V, 2*s_eps, xp, yp);
		yp += s_eps; xp += (k+1);
                V += s_eps*kk; kk--; /****/
		gemv_trans(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
		yp += s_eps; xp += k;
		gemv_trans(k, s_eps, V, 2*s_eps, xp, yp);
		yp += s_eps; xp += k;
                V += s_eps*kk; kk--; /****/
		for (k=(p-1)/2; k>=1; k--) {
			gemv_trans(k+1 , s_eps, V, 2*s_eps, xp, yp);
			yp += s_eps; xp += k+1;
			gemv_trans(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			yp += s_eps; xp += k;
                        V += s_eps*kk; kk--; /****/
			gemv_trans(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			yp += s_eps; xp += k;
			gemv_trans(k, s_eps, V, 2*s_eps, xp, yp);
			yp += s_eps; xp += k;
                        V += s_eps*kk; kk--; /****/
		}
		gemv_trans(1 , s_eps, V, 2*s_eps, xp, yp);
		memset(yp+s_eps, 0, s_eps*sizeof(_FLOAT_));
	}	
}
	

void M2X(FmmvHandle *FMMV, int dir, Box *box, _FLOAT_ *X1, _FLOAT_ *X2)
{
	int p = FMMV->pM;
	int len0 = (p+1)*(p+1); 
	int len0m = len0/2;
	int len0p = len0 - len0m;
	int len2 = (2*p+1)*FMMV->s_eps; 
	int s_exp = FMMV->s_exp;

        _FLOAT_ x1[(FMM_P_MAX+1)*(FMM_P_MAX+2)];
        _FLOAT_ x2[(FMM_P_MAX+1)*(FMM_P_MAX+2)];

        _FLOAT_ y1[2*(FMM_P_MAX+1)*FMM_S_EPS_MAX];
        _FLOAT_ y2[2*(FMM_P_MAX+1)*FMM_S_EPS_MAX];
        _FLOAT_ y3[2*(FMM_P_MAX+1)*FMM_S_EPS_MAX];

        _FLOAT_ z1[FMM_S_EXP_MAX];
        int i;


        for (i=0; i<8; i++) if (isSource(box->child[i])&&box->child[i]->M) {
		switch (dir) {
		case 0:	/* UD */
                	perm(len0, FMMV->P_Mriri2rrii, box->child[i]->M, x1);

			break;
		case 1:	/* NS */
			Rz_minus_pi2(p, box->child[i]->M, x1);
			perm(len0, FMMV->P_Mriri2rrii, x1, x2);
			Ry(p, FMMV->Ry_pi2, x2, x1);
			break;
		case 2:	/* EW */
                	perm(len0, FMMV->P_Mriri2rrii, box->child[i]->M, x2);
			Ry(p, FMMV->Ry_minus_pi2, x2, x1);
			break;
		}	

                perm(len0p, FMMV->P_MRT_plus, x1, x2);
                CMX_plus(FMMV, x2, y1);
                perm(len0m, FMMV->P_MRT_minus, x1, x2);
                CMX_minus(FMMV, x2, y2);

                VEC_COPY(len2, y1, y3);
                VEC_ADD(len2, y2, y1, y1);
                VEC_SUB(len2, y3, y2, y3);
                

	        y1[len2] = 0.0;	
                perm(2*FMMV->len_F, FMMV->P_VF, y1, y2);
                neg(FMMV->len_neg_F, FMMV->neg_F, y2);
#ifdef X_RRII
		F_M2X(FMMV, y2, z1);
		perm(s_exp, FMMV->P_X_riri2rrii, z1, X1 + i*s_exp);
#else
		F_M2X(FMMV, y2, X1 + i*s_exp);
#endif
		
	        y3[len2] = 0.0;	
                perm(2*FMMV->len_F, FMMV->P_VF, y3, y2);
                neg(FMMV->len_neg_F, FMMV->neg_F, y2);
#ifdef X_RRII
		F_M2X(FMMV, y2, z1);
		perm(s_exp, FMMV->P_X_riri2rrii, z1, X2 + i*s_exp);
#else
		F_M2X(FMMV, y2, X2 + i*s_exp);
#endif
        }
}
		



void X2L(FmmvHandle *FMMV, int dir, Box *box, int reduced) 
{
	int p = FMMV->pL;
	int len0 = (p+1)*(p+1); 
        int len0m = len0/2;
        int len0p = len0 - len0m;
	int len1 = (p+1)*(p+2);
	int len2 = (2*p+1)*FMMV->s_eps; 
        int s_exp = FMMV->s_exp;

	_FLOAT_ x1[(FMM_P_MAX+1)*(FMM_P_MAX+2)];
	_FLOAT_ x2[(FMM_P_MAX+1)*(FMM_P_MAX+2)];
	_FLOAT_ xx[(FMM_P_MAX+1)*(FMM_P_MAX+2)];
	_FLOAT_ y1[2*(FMM_P_MAX+1)*FMM_S_EPS_MAX];
	_FLOAT_ y2[4*(FMM_P_MAX+1)*FMM_S_EPS_MAX]; 
//CHECK: how long must y2 be eaxtly???
	_FLOAT_ y3[2*(FMM_P_MAX+1)*FMM_S_EPS_MAX];
#ifdef X_RRII
	_FLOAT_ z1[FMM_S_EXP_MAX];
#endif

	Box* child;
	int i;
	_FLOAT_ *X1;
	_FLOAT_ *X2;
	int n, j, jj;
        int level = box->level+1;
		
	if (!isTarget(box)) return;

	for (i=0; i<8; i++) {
	    child = box->child[i];
	    if (isTarget(child)) {
		if (reduced) {
			switch (dir) {
			case 0:	/* UD */
				X1 = child->X2[XU];
				X2 = child->X2[XD];
				break;
			case 1:	/* NS */
				X1 = child->X2[XN];
				X2 = child->X2[XS];
				break;
			case 2:	/* EW */
				X1 = child->X2[XE];
				X2 = child->X2[XW];
				break;
			}
		}
		else {
			switch (dir) {
			case 0:	/* UD */
				X1 = child->X[XU];
				X2 = child->X[XD];
				break;
			case 1:	/* NS */
				X1 = child->X[XN];
				X2 = child->X[XS];
				break;
			case 2:	/* EW */
				X1 = child->X[XE];
				X2 = child->X[XW];
				break;
			}
		}
		if (X1&&X2) {
			#ifdef X_RRII
			perm_inv(s_exp, FMMV->P_X_riri2rrii, X1, z1);
			F_X2L(FMMV, z1, y2);
			#else
			F_X2L(FMMV, X1, y2);
			#endif
			perm(len2, FMMV->P_FV, y2, y1);
			neg(FMMV->len_neg_V, FMMV->neg_V, y1);
                        
			#ifdef X_RRII
			perm_inv(s_exp, FMMV->P_X_riri2rrii, X2, z1);
			F_X2L(FMMV, z1, y2);
			#else
			F_X2L(FMMV, X2, y2);
			#endif
			perm(len2, FMMV->P_FV, y2, y3);
			neg(FMMV->len_neg_V, FMMV->neg_V, y3);
                        
			VEC_COPY(len2, y1, y2);
			VEC_ADD(len2, y3, y1, y1);
			VEC_SUB(len2, y2, y3, y2);
			CXL_plus(FMMV, y1, x1, reduced);
			CXL_minus(FMMV, y2, x2, reduced);
			perm_inv(len0p, FMMV->P_LRT_plus, x1, xx);
			perm_inv(len0m, FMMV->P_LRT_minus, x2, xx);
		}
		else if (X1) {
			#ifdef X_RRII
			perm_inv(s_exp, FMMV->P_X_riri2rrii, X1, z1);
			F_X2L(FMMV, z1, y1);
			#else
			F_X2L(FMMV, X1, y1);
			#endif
			perm(len2, FMMV->P_FV, y1, y2);
			neg(FMMV->len_neg_V, FMMV->neg_V, y2);
			CXL(FMMV, y2, x1, reduced);
			perm_inv(len0, FMMV->P_LRT, x1, xx);
		}
		else if (X2) {
			#ifdef X_RRII
			perm_inv(s_exp, FMMV->P_X_riri2rrii, X2, z1);
			F_X2L(FMMV, z1, y1);
			#else
			F_X2L(FMMV, X2, y1);
			#endif
			perm(len2, FMMV->P_FV, y1, y2);
			neg(FMMV->len_neg_V, FMMV->neg_V, y2);
			CXL(FMMV, y2, x1, reduced);
			perm_inv(len0, FMMV->P_LRT, x1, xx);
			Ry_pi(p, xx);
		}
		if (X1||X2) {

                        jj = 0;
                        for (n=0;n<=p;n++) {
                           if(n&1) { /* n odd */
                               for(j=0;j<2*n+1;j++) { // NOTE: Check this for reduced_scheme !!!
                                   xx[jj] = -ldexp(xx[jj], level);
                                   jj++;
                               }
                           }
                           else {
                               for(j=0;j<2*n+1;j++) {
                                   xx[jj] = +ldexp(xx[jj], level);
                                   jj++;
                               }
                           }
                        } 
                          
			switch (dir) {
			case 0:	 /* UD */
				perm_inv(len0, FMMV->P_Lriri2rrii, xx, x2);
				jj = 1;
				for (j=0; j<=p; j++) {
					x2[jj] = 0.0;
					jj += 2*j;
				}	
				child->L = VEC_ADD2(FMMV, len1, x2, child->L);
				break;
			case 1:	/* NS */
				Ry(p, FMMV->Ry_minus_pi2, xx, x2);
				perm_inv(len0, FMMV->P_Lriri2rrii, x2, x1);
				Rz_pi2(p, x1, x2);
				jj = 1;
				for (j=0; j<=p; j++) {
					x2[jj] = 0.0;
					jj += 2*j;
				}	
				child->L = VEC_ADD2(FMMV, len1, x2, child->L);
				break;
			case 2: /* EW */	
				Ry(p, FMMV->Ry_pi2, xx, x2);
				perm_inv(len0, FMMV->P_Lriri2rrii, x2, x1);
				jj = 1;
				for (j=0; j<=p; j++) {
					x1[jj] = 0.0;
					jj += 2*j;
				}	
				child->L = VEC_ADD2(FMMV, len1, x1, child->L);
				break;
			}	
		}
	    }
	}
}

