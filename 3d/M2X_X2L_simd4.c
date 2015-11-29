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


extern _FLOAT_ zeros[];

static void CMX_simd4(FmmvHandle *FMMV, _FLOAT_ *x, _FLOAT_ *y)
{
	int p = FMMV->pM;
	int s_eps = FMMV->s_eps;
	_FLOAT_ *V = FMMV->CMX_coeffs;
	_FLOAT_ *xp = x;
	_FLOAT_ *yp = y;
	int m, k;
	
	k = p+1;
	gemv_trans_simd4(k, s_eps, V, s_eps, xp, yp);
	xp += 4*k;
	yp += 4*s_eps;
        V += s_eps*k;
	
	for (m=1; m<=p; m++) {
		k--;
		gemv_trans_simd4(k, s_eps, V, s_eps, xp, yp);
		xp += 4*k;
		yp += 4*s_eps;
		gemv_trans_simd4(k, s_eps, V, s_eps, xp, yp);
		xp += 4*k;
		yp += 4*s_eps;
                V += s_eps*k;
	}	
}


static void CXL_simd4(FmmvHandle *FMMV, _FLOAT_ *x, _FLOAT_ *y, int reduced)
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
	gemv_simd4(k, s_eps, V, s_eps, xp, yp);
	xp += 4*s_eps;
	yp += 4*k;
        V += s_eps*k;

	for (m=1; m<=p; m++) {
		k--;
		gemv_simd4(k, s_eps, V, s_eps, xp, yp);
		xp += 4*s_eps;
		yp += 4*k;
		gemv_simd4(k, s_eps, V, s_eps, xp, yp);
		xp += 4*s_eps;
		yp += 4*k;
                V += s_eps*k;
	}	
}	


static void CXL_minus_simd4(FmmvHandle *FMMV,  _FLOAT_ *x, _FLOAT_ *y, int reduced)
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
		gemv_simd4(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
		xp += 4*s_eps; yp += 4*k;
                V += s_eps*kk; kk--; /****/
		for (k=p/2; k>=1; k--) {
			gemv_simd4(k+1, s_eps, V, 2*s_eps, xp, yp);
			xp += 4*s_eps; yp += 4*(k+1);
			gemv_simd4(k , s_eps, V+s_eps, 2*s_eps, xp, yp);
			xp += 4*s_eps; yp += 4*k;
                        V += s_eps*kk; kk--; /****/
			gemv_simd4(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			xp += 4*s_eps; yp += 4*k;
			gemv_simd4(k, s_eps, V, 2*s_eps, xp, yp);
			xp += 4*s_eps; yp += 4*k;
                        V += s_eps*kk; kk--; /****/
		}
		
	}
	else {
		k = p/2;
		gemv_simd4(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
		xp += 4*s_eps; yp += 4*k;
                V += s_eps*kk; kk--; /****/
		gemv_simd4(k, s_eps, V, 2*s_eps, xp, yp);
		xp += 4*s_eps; yp += 4*k;
		gemv_simd4(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
		xp += 4*s_eps; yp += 4*k;
                V += s_eps*kk; kk--; /****/
		for (k=(p-1)/2; k>=1; k--) {
			gemv_simd4(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			xp += 4*s_eps; yp += 4*k;
			gemv_simd4(k+1, s_eps, V, 2*s_eps, xp, yp);
			xp += 4*s_eps; yp += 4*(k+1);
                        V += s_eps*kk; kk--; /****/
			gemv_simd4(k, s_eps, V, 2*s_eps, xp, yp);
			xp += 4*s_eps; yp += 4*k;
			gemv_simd4(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			xp += 4*s_eps; yp += 4*k;
                        V += s_eps*kk; kk--; /****/
		}
	}	
	gemv_simd4(1 , s_eps, V, 2*s_eps, xp, yp);

}



static void CXL_plus_simd4(FmmvHandle *FMMV, _FLOAT_ *x, _FLOAT_ *y, int reduced)
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
		gemv_simd4(k, s_eps, V, 2*s_eps, xp, yp);
		xp += 4*s_eps; yp += 4*k;
                V += s_eps*kk; kk--; /****/
		for (k=p/2; k>=1; k--) {
			gemv_simd4(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			xp += 4*s_eps; yp += 4*k;
			gemv_simd4(k+1 , s_eps, V, 2*s_eps, xp, yp);
			xp += 4*s_eps; yp += 4*(k+1);
                        V += s_eps*kk; kk--; /****/
			gemv_simd4(k, s_eps, V, 2*s_eps, xp, yp);
			xp += 4*s_eps; yp += 4*k;
			gemv_simd4(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			xp += 4*s_eps; yp += 4*k;
                        V += s_eps*kk; kk--; /****/
		}
	}
	else {
		k = p/2;
		gemv_simd4(k+1, s_eps, V, 2*s_eps, xp, yp);
		xp += 4*s_eps; yp += 4*(k+1);
                V += s_eps*kk; kk--; /****/
		gemv_simd4(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
		xp += 4*s_eps; yp += 4*k;
		gemv_simd4(k, s_eps, V, 2*s_eps, xp, yp);
		xp += 4*s_eps; yp += 4*k;
                V += s_eps*kk; kk--; /****/
		for (k=(p-1)/2; k>=1; k--) {
			gemv_simd4(k+1 , s_eps, V, 2*s_eps, xp, yp);
			xp += 4*s_eps; yp += 4*(k+1);
			gemv_simd4(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			xp += 4*s_eps; yp += 4*k;
                        V += s_eps*kk; kk--; /****/
			gemv_simd4(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			xp += 4*s_eps; yp += 4*k;
			gemv_simd4(k, s_eps, V, 2*s_eps, xp, yp);
			xp += 4*s_eps; yp += 4*k;
                        V += s_eps*kk; kk--; /****/
		}
	}	
	gemv_simd4(1 , s_eps, V, 2*s_eps, xp, yp);
}
	
static void CMX_minus_simd4(FmmvHandle *FMMV,  _FLOAT_ *x, _FLOAT_ *y)
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
		gemv_trans_simd4(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
		yp += 4*s_eps; xp += 4*k;
                V += s_eps*kk; kk--; /****/
		for (k=p/2; k>=1; k--) {
			gemv_trans_simd4(k+1, s_eps, V, 2*s_eps, xp, yp);
			yp += 4*s_eps; xp += 4*(k+1);
			gemv_trans_simd4(k , s_eps, V+s_eps, 2*s_eps, xp, yp);
			yp += 4*s_eps; xp += 4*k;
                        V += s_eps*kk; kk--; /****/
			gemv_trans_simd4(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			yp += 4*s_eps; xp += 4*k;
			gemv_trans_simd4(k, s_eps, V, 2*s_eps, xp, yp);
			yp += 4*s_eps; xp += 4*k;
                        V += s_eps*kk; kk--; /****/
		}
		gemv_trans_simd4(1, s_eps, V, 2*s_eps, xp, yp);
		memset(yp+4*s_eps, 0, 4*s_eps*sizeof(_FLOAT_));
		
	}
	else {
		k = p/2;
		gemv_trans_simd4(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
		yp += 4*s_eps; xp += 4*k;
                V += s_eps*kk; kk--; /****/
		gemv_trans_simd4(k, s_eps, V, 2*s_eps, xp, yp);
		yp += 4*s_eps; xp += 4*k;
		gemv_trans_simd4(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
		yp += 4*s_eps; xp += 4*k;
                V += s_eps*kk; kk--; /****/
		for (k=(p-1)/2; k>=1; k--) {
			gemv_trans_simd4(k , s_eps, V+s_eps, 2*s_eps, xp, yp);
			yp += 4*s_eps; xp += 4*k;
			gemv_trans_simd4(k+1, s_eps, V, 2*s_eps, xp, yp);
			yp += 4*s_eps; xp += 4*(k+1);
                        V += s_eps*kk; kk--; /****/
			gemv_trans_simd4(k, s_eps, V, 2*s_eps, xp, yp);
			yp += 4*s_eps; xp += 4*k;
			gemv_trans_simd4(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			yp += 4*s_eps; xp += 4*k;
                        V += s_eps*kk; kk--; /****/
		}
		memset(yp, 0, 4*s_eps*sizeof(_FLOAT_));
		gemv_trans_simd4(1, s_eps, V, 2*s_eps, xp, yp+4*s_eps);
	}	

}


	
static void CMX_plus_simd4(FmmvHandle *FMMV, _FLOAT_ *x, _FLOAT_ *y)
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
		gemv_trans_simd4(k, s_eps, V, 2*s_eps, xp, yp);
		yp += 4*s_eps; xp += 4*k;
                V += s_eps*kk; kk--; /****/
		for (k=p/2; k>=1; k--) {
			gemv_trans_simd4(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			yp += 4*s_eps; xp += 4*k;
			gemv_trans_simd4(k+1 , s_eps, V, 2*s_eps, xp, yp);
			yp += 4*s_eps; xp += 4*(k+1);
                        V += s_eps*kk; kk--; /****/
			gemv_trans_simd4(k, s_eps, V, 2*s_eps, xp, yp);
			yp += 4*s_eps; xp += 4*k;
			gemv_trans_simd4(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			yp += 4*s_eps; xp += 4*k;
                        V += s_eps*kk; kk--; /****/
		}
		memset(yp, 0, 4*s_eps*sizeof(_FLOAT_));
		gemv_trans_simd4(1, s_eps, V, 2*s_eps, xp, yp+4*s_eps);
	}
	else {
		k = p/2;
		gemv_trans_simd4(k+1, s_eps, V, 2*s_eps, xp, yp);
		yp += 4*s_eps; xp += 4*(k+1);
                V += s_eps*kk; kk--; /****/
		gemv_trans_simd4(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
		yp += 4*s_eps; xp += 4*k;
		gemv_trans_simd4(k, s_eps, V, 2*s_eps, xp, yp);
		yp += 4*s_eps; xp += 4*k;
                V += s_eps*kk; kk--; /****/
		for (k=(p-1)/2; k>=1; k--) {
			gemv_trans_simd4(k+1 , s_eps, V, 2*s_eps, xp, yp);
			yp += 4*s_eps; xp += 4*(k+1);
			gemv_trans_simd4(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			yp += 4*s_eps; xp += 4*k;
                        V += s_eps*kk; kk--; /****/
			gemv_trans_simd4(k, s_eps, V+s_eps, 2*s_eps, xp, yp);
			yp += 4*s_eps; xp += 4*k;
			gemv_trans_simd4(k, s_eps, V, 2*s_eps, xp, yp);
			yp += 4*s_eps; xp += 4*k;
                        V += s_eps*kk; kk--; /****/
		}
		gemv_trans_simd4(1 , s_eps, V, 2*s_eps, xp, yp);
		memset(yp+4*s_eps, 0, 4*s_eps*sizeof(_FLOAT_));
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
	int s_exp2 = s_exp/2;

	SIMD_ALIGN _FLOAT_ x1[4*(FMM_P_MAX+1)*(FMM_P_MAX+2)];
	SIMD_ALIGN _FLOAT_ x2[4*(FMM_P_MAX+1)*(FMM_P_MAX+2)];
	
	SIMD_ALIGN _FLOAT_ y1[4*2*(FMM_P_MAX+1)*FMM_S_EPS_MAX]; 
	SIMD_ALIGN _FLOAT_ y2[4*4*(FMM_P_MAX+1)*FMM_S_EPS_MAX];
	SIMD_ALIGN _FLOAT_ y3[4*2*(FMM_P_MAX+1)*FMM_S_EPS_MAX];
	
	SIMD_ALIGN _FLOAT_ z1[4*FMM_S_EXP_MAX];

	_FLOAT_ *XX1[8];
	_FLOAT_ *XX2[8];
	_FLOAT_ **X1_p[8] = {0,0,0,0,0,0,0,0};
	_FLOAT_ **X2_p[8] = {0,0,0,0,0,0,0,0};
	_FLOAT_ *M[8] = {zeros, zeros, zeros, zeros, zeros, zeros, zeros, zeros};

	int i, k, missing;
	int q[8]; 

	y2[4*len2] = 0.0;	
	y2[4*len2+1] = 0.0;		
	y2[4*len2+2] = 0.0;		
	y2[4*len2+3] = 0.0;		

	k = 0;
	for (i=0; i<8; i++) {
	   if (isSource(box->child[i])&&box->child[i]->M) {
		M[i] = box->child[i]->M;
		XX1[i] = X1 + i*s_exp;
		XX2[i] = X2 + i*s_exp;
		X1_p[i] = XX1+i;
		X2_p[i] = XX2+i;
		memset(X1 + i*s_exp, 0, s_exp*sizeof(_FLOAT_));
		memset(X2 + i*s_exp, 0, s_exp*sizeof(_FLOAT_));
		q[k]=i;
		k++;
	    }
	    else missing=i;
	}
	if (k==0) return;
	for (i=k; i<8; i++) q[i] = missing;
	
	for (i=0; i<=(k<=4?0:1); i++) {
		switch (dir) {
		case 0: /* UD */			
		        P_4riri2rrii_simd4(FMMV->pM, M[q[4*i]], M[q[4*i+1]], M[q[4*i+2]], M[q[4*i+3]], x1);
			break;
		case 1: /* NS */	
		        P_4riri2rrii_simd4(FMMV->pM, M[q[4*i]], M[q[4*i+1]], M[q[4*i+2]], M[q[4*i+3]], x1);
			Rz_minus_pi2_rrii_simd4(p, x1, x2);
			Ry_simd4(p, FMMV->Ry_pi2, x2, x1);
			break;
		case 2: /* EW */
		        P_4riri2rrii_simd4(FMMV->pM, M[q[4*i]], M[q[4*i+1]], M[q[4*i+2]], M[q[4*i+3]], x2);
			Ry_simd4(p, FMMV->Ry_minus_pi2, x2, x1);
			break;
		}	

                perm_simd4(len0p, FMMV->P_MRT_plus, x1, x2);
		CMX_plus_simd4(FMMV, x2, y1);
                perm_simd4(len0m, FMMV->P_MRT_minus, x1, x2);
		CMX_minus_simd4(FMMV, x2, y2);

		VEC_COPY(4*len2, y1, y3);
		VEC_ADD(4*len2, y2, y1, y1);
		VEC_SUB(4*len2, y3, y2, y3);
	

                y1[4*len2]=0;   y1[4*len2+1]=0;
                y1[4*len2+2]=0; y1[4*len2+3]=0;
		perm_simd4(2*FMMV->len_F, FMMV->P_VF, y1, y2);
		neg_simd4(FMMV->len_neg_F, FMMV->neg_F, y2);
		F_M2X_simd4(FMMV, y2, z1);
		P_X_riri2rrii4_simd4(s_exp2, z1,
			X1_p[q[4*i]],
			X1_p[q[4*i+1]],
			X1_p[q[4*i+2]],
			X1_p[q[4*i+3]]);
		

                y3[4*len2]=0;   y3[4*len2+1]=0;
                y3[4*len2+2]=0; y3[4*len2+3]=0;
		perm_simd4(2*FMMV->len_F, FMMV->P_VF, y3, y2);
		neg_simd4(FMMV->len_neg_F, FMMV->neg_F, y2);
		F_M2X_simd4(FMMV, y2, z1);
                P_X_riri2rrii4_simd4(s_exp2, z1,
			X2_p[q[4*i]],
			X2_p[q[4*i+1]],
			X2_p[q[4*i+2]],
			X2_p[q[4*i+3]]);
	}
}

void X2L(FmmvHandle *FMMV, int dir, Box *box, int reduced) 
{
	int p = FMMV->pL;
	int len0 = (p+1)*(p+1); 
        int len0m = len0/2;
        int len0p = len0 - len0m;
	int len2 = (2*p+1)*FMMV->s_eps; 
	int s_exp2 = FMMV->s_exp/2;

	SIMD_ALIGN _FLOAT_ x1[4*(FMM_P_MAX+1)*(FMM_P_MAX+2)];
	SIMD_ALIGN _FLOAT_ x2[4*(FMM_P_MAX+1)*(FMM_P_MAX+2)+12];
	SIMD_ALIGN _FLOAT_ xx[4*(FMM_P_MAX+1)*(FMM_P_MAX+2)];
	SIMD_ALIGN _FLOAT_ y1[4*4*(FMM_P_MAX+1)*FMM_S_EPS_MAX];
	SIMD_ALIGN _FLOAT_ y2[4*4*(FMM_P_MAX+1)*FMM_S_EPS_MAX];
	SIMD_ALIGN _FLOAT_ y3[4*2*(FMM_P_MAX+1)*FMM_S_EPS_MAX];
	SIMD_ALIGN _FLOAT_ z1[4*FMM_S_EXP_MAX];

	Box* child;
	_FLOAT_ *X1[8] = {zeros, zeros, zeros, zeros, zeros, zeros, zeros, zeros};
	_FLOAT_ *X2[8] = {zeros, zeros, zeros, zeros, zeros, zeros, zeros, zeros};
	_FLOAT_ **L_p[8] = {0,0,0,0,0,0,0,0};
	_FLOAT_ *Ad;	
	

	int i, k, k1, k2, missing, f;
	int q[8]; 

	if (!isTarget(box)) return;

	k = 0;
	for (i=0; i<8; i++) { 
	   child = box->child[i];
	   if (isTarget(child)) {
		L_p[i] = &(child->L);
		f = 0;
		if (reduced) {
			switch (dir) {
			case 0:	/* UD */
				if (child->X2[XU]) {X1[i] = child->X2[XU]; f=1;}
				if (child->X2[XD]) {X2[i] = child->X2[XD]; f=1;}
				break;
			case 1:	/* NS */
				if (child->X2[XN]) {X1[i] = child->X2[XN]; f=1;}
				if (child->X2[XS]) {X2[i] = child->X2[XS]; f=1;}
				break;
			case 2:	/* EW */
				if (child->X2[XE]) {X1[i] = child->X2[XE]; f=1;}
				if (child->X2[XW]) {X2[i] = child->X2[XW]; f=1;}
				break;
			}
		}
		else {
			switch (dir) {
			case 0:	/* UD */
				if (child->X[XU]) {X1[i] = child->X[XU]; f=1;}
				if (child->X[XD]) {X2[i] = child->X[XD]; f=1;}
				break;
			case 1:	/* NS */
				if (child->X[XN]) {X1[i] = child->X[XN]; f=1;}
				if (child->X[XS]) {X2[i] = child->X[XS]; f=1;}
				break;
			case 2:	/* EW */
				if (child->X[XE]) {X1[i] = child->X[XE]; f=1;}
				if (child->X[XW]) {X2[i] = child->X[XW]; f=1;}
				break;
			}
		}
		if (f) {
			q[k] = i;
			k++;
		}	
		else {
			missing = i;
		}	
	    }
   	    else missing=i;	   
	}

	if (k==0) return;
	for (i=k; i<8; i++) q[i] = missing;
	
	for (i=0; i<=(k<=4?0:1); i++) {
		k1 = (X1[q[4*i]]!=zeros)||(X1[q[4*i+1]]!=zeros)||(X1[q[4*i+2]]!=zeros)||(X1[q[4*i+3]]!=zeros);
		k2 = (X2[q[4*i]]!=zeros)||(X2[q[4*i+1]]!=zeros)||(X2[q[4*i+2]]!=zeros)||(X2[q[4*i+3]]!=zeros);
		if (k1&&k2) {
			P_X_4rrii2riri_simd4(s_exp2, X1[q[4*i]], X1[q[4*i+1]], X1[q[4*i+2]], X1[q[4*i+3]],z1);
			F_X2L_simd4(FMMV, z1, y2);
			perm_simd4(len2, FMMV->P_FV, y2, y1);
			neg_simd4(FMMV->len_neg_V, FMMV->neg_V, y1);
                        
                        P_X_4rrii2riri_simd4(s_exp2, X2[q[4*i]], X2[q[4*i+1]], X2[q[4*i+2]], X2[q[4*i+3]],z1);
			F_X2L_simd4(FMMV, z1, y2);
			perm_simd4(len2, FMMV->P_FV, y2, y3);
			neg_simd4(FMMV->len_neg_V, FMMV->neg_V, y3);
                        
                        
			VEC_COPY(4*len2, y1, y2);
			VEC_ADD(4*len2, y3, y1, y1);
			VEC_SUB(4*len2, y2, y3, y2);
			CXL_plus_simd4(FMMV, y1, x1, reduced);
			CXL_minus_simd4(FMMV, y2, x2, reduced);
			perm_inv_simd4(len0p, FMMV->P_LRT_plus, x1, xx);
			perm_inv_simd4(len0m, FMMV->P_LRT_minus, x2, xx);
		}
		else if (k1) {
			P_X_4rrii2riri_simd4(s_exp2, X1[q[4*i]], X1[q[4*i+1]], X1[q[4*i+2]], X1[q[4*i+3]],z1);
			F_X2L_simd4(FMMV, z1, y1);
			perm_simd4(len2, FMMV->P_FV, y1, y2);
			neg_simd4(FMMV->len_neg_V, FMMV->neg_V, y2);
			CXL_simd4(FMMV, y2, x1, reduced);
			perm_inv_simd4(len0, FMMV->P_LRT, x1, xx);
		}
		else if (k2) {
			P_X_4rrii2riri_simd4(s_exp2, X2[q[4*i]], X2[q[4*i+1]], X2[q[4*i+2]], X2[q[4*i+3]],z1);
			F_X2L_simd4(FMMV, z1, y1);
			perm_simd4(len2, FMMV->P_FV, y1, y2);
			neg_simd4(FMMV->len_neg_V, FMMV->neg_V, y2);
			CXL_simd4(FMMV, y2, x1, reduced);
			perm_inv_simd4(len0, FMMV->P_LRT, x1, xx);
			Ry_pi_simd4(p, xx);
		}
                if (k1||k2) {
                        scale_X_simd4(p, xx, box->level+1);

		        switch (dir) {
		        case 0:	 /* UD */
		                P_rrii2riri4_simd4(FMMV, p, xx, L_p[q[4*i]], L_p[q[4*i+1]], L_p[q[4*i+2]], L_p[q[4*i+3]], x1);
			        break;
		        case 1:	/* NS */
			        Ry_simd4(p, FMMV->Ry_minus_pi2,xx, x1);
			        Rz_pi2_rrii_simd4(p, x1, xx);
		                P_rrii2riri4_simd4(FMMV, p, xx, L_p[q[4*i]], L_p[q[4*i+1]], L_p[q[4*i+2]], L_p[q[4*i+3]], x1);
			        break;
		        case 2: /* EW */	
			        Ry_simd4(p, FMMV->Ry_pi2, xx, x1);
		                P_rrii2riri4_simd4(FMMV, p, x1, L_p[q[4*i]], L_p[q[4*i+1]], L_p[q[4*i+2]], L_p[q[4*i+3]], xx);
			        break;
		        }
                }
                
	}	
}
