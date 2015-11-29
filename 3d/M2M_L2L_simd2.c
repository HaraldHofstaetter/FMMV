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

#include "simd.h"
#if (FMM_PRECISION==0)
   #include"simd2s.h"
#else
   #include"simd2d.h"
#endif

extern _FLOAT_ zeros[];

static void Tz_M2M_simd2(FmmvHandle *FMMV, _FLOAT_ *x)
{
	int p = FMMV->pM; 
	_FLOAT_ **Tz = FMMV->Tz_M2M;
	int k, dk, kk;

        if (FMMV->beta==0) {
	    dk = p+1;
	    kk = p+1;
	    tpmv_lower_simd2(dk, Tz[0], x);
	    for (k=1; k<=p; k++) {
		    dk--;
		    tpmv_lower_simd2(dk, Tz[k], x + 2*kk);	
		    kk += dk;
		    tpmv_lower_simd2(dk, Tz[k], x + 2*kk);	
		    kk += dk;
	    }	
        }
        else {
            _FLOAT_ y[2*(FMM_P_MAX+1)*(FMM_P_MAX+2)];
	    dk = p+1;
	    kk = p+1;
	    gemv_simd2(dk, dk, Tz[0], dk, x, y);
	    for (k=1; k<=p; k++) {
		    dk--;
		    gemv_simd2(dk, dk, Tz[k], dk, x + 2*kk, y + 2*kk);	
		    kk += dk;
		    gemv_simd2(dk, dk, Tz[k], dk, x + 2*kk, y + 2*kk);	
		    kk += dk;
	    }	
            memcpy(x, y, 2*sizeof(_FLOAT_)*(p+1)*(p+2));
        }
}	

static void Tz_L2L_simd2(FmmvHandle *FMMV, _FLOAT_ *x)
{
	int p = FMMV->pL; 
	_FLOAT_ **Tz = FMMV->Tz_L2L;
	int k, dk, kk;

        if (FMMV->beta==0) {
	    dk = p+1;
	    kk = p+1;
	    tpmv_upper_simd2(dk, Tz[0], x);
	    for (k=1; k<=p; k++) {
		    dk--;
		    tpmv_upper_simd2(dk, Tz[k], x + 2*kk);	
		    kk += dk;
		    tpmv_upper_simd2(dk, Tz[k], x + 2*kk);	
		    kk += dk;
	    }	
	}	
        else {
            _FLOAT_ y[2*(FMM_P_MAX+1)*(FMM_P_MAX+2)];
	    dk = p+1;
	    kk = p+1;
	    gemv_simd2(dk, dk, Tz[0], dk, x, y);
	    for (k=1; k<=p; k++) {
		    dk--;
		    gemv_simd2(dk, dk, Tz[k], dk, x + 2*kk, y + 2*kk);	
		    kk += dk;
		    gemv_simd2(dk, dk, Tz[k], dk, x + 2*kk, y + 2*kk);	
		    kk += dk;
	    }	
            memcpy(x, y, 2*sizeof(_FLOAT_)*(p+1)*(p+2));
        }
}


void M2M(FmmvHandle *FMMV, Box *box)
{
        int p = FMMV->pM;
	int len0 = (p+1)*(p+1);
	int len = (p+1)*(p+2);
        int *P_RT = FMMV->P_MRT;
        int *P_riri2rrii = FMMV->P_Mriri2rrii;
	
        SIMD_ALIGN _FLOAT_ x1[2*(FMM_P_MAX+1)*(FMM_P_MAX+2)];
        SIMD_ALIGN _FLOAT_ x2[2*(FMM_P_MAX+1)*(FMM_P_MAX+2)];
        SIMD_ALIGN _FLOAT_ xx[2*(FMM_P_MAX+1)*(FMM_P_MAX+2)];

	_FLOAT_ *M1[4] = {zeros, zeros, zeros, zeros};
	_FLOAT_ *M2[4] = {zeros, zeros, zeros, zeros};

	int i, k, k1, k2, missing, f;
	int q[4];

        if (!isSource(box)) return;

	k=0;
	for (i=0; i<4; i++) {
		f = 0;
		if (isSource(box->child[i])&&box->child[i]->M) {
			f = 1;
			M1[i] = box->child[i]->M;
		}	
		if (isSource(box->child[7-i])&&box->child[7-i]->M) {
			/* i-7 ... child diagonally opposite to child i */
			f = 1;
			M2[i] = box->child[7-i]->M;
		}	
		if (f) {
			q[k] = i;
			k++;
		}
		else {	
			missing = i;
		}
	}
	if (k==0) return;
	for (i=k; i<4; i++) q[i] = missing;

	for (i=0; i<=(k<=2?0:1); i++) {
		k1 = (M1[q[2*i]]!=zeros)||(M1[q[2*i+1]]!=zeros);
		k2 = (M2[q[2*i]]!=zeros)||(M2[q[2*i+1]]!=zeros);
		if (k1) { 
			vec2_copy_simd2(len, M1[q[2*i]], M1[q[2*i+1]], x1); 
			Rz_simd2(p, x1, x2, q[2*i], q[2*i+1], 0);
			perm_simd2(len0, P_riri2rrii, x2, x1);
			Ry_simd2(p, FMMV->Ry_pi_minus_theta, x1, x2);
			perm_simd2(len0, P_RT, x2, x1);
			Tz_M2M_simd2(FMMV, x1);
			perm_inv_simd2(len0, P_RT, x1,xx);
			if (k2) { 
				vec2_copy_simd2(len, M2[q[2*i]], M2[q[2*i+1]], x1); 
				Rz_simd2(p, x1, x2, q[2*i], q[2*i+1], 0);
				perm_simd2(len0, P_riri2rrii, x2, x1);
				Ry_simd2(p, FMMV->Ry_minus_theta, x1, x2);
				perm_simd2(len0, P_RT, x2,x1);
				Tz_M2M_simd2(FMMV, x1);
				perm_inv_simd2(len0, P_RT, x1,x2);
				Ry_pi_simd2(p, x2);
				VEC_ADD(2*len, x2, xx, xx);
			}
			Ry_simd2(p, FMMV->Ry_minus_pi_minus_theta, xx, x1);
			perm_inv_simd2(len0, P_riri2rrii, x1, x2);
			Rz_simd2(p, x2, x1, q[2*i], q[2*i+1], 1);
			box->M = vec2_add2_simd2(FMMV, len, x1, box->M);
		}
		else if (k2) {
			vec2_copy_simd2(len, M2[q[2*i]], M2[q[2*i+1]], x1); 
			Rz_simd2(p, x1, x2, q[2*i], q[2*i+1], 0);
			perm_simd2(len0, P_riri2rrii, x2, x1);
			Ry_simd2(p, FMMV->Ry_minus_theta, x1, x2);
			perm_simd2(len0, P_RT, x2,x1);
			Tz_M2M_simd2(FMMV, x1);
			perm_inv_simd2(len0, P_RT, x1,x2);
			Ry_simd2(p, FMMV->Ry_theta, x2, x1);
			perm_inv_simd2(len0, P_riri2rrii, x1, x2);
			Rz_simd2(p, x2, x1, q[2*i], q[2*i+1], 1);
			box->M = vec2_add2_simd2(FMMV, len, x1, box->M);
		}	
	}	
}


void L2L(FmmvHandle *FMMV, Box *box)
{
	int p = FMMV->pL;
	int len0 = (p+1)*(p+1);
	int len = (p+1)*(p+2);
        int *P_RT = FMMV->P_LRT;
        int *P_riri2rrii = FMMV->P_Lriri2rrii;

        SIMD_ALIGN _FLOAT_ x1[2*(FMM_P_MAX+1)*(FMM_P_MAX+2)];
        SIMD_ALIGN _FLOAT_ x2[2*(FMM_P_MAX+1)*(FMM_P_MAX+2)]; 
        SIMD_ALIGN _FLOAT_ xx[2*(FMM_P_MAX+1)*(FMM_P_MAX+2)];

        _FLOAT_ **L1_p[4] = {0,0,0,0};
        _FLOAT_ **L2_p[4] = {0,0,0,0};

	int i, k, k1, k2, missing, f;
	int q[4]; 

        if (!isTarget(box)||!box->L) return;

	k=0;
	for (i=0; i<4; i++) {
		f = 0;
		if (isTarget(box->child[i])) {
			f = 1;
			L1_p[i] = &(box->child[i]->L);
		}	
		if (isTarget(box->child[7-i])) {
			/* i-7 ... child diagonally opposite to child i */
			f = 1;
			L2_p[i] = &(box->child[7-i]->L);
		}	
		if (f) {
			q[k] = i;
			k++;
		}
		else {	
			missing = i;
		}
	}
	if (k==0) return;
	for (i=k; i<4; i++) q[i] = missing;
	for (i=0; i<=(k<=2?0:1); i++) {
		k1 = (L1_p[q[2*i]]!=0)||(L1_p[q[2*i+1]]!=0);
		k2 = (L2_p[q[2*i]]!=0)||(L2_p[q[2*i+1]]!=0);
		if (k2) {
			Rz1_simd2(p, box->L, x2, q[2*i], q[2*i+1], 0);
			perm_simd2(len0, P_riri2rrii,x2, x1);
			Ry_simd2(p, FMMV->Ry_pi_minus_theta, x1, xx);
			perm_simd2(len0, P_RT, xx, x1);
			Tz_L2L_simd2(FMMV, x1);
			perm_inv_simd2(len0, P_RT, x1, x2);
			Ry_simd2(p, FMMV->Ry_minus_pi_minus_theta, x2, x1);
			perm_inv_simd2(len0, P_riri2rrii, x1, x2);
			Rz_simd2(p, x2, x1, q[2*i], q[2*i+1], 1);
			VEC_2add2_simd2(FMMV, len, x1, L2_p[q[2*i]], L2_p[q[2*i+1]]); 
			if (k1) {
				Ry_pi_simd2(p, xx);
				perm_simd2(len0, P_RT,xx,x1);
				Tz_L2L_simd2(FMMV, x1);
				perm_inv_simd2(len0, P_RT, x1, x2);
				Ry_simd2(p, FMMV->Ry_theta, x2, x1);
				perm_inv_simd2(len0, P_riri2rrii, x1, x2);
				Rz_simd2(p, x2, x1, q[2*i], q[2*i+1], 1);
				VEC_2add2_simd2(FMMV, len, x1, L1_p[q[2*i]], L1_p[q[2*i+1]]); 
			}
		}
		else if (k1) {
			Rz1_simd2(p, box->L, x2, q[2*i], q[2*i+1], 0);
			perm_simd2(len0, P_riri2rrii, x2, x1);
			Ry_simd2(p, FMMV->Ry_minus_theta, x1, x2);
			perm_simd2(len0, P_RT, x2,x1);
			Tz_L2L_simd2(FMMV, x1);
			perm_inv_simd2(len0, P_RT, x1, x2);
			Ry_simd2(p, FMMV->Ry_theta, x2, x1);
			perm_inv_simd2(len0, P_riri2rrii, x1, x2);
			Rz_simd2(p, x2, x1, q[2*i], q[2*i+1],1);
			VEC_2add2_simd2(FMMV, len, x1, L1_p[q[2*i]], L1_p[q[2*i+1]]); 
		}
	}
}
