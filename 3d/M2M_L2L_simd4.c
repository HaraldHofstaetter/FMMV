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
   #include"simd4s.h"
#else
   #include"simd4d.h"
#endif

static void Tz_M2M_simd4(FmmvHandle *FMMV, _FLOAT_ *x)
{
	int p = FMMV->pM; 
	_FLOAT_ **Tz = FMMV->Tz_M2M;
	int k, dk, kk;

        if (FMMV->beta==0) {
	    dk = p+1;
	    kk = p+1;
	    tpmv_lower_simd4(dk, Tz[0], x);
	    for (k=1; k<=p; k++) {
		    dk--;
		    tpmv_lower_simd4(dk, Tz[k], x + 4*kk);	
		    kk += dk;
		    tpmv_lower_simd4(dk, Tz[k], x + 4*kk);	
		    kk += dk;
	    }	
        }
        else {
            _FLOAT_ y[4*(FMM_P_MAX+1)*(FMM_P_MAX+2)];
	    dk = p+1;
	    kk = p+1;
	    gemv_simd4(dk, dk, Tz[0], dk, x, y);
	    for (k=1; k<=p; k++) {
		    dk--;
		    gemv_simd4(dk, dk, Tz[k], dk, x + 4*kk, y + 4*kk);	
		    kk += dk;
		    gemv_simd4(dk, dk, Tz[k], dk, x + 4*kk, y + 4*kk);	
		    kk += dk;
	    }	
            memcpy(x, y, 4*sizeof(_FLOAT_)*(p+1)*(p+2));
        }
}	

static void Tz_L2L_simd4(FmmvHandle *FMMV, _FLOAT_ *x)
{
	int p = FMMV->pL; 
	_FLOAT_ **Tz = FMMV->Tz_L2L;
	int k, dk, kk;

        if (FMMV->beta==0) {
	    dk = p+1;
	    kk = p+1;
	    tpmv_upper_simd4(dk, Tz[0], x);
	    for (k=1; k<=p; k++) {
		    dk--;
		    tpmv_upper_simd4(dk, Tz[k], x + 4*kk);	
		    kk += dk;
		    tpmv_upper_simd4(dk, Tz[k], x + 4*kk);	
		    kk += dk;
	    }	
        }
        else {
            _FLOAT_ y[4*(FMM_P_MAX+1)*(FMM_P_MAX+2)];
	    dk = p+1;
	    kk = p+1;
	    gemv_simd4(dk, dk, Tz[0], dk, x, y);
	    for (k=1; k<=p; k++) {
		    dk--;
		    gemv_simd4(dk, dk, Tz[k], dk, x + 4*kk, y + 4*kk);	
		    kk += dk;
		    gemv_simd4(dk, dk, Tz[k], dk, x + 4*kk, y + 4*kk);	
		    kk += dk;
	    }	
            memcpy(x, y, 4*sizeof(_FLOAT_)*(p+1)*(p+2));
        }
}

extern _FLOAT_ zeros[];

void M2M(FmmvHandle *FMMV, Box *box)
{
        int p = FMMV->pM;
	int len0 = (p+1)*(p+1);
	int len = (p+1)*(p+2);
        int *P_RT = FMMV->P_MRT;
        int *P_riri2rrii = FMMV->P_Mriri2rrii;
	int no_of_down_childs;
	int no_of_up_childs;

        SIMD_ALIGN _FLOAT_ x1[4*(FMM_P_MAX+1)*(FMM_P_MAX+2)];
        SIMD_ALIGN _FLOAT_ x2[4*(FMM_P_MAX+1)*(FMM_P_MAX+2)];
        SIMD_ALIGN _FLOAT_ xx[4*(FMM_P_MAX+1)*(FMM_P_MAX+2)];

	_FLOAT_ *M_SWD=zeros, *M_NWD=zeros, *M_SED=zeros, *M_NED=zeros;
	_FLOAT_ *M_NEU=zeros, *M_SEU=zeros, *M_NWU=zeros, *M_SWU=zeros;

        if (!isSource(box)) return;

	no_of_down_childs = 0;
	if (isSource(box->child[SWD])&&box->child[SWD]->M) {M_SWD = box->child[SWD]->M; ++no_of_down_childs; }
	if (isSource(box->child[NWD])&&box->child[NWD]->M) {M_NWD = box->child[NWD]->M; ++no_of_down_childs; }
	if (isSource(box->child[SED])&&box->child[SED]->M) {M_SED = box->child[SED]->M; ++no_of_down_childs; }
	if (isSource(box->child[NED])&&box->child[NED]->M) {M_NED = box->child[NED]->M; ++no_of_down_childs; }

	no_of_up_childs = 0;
	if (isSource(box->child[NEU])&&box->child[NEU]->M) {M_NEU = box->child[NEU]->M; ++no_of_up_childs; }
	if (isSource(box->child[SEU])&&box->child[SEU]->M) {M_SEU = box->child[SEU]->M; ++no_of_up_childs; }
	if (isSource(box->child[NWU])&&box->child[NWU]->M) {M_NWU = box->child[NWU]->M; ++no_of_up_childs; }
	if (isSource(box->child[SWU])&&box->child[SWU]->M) {M_SWU = box->child[SWU]->M; ++no_of_up_childs; }

        if (no_of_down_childs>0) { 
		vec4_copy_simd4(len, M_SWD, M_NWD, M_SED, M_NED, x1); 
		Rz_simd4(p, x1, x2, 0);
        	perm_simd4(len0, P_riri2rrii, x2, x1);
        	Ry_simd4(p, FMMV->Ry_pi_minus_theta, x1, x2);
        	perm_simd4(len0, P_RT, x2, x1);
        	Tz_M2M_simd4(FMMV, x1);
        	perm_inv_simd4(len0, P_RT, x1, xx);
        	if (no_of_up_childs>0) { 
			vec4_copy_simd4(len, M_NEU, M_SEU, M_NWU, M_SWU, x1); 
        		Rz_simd4(p, x1, x2, 0);
	        	perm_simd4(len0, P_riri2rrii, x2, x1);
	        	Ry_simd4(p, FMMV->Ry_minus_theta, x1, x2);
	        	perm_simd4(len0, P_RT, x2, x1);
	        	Tz_M2M_simd4(FMMV, x1);
	        	perm_inv_simd4(len0, P_RT, x1, x2);
        		Ry_pi_simd4(p, x2);
                	VEC_ADD(4*len, x2, xx, xx);
        	}
        	Ry_simd4(p, FMMV->Ry_minus_pi_minus_theta, xx, x1);
        	perm_inv_simd4(len0, P_riri2rrii, x1, x2);
        	Rz_simd4(p, x2, x1, 1);
        	box->M = vec4_add2_simd4(FMMV, len, x1, box->M);
	}
        else if (no_of_up_childs>0) {
		vec4_copy_simd4(len, M_NEU, M_SEU, M_NWU, M_SWU, x1); 
                Rz_simd4(p, x1, x2, 0);
                perm_simd4(len0, P_riri2rrii, x2, x1);
                Ry_simd4(p, FMMV->Ry_minus_theta, x1, x2);
                perm_simd4(len0, P_RT, x2,x1);
                Tz_M2M_simd4(FMMV, x1);
                perm_inv_simd4(len0, P_RT, x1,x2);
                Ry_simd4(p, FMMV->Ry_theta, x2, x1);
                perm_inv_simd4(len0, P_riri2rrii, x1, x2);
                Rz_simd4(p, x2, x1, 1);
        	box->M = vec4_add2_simd4(FMMV, len, x1, box->M);
        }	
}


void L2L(FmmvHandle *FMMV, Box *box)
{
	int p = FMMV->pL;
	int len0 = (p+1)*(p+1);
	int len = (p+1)*(p+2);
        int *P_RT = FMMV->P_LRT;
        int *P_riri2rrii = FMMV->P_Lriri2rrii;
        int no_of_down_childs;
        int no_of_up_childs;

        SIMD_ALIGN _FLOAT_ x1[4*(FMM_P_MAX+1)*(FMM_P_MAX+2)];
        SIMD_ALIGN _FLOAT_ x2[4*(FMM_P_MAX+1)*(FMM_P_MAX+2)]; 
        SIMD_ALIGN _FLOAT_ xx[4*(FMM_P_MAX+1)*(FMM_P_MAX+2)];

        _FLOAT_ **L_SWD_p=0, **L_NWD_p=0, **L_SED_p=0, **L_NED_p=0;
        _FLOAT_ **L_NEU_p=0, **L_SEU_p=0, **L_NWU_p=0, **L_SWU_p=0;


        if (!isTarget(box)||!box->L) return;

        no_of_down_childs = 0;
        if (isTarget(box->child[SWD])) {L_SWD_p = &(box->child[SWD]->L); ++no_of_down_childs; }
        if (isTarget(box->child[NWD])) {L_NWD_p = &(box->child[NWD]->L); ++no_of_down_childs; }
        if (isTarget(box->child[SED])) {L_SED_p = &(box->child[SED]->L); ++no_of_down_childs; }
        if (isTarget(box->child[NED])) {L_NED_p = &(box->child[NED]->L); ++no_of_down_childs; }

        no_of_up_childs = 0;
        if (isTarget(box->child[NEU])) {L_NEU_p = &(box->child[NEU]->L); ++no_of_up_childs; }
        if (isTarget(box->child[SEU])) {L_SEU_p = &(box->child[SEU]->L); ++no_of_up_childs; }
        if (isTarget(box->child[NWU])) {L_NWU_p = &(box->child[NWU]->L); ++no_of_up_childs; }
	if (isTarget(box->child[SWU])) {L_SWU_p = &(box->child[SWU]->L); ++no_of_up_childs; }

        if (no_of_up_childs>0) {
                Rz1_simd4(p, box->L, x2, 0);
                perm_simd4(len0, P_riri2rrii, x2, x1);
                Ry_simd4(p, FMMV->Ry_pi_minus_theta, x1, xx);
                perm_simd4(len0, P_RT, xx, x1);
                Tz_L2L_simd4(FMMV, x1);
                perm_inv_simd4(len0, P_RT, x1, x2);
                Ry_simd4(p, FMMV->Ry_minus_pi_minus_theta, x2, x1);
                perm_inv_simd4(len0, P_riri2rrii, x1, x2);
                Rz_simd4(p, x2, x1, 1);
		VEC_4add2_simd4(FMMV, len, x1, L_NEU_p, L_SEU_p, L_NWU_p, L_SWU_p); 
                if (no_of_down_childs>0) {
                        Ry_pi_simd4(p, xx);
                        perm_simd4(len0, P_RT, xx, x1);
                        Tz_L2L_simd4(FMMV, x1);
                        perm_inv_simd4(len0, P_RT, x1, x2);
                        Ry_simd4(p, FMMV->Ry_theta,x2, x1);
                        perm_inv_simd4(len0, P_riri2rrii, x1, x2);
                        Rz_simd4(p, x2, x1, 1);
			VEC_4add2_simd4(FMMV, len, x1, L_SWD_p, L_NWD_p, L_SED_p, L_NED_p); 
                }
        }
        else if (no_of_down_childs>0) {
                Rz1_simd4(p, box->L, x2, 0);
                perm_simd4(len0, P_riri2rrii, x2, x1);
                Ry_simd4(p, FMMV->Ry_minus_theta,x1, x2);
                perm_simd4(len0, P_RT, x2, x1);
                Tz_L2L_simd4(FMMV, x1);
                perm_inv_simd4(len0, P_RT, x1, x2);
                Ry_simd4(p, FMMV->Ry_theta,x2, x1);
                perm_inv_simd4(len0, P_riri2rrii, x1, x2);
                Rz_simd4(p, x2, x1, 1);
		VEC_4add2_simd4(FMMV, len, x1, L_SWD_p, L_NWD_p, L_SED_p, L_NED_p); 
        }
}

