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

static void Tz_M2M(FmmvHandle *FMMV, _FLOAT_ *x)
{
	int p = FMMV->pM; 
	_FLOAT_ **Tz = FMMV->Tz_M2M;
	int k, dk, kk;
        
        if (FMMV->beta==0) {
	    dk = p+1;
	    kk = p+1;
	    tpmv_lower(dk, Tz[0], x);
	    for (k=1; k<=p; k++) {
		    dk--;
		    tpmv_lower(dk, Tz[k], x + kk);	
		    kk += dk;
		    tpmv_lower(dk, Tz[k], x + kk);	
		    kk += dk;
	    }	
        }
        else {
            _FLOAT_ y[(FMM_P_MAX+1)*(FMM_P_MAX+2)];
	    dk = p+1;
	    kk = p+1;
	    gemv(dk, dk, Tz[0], dk, x, y);
	    for (k=1; k<=p; k++) {
		    dk--;
		    gemv(dk, dk, Tz[k], dk, x + kk, y + kk);	
		    kk += dk;
		    gemv(dk, dk, Tz[k], dk, x + kk, y + kk);	
		    kk += dk;
	    }	
            memcpy(x, y, sizeof(_FLOAT_)*(p+1)*(p+2));
        }
}	


static void Tz_L2L(FmmvHandle *FMMV, _FLOAT_ *x)
{
	int p = FMMV->pL; 
	_FLOAT_ **Tz = FMMV->Tz_L2L;
	int k, dk, kk;

        if (FMMV->beta==0) {
	    dk = p+1;
	    kk = p+1;
	    tpmv_upper(dk, Tz[0], x);
	    for (k=1; k<=p; k++) {
		    dk--;
		    tpmv_upper(dk, Tz[k], x + kk);	
		    kk += dk;
		    tpmv_upper(dk, Tz[k], x + kk);	
		    kk += dk;
	    }	
        }
        else {
            _FLOAT_ y[(FMM_P_MAX+1)*(FMM_P_MAX+2)];
	    dk = p+1;
	    kk = p+1;
	    gemv(dk, dk, Tz[0], dk, x, y);
	    for (k=1; k<=p; k++) {
		    dk--;
		    gemv(dk, dk, Tz[k], dk, x + kk, y + kk);	
		    kk += dk;
		    gemv(dk, dk, Tz[k], dk, x + kk, y + kk);	
		    kk += dk;
	    }	
            memcpy(x, y, sizeof(_FLOAT_)*(p+1)*(p+2));
        }
}


void M2M(FmmvHandle *FMMV, Box *box)
{
        int p = FMMV->pM;
	int len0 = (p+1)*(p+1);
	int len = (p+1)*(p+2);
        int *P_RT = FMMV->P_MRT;
	int *P_riri2rrii = FMMV->P_Mriri2rrii;

	_FLOAT_ x1[(FMM_P_MAX+1)*(FMM_P_MAX+2)];
	_FLOAT_ x2[(FMM_P_MAX+1)*(FMM_P_MAX+2)];
	_FLOAT_ xx[(FMM_P_MAX+1)*(FMM_P_MAX+2)];

	if (!isSource(box)) return;
	
	if (isSource(box->child[SWD])&&box->child[SWD]->M) {
		Rz_pi4(p, box->child[SWD]->M, x2);
		perm(len0, P_riri2rrii, x2, x1);
		Ry(p, FMMV->Ry_pi_minus_theta, x1, x2);
		perm(len0, P_RT, x2, x1);
		Tz_M2M(FMMV, x1);
		perm_inv(len0, P_RT, x1, xx);
		if (isSource(box->child[NEU])&&box->child[NEU]->M) {
			Rz_pi4(p, box->child[NEU]->M, x2);
			perm(len0, P_riri2rrii, x2, x1);
			Ry(p, FMMV->Ry_minus_theta, x1, x2);
			perm(len0, P_RT, x2, x1);
			Tz_M2M(FMMV, x1);
			perm_inv(len0, P_RT, x1, x2);
			Ry_pi(p, x2);
			VEC_ADD(len, x2, xx, xx);
		}
		Ry(p, FMMV->Ry_minus_pi_minus_theta, xx, x1);
		perm_inv(len0, P_riri2rrii, x1, x2);
		Rz_minus_pi4(p, x2, x1);
		box->M = VEC_ADD2(FMMV, len, x1, box->M);
	}
	else if (isSource(box->child[NEU])&&box->child[NEU]->M) {
		Rz_pi4(p, box->child[NEU]->M, x2);
		perm(len0, P_riri2rrii, x2, x1);
		Ry(p, FMMV->Ry_minus_theta,x1, x2);
		perm(len0, P_RT, x2, x1);
		Tz_M2M(FMMV, x1);
		perm_inv(len0, P_RT, x1, x2);
		Ry(p, FMMV->Ry_theta, x2, x1);
		perm_inv(len0, P_riri2rrii, x1, x2);
		Rz_minus_pi4(p, x2, x1);
		box->M = VEC_ADD2(FMMV, len, x1, box->M);
	}
	if (isSource(box->child[NWD])&&box->child[NWD]->M) {
		Rz_minus_pi4(p, box->child[NWD]->M, x2);
		perm(len0, P_riri2rrii, x2, x1);
		Ry(p, FMMV->Ry_pi_minus_theta, x1, x2);
		perm(len0, P_RT, x2, x1);
		Tz_M2M(FMMV, x1);
		perm_inv(len0, P_RT, x1, xx);
		if (isSource(box->child[SEU])&&box->child[SEU]->M) {
			Rz_minus_pi4(p, box->child[SEU]->M, x2);
			perm(len0, P_riri2rrii, x2, x1);
			Ry(p, FMMV->Ry_minus_theta, x1, x2);
			perm(len0, P_RT, x2, x1);
			Tz_M2M(FMMV, x1);
			perm_inv(len0, P_RT, x1, x2);
			Ry_pi(p, x2);
			VEC_ADD(len, x2, xx, xx);
		}
		Ry(p, FMMV->Ry_minus_pi_minus_theta, xx, x1);
		perm_inv(len0, P_riri2rrii, x1, x2);
		Rz_pi4(p, x2, x1);
		box->M = VEC_ADD2(FMMV, len, x1, box->M);
	}
	else if (isSource(box->child[SEU])&&box->child[SEU]->M) {
		Rz_minus_pi4(p, box->child[SEU]->M, x2);
		perm(len0, P_riri2rrii, x2, x1);
		Ry(p, FMMV->Ry_minus_theta, x1, x2);
		perm(len0, P_RT, x2, x1);
		Tz_M2M(FMMV, x1);
		perm_inv(len0, P_RT, x1, x2);
		Ry(p, FMMV->Ry_theta, x2, x1);
		perm_inv(len0, P_riri2rrii, x1, x2);
		Rz_pi4(p, x2, x1);
		box->M = VEC_ADD2(FMMV, len, x1, box->M);
	}
	if (isSource(box->child[SED])&&box->child[SED]->M) {
		Rz_3pi4(p, box->child[SED]->M, x2);
		perm(len0, P_riri2rrii, x2, x1);
		Ry(p, FMMV->Ry_pi_minus_theta, x1, x2);
		perm(len0, P_RT, x2, x1);
		Tz_M2M(FMMV, x1);
		perm_inv(len0, P_RT, x1, xx);
		if (isSource(box->child[NWU])&&box->child[NWU]->M) {
			Rz_3pi4(p, box->child[NWU]->M, x2);
			perm(len0, P_riri2rrii, x2, x1);
			Ry(p, FMMV->Ry_minus_theta, x1, x2);
			perm(len0, P_RT, x2, x1);
			Tz_M2M(FMMV, x1);
			perm_inv(len0, P_RT, x1, x2);
			Ry_pi(p, x2);
			VEC_ADD(len, x2, xx, xx);
		}
		Ry(p, FMMV->Ry_minus_pi_minus_theta, xx, x1);
		perm_inv(len0, P_riri2rrii, x1, x2);
		Rz_minus_3pi4(p, x2, x1);
		box->M = VEC_ADD2(FMMV, len, x1, box->M);
	}
	else if (isSource(box->child[NWU])&&box->child[NWU]->M) {
		Rz_3pi4(p, box->child[NWU]->M, x2);
		perm(len0, P_riri2rrii, x2, x1);
		Ry(p, FMMV->Ry_minus_theta, x1, x2);
		perm(len0, P_RT, x2, x1);
		Tz_M2M(FMMV, x1);
		perm_inv(len0, P_RT, x1, x2);
		Ry(p, FMMV->Ry_theta, x2, x1);
		perm_inv(len0, P_riri2rrii, x1, x2);
		Rz_minus_3pi4(p, x2, x1);
		box->M = VEC_ADD2(FMMV, len, x1, box->M);
	}
	if (isSource(box->child[NED])&&box->child[NED]->M) {
		Rz_minus_3pi4(p, box->child[NED]->M, x2);
		perm(len0, P_riri2rrii, x2, x1);
		Ry(p, FMMV->Ry_pi_minus_theta, x1, x2);
		perm(len0, P_RT, x2, x1);
		Tz_M2M(FMMV, x1);
		perm_inv(len0, P_RT, x1, xx);
		if (isSource(box->child[SWU])&&box->child[SWU]->M) {
			Rz_minus_3pi4(p, box->child[SWU]->M, x2);
			perm(len0, P_riri2rrii, x2, x1);
			Ry(p, FMMV->Ry_minus_theta, x1, x2);
			perm(len0, P_RT, x2, x1);
			Tz_M2M(FMMV, x1);
			perm_inv(len0, P_RT, x1, x2);
			Ry_pi(p, x2);
			VEC_ADD(len, x2, xx, xx);
		}
		Ry(p, FMMV->Ry_minus_pi_minus_theta, xx, x1);
		perm_inv(len0, P_riri2rrii, x1, x2);
		Rz_3pi4(p, x2, x1);
		box->M = VEC_ADD2(FMMV, len, x1, box->M);
	}
	else if (isSource(box->child[SWU])&&box->child[SWU]->M) {
		Rz_minus_3pi4(p, box->child[SWU]->M, x2);
		perm(len0, P_riri2rrii, x2, x1);
		Ry(p, FMMV->Ry_minus_theta, x1, x2);
		perm(len0, P_RT, x2, x1);
		Tz_M2M(FMMV, x1);
		perm_inv(len0, P_RT, x1, x2);
		Ry(p, FMMV->Ry_theta, x2, x1);
		perm_inv(len0, P_riri2rrii, x1, x2);
		Rz_3pi4(p, x2, x1);
		box->M = VEC_ADD2(FMMV, len, x1, box->M);
	}
}

void L2L(FmmvHandle *FMMV, Box *box)
{
	int p = FMMV->pL;
	int len0 = (p+1)*(p+1);
	int len = (p+1)*(p+2);
        int *P_RT = FMMV->P_LRT;
        int *P_riri2rrii = FMMV->P_Lriri2rrii;

	_FLOAT_ x1[(FMM_P_MAX+1)*(FMM_P_MAX+2)];
	_FLOAT_ x2[(FMM_P_MAX+1)*(FMM_P_MAX+2)];
	_FLOAT_ xx[(FMM_P_MAX+1)*(FMM_P_MAX+2)];

	if (!isTarget(box)||!box->L) return;
	if (isTarget(box->child[NEU])) {
		Rz_pi4(p, box->L, x2);
		perm(len0, P_riri2rrii, x2, x1);
		Ry(p, FMMV->Ry_pi_minus_theta, x1, xx);
		perm(len0, P_RT, xx, x1);
		Tz_L2L(FMMV, x1);
		perm_inv(len0, P_RT, x1, x2);
		Ry(p, FMMV->Ry_minus_pi_minus_theta, x2, x1);
		perm_inv(len0, P_riri2rrii, x1, x2);
		Rz_minus_pi4(p, x2, x1);
		box->child[NEU]->L = VEC_ADD2(FMMV, len, x1, box->child[NEU]->L);
		if (isTarget(box->child[SWD])) {
			Ry_pi(p, xx);
			perm(len0, P_RT, xx, x1);
			Tz_L2L(FMMV, x1);
			perm_inv(len0, P_RT, x1, x2);
			Ry(p, FMMV->Ry_theta, x2, x1);
			perm_inv(len0, P_riri2rrii, x1, x2);
			Rz_minus_pi4(p, x2, x1);
			box->child[SWD]->L = VEC_ADD2(FMMV, len, x1, box->child[SWD]->L);
		}
	}
	else if (isTarget(box->child[SWD])) {
		Rz_pi4(p, box->L, x2);
		perm(len0, P_riri2rrii, x2, x1);
		Ry(p, FMMV->Ry_minus_theta, x1, x2);
		perm(len0, P_RT, x2, x1);
		Tz_L2L(FMMV, x1);
		perm_inv(len0, P_RT, x1, x2);
		Ry(p, FMMV->Ry_theta, x2, x1);
		perm_inv(len0, P_riri2rrii, x1, x2);
		Rz_minus_pi4(p, x2, x1);
		box->child[SWD]->L = VEC_ADD2(FMMV, len, x1, box->child[SWD]->L);
	}
	if (isTarget(box->child[SEU])) {
		Rz_minus_pi4(p, box->L, x2);
		perm(len0, P_riri2rrii, x2, x1);
		Ry(p, FMMV->Ry_pi_minus_theta, x1, xx);
		perm(len0, P_RT, xx, x1);
		Tz_L2L(FMMV, x1);
		perm_inv(len0, P_RT, x1, x2);
		Ry(p, FMMV->Ry_minus_pi_minus_theta, x2, x1);
		perm_inv(len0, P_riri2rrii, x1, x2);
		Rz_pi4(p, x2, x1);
		box->child[SEU]->L = VEC_ADD2(FMMV, len, x1, box->child[SEU]->L);
		if (isTarget(box->child[NWD])) {
			Ry_pi(p, xx);
			perm(len0, P_RT, xx, x1);
			Tz_L2L(FMMV, x1);
			perm_inv(len0, P_RT, x1, x2);
			Ry(p, FMMV->Ry_theta, x2, x1);
			perm_inv(len0, P_riri2rrii, x1, x2);
			Rz_pi4(p, x2, x1);
			box->child[NWD]->L = VEC_ADD2(FMMV, len, x1, box->child[NWD]->L);
		}
	}
	else if (isTarget(box->child[NWD])) {
		Rz_minus_pi4(p, box->L, x2);
		perm(len0, P_riri2rrii, x2, x1);
		Ry(p, FMMV->Ry_minus_theta, x1, x2);
		perm(len0, P_RT, x2, x1);
		Tz_L2L(FMMV, x1);
		perm_inv(len0, P_RT, x1, x2);
		Ry(p, FMMV->Ry_theta, x2, x1);
		perm_inv(len0, P_riri2rrii, x1, x2);
		Rz_pi4(p, x2, x1);
		box->child[NWD]->L = VEC_ADD2(FMMV, len, x1, box->child[NWD]->L);
	}
	if (isTarget(box->child[NWU])) {
		Rz_3pi4(p, box->L, x2);
		perm(len0, P_riri2rrii, x2, x1);
		Ry(p, FMMV->Ry_pi_minus_theta, x1, xx);
		perm(len0, P_RT, xx, x1);
		Tz_L2L(FMMV, x1);
		perm_inv(len0, P_RT, x1, x2);
		Ry(p, FMMV->Ry_minus_pi_minus_theta, x2, x1);
		perm_inv(len0, P_riri2rrii, x1, x2);
		Rz_minus_3pi4(p, x2, x1);
		box->child[NWU]->L = VEC_ADD2(FMMV, len, x1, box->child[NWU]->L);
		if (isTarget(box->child[SED])) {
			Ry_pi(p, xx);
			perm(len0, P_RT, xx, x1);
			Tz_L2L(FMMV, x1);
			perm_inv(len0, P_RT, x1, x2);
			Ry(p, FMMV->Ry_theta, x2, x1);
			perm_inv(len0, P_riri2rrii, x1, x2);
			Rz_minus_3pi4(p, x2, x1);
			box->child[SED]->L = VEC_ADD2(FMMV, len, x1, box->child[SED]->L);
		}
	}
	else if (isTarget(box->child[SED])) {
		Rz_3pi4(p, box->L, x2);
		perm(len0, P_riri2rrii, x2, x1);
		Ry(p, FMMV->Ry_minus_theta, x1, x2);
		perm(len0, P_RT, x2, x1);
		Tz_L2L(FMMV, x1);
		perm_inv(len0, P_RT, x1, x2);
		Ry(p, FMMV->Ry_theta, x2, x1);
		perm_inv(len0, P_riri2rrii, x1, x2);
		Rz_minus_3pi4(p, x2, x1);
		box->child[SED]->L = VEC_ADD2(FMMV, len, x1, box->child[SED]->L);
	}
	if (isTarget(box->child[SWU])) {
		Rz_minus_3pi4(p, box->L, x2);
		perm(len0, P_riri2rrii, x2, x1);
		Ry(p, FMMV->Ry_pi_minus_theta, x1, xx);
		perm(len0, P_RT, xx, x1);
		Tz_L2L(FMMV, x1);
		perm_inv(len0, P_RT, x1, x2);
		Ry(p, FMMV->Ry_minus_pi_minus_theta, x2, x1);
		perm_inv(len0, P_riri2rrii, x1, x2);
		Rz_3pi4(p, x2, x1);
		box->child[SWU]->L = VEC_ADD2(FMMV, len, x1, box->child[SWU]->L);
		if (isTarget(box->child[NED])) {
			Ry_pi(p, xx);
			perm(len0, P_RT, xx, x1);
			Tz_L2L(FMMV, x1);
			perm_inv(len0, P_RT, x1, x2);
			Ry(p, FMMV->Ry_theta, x2, x1);
			perm_inv(len0, P_riri2rrii, x1, x2);
			Rz_3pi4(p, x2, x1);
			box->child[NED]->L = VEC_ADD2(FMMV, len, x1, box->child[NED]->L);
		}
	}
	else if (isTarget(box->child[NED])) {
		Rz_minus_3pi4(p, box->L, x2);
		perm(len0, P_riri2rrii, x2, x1);
		Ry(p, FMMV->Ry_minus_theta, x1, x2);
		perm(len0, P_RT, x2, x1);
		Tz_L2L(FMMV, x1);
		perm_inv(len0, P_RT, x1, x2);
		Ry(p, FMMV->Ry_theta, x2, x1);
		perm_inv(len0, P_riri2rrii, x1, x2);
		Rz_3pi4(p, x2, x1);
		box->child[NED]->L = VEC_ADD2(FMMV, len, x1, box->child[NED]->L);
	}
}

