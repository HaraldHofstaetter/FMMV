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

#include<math.h>
#include"_fmmv.h"

#include"simd.h"
#if (FMM_PRECISION==0)
   #include"simd4s.h"
#else
   #include"simd4d.h"
#endif

#if (defined(STANDARD))
	#define FMM_KIND FMM_STANDARD	
	#define gen_M gen_M_standard
	#define eval_L eval_L_standard
	#define gen_L gen_L_standard
	#define eval_M eval_M_standard
	#define gen_L_eval_M gen_L_eval_M_standard
        #define DIPOLE 0	
        #define GRAD 0	
#elif (defined(GRAD))
	#define FMM_KIND FMM_GRAD	
	#define gen_M gen_M_grad
	#define eval_L eval_L_grad
	#define gen_L gen_L_grad
	#define eval_M eval_M_grad
	#define gen_L_eval_M gen_L_eval_M_grad
        #define DIPOLE 0	
        #define GRAD 1	
#elif (defined(DIPOLE))
	#define FMM_KIND FMM_DIPOLE	
	#define gen_M gen_M_dipole
	#define eval_L eval_L_dipole
	#define gen_L gen_L_dipole
	#define eval_M eval_M_dipole
	#define gen_L_eval_M gen_L_eval_M_dipole
        #define DIPOLE 1	
        #define GRAD 0	
#elif (defined(DIPOLE_GRAD))
	#define FMM_KIND FMM_DIPOLE_GRAD	
	#define gen_M gen_M_dipole_grad
	#define eval_L eval_L_dipole_grad
	#define gen_L gen_L_dipole_grad
	#define eval_M eval_M_dipole_grad
	#define gen_L_eval_M gen_L_eval_M_dipole_grad
        #define DIPOLE 1	
        #define GRAD 1	
#elif (defined(ST_STANDARD))
	#define FMM_KIND FMM_ST_STANDARD	
	#define gen_M gen_M_ST_standard
	#define eval_L eval_L_ST_standard
	#define gen_L gen_L_ST_standard
	#define eval_M eval_M_ST_standard
	#define gen_L_eval_M gen_L_eval_M_ST_standard
        #define DIPOLE 0	
        #define GRAD 0	
#elif (defined(ST_GRAD))
	#define FMM_KIND FMM_ST_GRAD	
	#define gen_M gen_M_ST_grad
	#define eval_L eval_L_ST_grad
	#define gen_L gen_L_ST_grad
	#define eval_M eval_M_ST_grad
	#define gen_L_eval_M gen_L_eval_M_ST_grad
        #define DIPOLE 0	
        #define GRAD 1	
#elif (defined(ST_DIPOLE))
	#define FMM_KIND FMM_ST_DIPOLE	
	#define gen_M gen_M_ST_dipole
	#define eval_L eval_L_ST_dipole
	#define gen_L gen_L_ST_dipole
	#define eval_M eval_M_ST_dipole
	#define gen_L_eval_M gen_L_eval_M_ST_dipole
        #define DIPOLE 1	
        #define GRAD 0	
#elif (defined(ST_DIPOLE_GRAD))
	#define FMM_KIND FMM_ST_DIPOLE_GRAD	
	#define gen_M gen_M_ST_dipole_grad
	#define eval_L eval_L_ST_dipole_grad
	#define gen_L gen_L_ST_dipole_grad
	#define eval_M eval_M_ST_dipole_grad
	#define gen_L_eval_M gen_L_eval_M_ST_dipole_grad
        #define DIPOLE 1	
        #define GRAD 1	
#endif

#include "fmmv_access.h"

void gen_M_base_simd4(FmmvHandle *FMMV, Box *box, V4_TYPE x, V4_TYPE y, V4_TYPE z, V4_TYPE q, V4_TYPE mx, V4_TYPE my, V4_TYPE mz, int dipole);
void gen_L_base_simd4(FmmvHandle *FMMV, Box *target, V4_TYPE x, V4_TYPE y, V4_TYPE z, V4_TYPE q, V4_TYPE mx, V4_TYPE my, V4_TYPE mz, int dipole);
void eval_L_base_simd4(FmmvHandle *FMMV, Box *box, V4_TYPE x, V4_TYPE y, V4_TYPE z, V4_TYPE *pot, V4_TYPE *dx, V4_TYPE *dy, V4_TYPE *dz, int grad);
void eval_M_base_simd4(FmmvHandle *FMMV, Box *source, V4_TYPE x, V4_TYPE y, V4_TYPE z, V4_TYPE *pot, V4_TYPE *dx, V4_TYPE *dy, V4_TYPE *dz, int grad);
void gen_L_eval_M_base_simd4(FmmvHandle *FMMV, Box *list3, 
                        V4_TYPE x, V4_TYPE y, V4_TYPE z, V4_TYPE q, V4_TYPE mx, V4_TYPE my, V4_TYPE mz,
                        V4_TYPE *pot, V4_TYPE *dx, V4_TYPE *dy, V4_TYPE *dz,
                        int dipole, int grad);

extern void gen_M_base(FmmvHandle *FMMV, Box *box, _FLOAT_ x, _FLOAT_ y, _FLOAT_ z, _FLOAT_ q, _FLOAT_ mx, _FLOAT_ my, _FLOAT_ mz, int dipole);
extern void gen_L_base(FmmvHandle *FMMV, Box *target, _FLOAT_ x, _FLOAT_ y, _FLOAT_ z, _FLOAT_ q, _FLOAT_ mx, _FLOAT_ my, _FLOAT_ mz, int dipole);
extern void eval_L_base(FmmvHandle *FMMV, Box *box,  _FLOAT_ x, _FLOAT_ y, _FLOAT_ z, _FLOAT_ *pot, _FLOAT_ *dx, _FLOAT_ *dy, _FLOAT_ *dz, int grad);
extern void eval_M_base(FmmvHandle *FMMV, Box *source, _FLOAT_ x, _FLOAT_ y, _FLOAT_ z, _FLOAT_ *pot, _FLOAT_ *dx, _FLOAT_ *dy, _FLOAT_ *dz, int grad);
extern void gen_L_eval_M_base(FmmvHandle *FMMV, Box *list3, 
                        _FLOAT_ x, _FLOAT_ y, _FLOAT_ z, 
                        _FLOAT_ q, _FLOAT_ mx, _FLOAT_ my, _FLOAT_ mz, 
                        _FLOAT_ *pot, _FLOAT_ *dx, _FLOAT_ *dy, _FLOAT_ *dz, int dipole, int grad);






void gen_M(FmmvHandle *FMMV, Box *box)
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	V4_TYPE x, y, z, q;
	V4_BASETYPE x_b, y_b, z_b, q_b;
	V4_TYPE bx, by, bz;
	int i, i0, i1, ii, ii0, ii1, iii, m;
	V4_TYPE mx, my, mz;
	V4_BASETYPE mx_b, my_b, mz_b;

	if (!box) return;

	if (!box->M) {
            int p = FMMV->pM;
            int len = (p+1)*(p+2);
	    box->M = (V4_BASETYPE*) FMMV_MALLOC(FMMV, len*sizeof(V4_BASETYPE));
	    memset(box->M, 0, len*sizeof(V4_BASETYPE));
	}	

	bx = V4_SET1(box->x);
	by = V4_SET1(box->y);
	bz = V4_SET1(box->z);

	i0 = box->firstParticle;
	i1 = i0 + box->noOfParticles;
	ii0 = i0>>2;
	if (i0&3) ii0++;	
	ii1 = i1>>2;

	for (ii=ii0; ii<ii1; ii++ ) {
    	    x = V4_SUB(V4_LOAD(access_x(ii)), bx);
	    y = V4_SUB(V4_LOAD(access_y(ii)), by);
	    z = V4_SUB(V4_LOAD(access_z(ii)), bz);
	    q = V4_LOAD(access_q(ii));
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	    mx = V4_LOAD(access_mx(ii));
	    my = V4_LOAD(access_my(ii));
	    mz = V4_LOAD(access_mz(ii));
        #endif
            gen_M_base_simd4(FMMV, box, x, y, z, q, mx, my, mz, DIPOLE);
	}
	
        m = i0&3; 
        ii = ii0-1;
        if (m) {
	    x = V4_SUB(V4_LOAD_LEFT(access_x(ii), m), bx);
	    y = V4_SUB(V4_LOAD_LEFT(access_y(ii), m), by);
	    z = V4_SUB(V4_LOAD_LEFT(access_z(ii), m), bz);
	    q = V4_LOAD_LEFT0(access_q(ii), m);
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	    mx = V4_LOAD_LEFT0(access_mx(ii), m);
	    my = V4_LOAD_LEFT0(access_my(ii), m);
	    mz = V4_LOAD_LEFT0(access_mz(ii), m);
        #endif
            gen_M_base_simd4(FMMV, box, x, y, z, q, mx, my, mz, DIPOLE);

	    if (ii0==ii1) return;
        }

        m = i1&3; 
        ii = ii1;
        if (m) {
	    x = V4_SUB(V4_LOAD_RIGHT(access_x(ii), m), bx);
	    y = V4_SUB(V4_LOAD_RIGHT(access_y(ii), m), by);
	    z = V4_SUB(V4_LOAD_RIGHT(access_z(ii), m), bz);
	    q = V4_LOAD_RIGHT0(access_q(ii), m);
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	    mx = V4_LOAD_RIGHT0(access_mx(ii), m);
	    my = V4_LOAD_RIGHT0(access_my(ii), m);
	    mz = V4_LOAD_RIGHT0(access_mz(ii), m);
        #endif
            gen_M_base_simd4(FMMV, box, x, y, z, q, mx, my, mz, DIPOLE);
        }

}	


void gen_L(FmmvHandle *FMMV, Box *target, Box *source)
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	V4_TYPE x, y, z, q;
	V4_BASETYPE x_b, y_b, z_b, q_b, m;
	V4_TYPE tx, ty, tz;
	int i, i0, i1, ii, ii0, ii1, iii;
	V4_TYPE mx, my, mz;
	V4_BASETYPE mx_b, my_b, mz_b;
	
	if (!target || !source) return;

	if (!target->L) {
            int p = FMMV->pL;
            int len = (p+1)*(p+2);
	    target->L = (V4_BASETYPE*) FMMV_MALLOC(FMMV, len*sizeof(V4_BASETYPE));
	    memset(target->L, 0, len*sizeof(V4_BASETYPE));
	}	

	tx = V4_SET1(target->x);
	ty = V4_SET1(target->y);
	tz = V4_SET1(target->z);

	i0 = source->firstParticle;
	i1 = i0 + source->noOfParticles;
	ii0 = i0>>2;
	if (i0&3) ii0++;	
	ii1 = i1>>2;

	for (ii=ii0; ii<ii1; ii++ ) {
		x = V4_SUB(V4_LOAD(access_x(ii)), tx);
		y = V4_SUB(V4_LOAD(access_y(ii)), ty);
		z = V4_SUB(V4_LOAD(access_z(ii)), tz);
		q = V4_LOAD(access_q(ii));
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		mx = V4_LOAD(access_mx(ii));
		my = V4_LOAD(access_my(ii));
		mz = V4_LOAD(access_mz(ii));
        #endif
                gen_L_base_simd4(FMMV, target, x, y, z, q, mx, my, mz, DIPOLE);
	}

        m = i0&3; 
        ii = ii0-1;
        if (m) {
		x = V4_SUB(V4_LOAD_LEFT(access_x(ii), m), tx);
		y = V4_SUB(V4_LOAD_LEFT(access_y(ii), m), ty);
		z = V4_SUB(V4_LOAD_LEFT(access_z(ii), m), tz);
		q = V4_LOAD_LEFT0(access_q(ii), m);
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		mx = V4_LOAD_LEFT0(access_mx(ii), m);
		my = V4_LOAD_LEFT0(access_my(ii), m);
		mz = V4_LOAD_LEFT0(access_mz(ii), m);
        #endif
                gen_L_base_simd4(FMMV, target, x, y, z, q, mx, my, mz, DIPOLE);

    	        if (ii0==ii1) return;
        }

        m = i1&3; 
        ii = ii1;
        if (m) {
		x = V4_SUB(V4_LOAD_RIGHT(access_x(ii), m), tx);
		y = V4_SUB(V4_LOAD_RIGHT(access_y(ii), m), ty);
		z = V4_SUB(V4_LOAD_RIGHT(access_z(ii), m), tz);
		q = V4_LOAD_RIGHT0(access_q(ii), m);
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		mx = V4_LOAD_RIGHT0(access_mx(ii), m);
		my = V4_LOAD_RIGHT0(access_my(ii), m);
		mz = V4_LOAD_RIGHT0(access_mz(ii), m);
        #endif
                gen_L_base_simd4(FMMV, target, x, y, z, q, mx, my, mz, DIPOLE);

        }
}        
        

void eval_L(FmmvHandle *FMMV, Box *box)
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	V4_TYPE x, y, z, pot;
	V4_BASETYPE x_b, y_b, z_b, pot_b;
	V4_TYPE bx, by, bz;
	int i, i0, i1, ii, ii0, ii1, iii;
	V4_TYPE dx, dy, dz;
	V4_BASETYPE dx_b, dy_b, dz_b;

	if ((!box)||(!box->L)) return;

	bx = V4_SET1(box->x);
	by = V4_SET1(box->y);
	bz = V4_SET1(box->z);

	i0 = box->firstTarget;
	i1 = i0 + box->noOfTargets;
	ii0 = i0>>2;
	if (i0&3) ii0++;	
	ii1 = i1>>2;

	for (ii=ii0; ii<ii1; ii++ ) {
		x = V4_SUB(V4_LOAD(access_tx(ii)), bx);
		y = V4_SUB(V4_LOAD(access_ty(ii)), by);
		z = V4_SUB(V4_LOAD(access_tz(ii)), bz);
                eval_L_base_simd4(FMMV, box, x, y, z, &pot, &dx, &dy, &dz, GRAD);
		V4_STORE(access_pot(ii), V4_ADD(V4_LOAD(access_pot(ii)),pot));
       	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		V4_STORE(access_gradx(ii), V4_ADD(dx, V4_LOAD(access_gradx(ii))));
		V4_STORE(access_grady(ii), V4_ADD(dy, V4_LOAD(access_grady(ii))));
		V4_STORE(access_gradz(ii), V4_ADD(dz, V4_LOAD(access_gradz(ii))));
        #endif
	}	
	
	if (i0&3) {
	    ii0--; 	
	    iii = (ii0==ii1 ? i1&3 : 4);

	    for (i=i0&3; i<iii; i++) {
		x_b = access_tx(ii0)[i] - box->x;
		y_b = access_ty(ii0)[i] - box->y;
		z_b = access_tz(ii0)[i] - box->z;
                eval_L_base(FMMV, box, x_b, y_b, z_b, &pot_b, &dx_b, &dy_b, &dz_b, GRAD);
		access_pot(ii0)[i] += pot_b;
       	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		access_gradx(ii0)[i] += dx_b;
		access_grady(ii0)[i] += dy_b;
		access_gradz(ii0)[i] += dz_b;
        #endif
	    }
	    if (ii0==ii1) return;
	}    

	for (i=0; i<(i1&3); i++) {
		x_b = access_tx(ii1)[i] - box->x;
		y_b = access_ty(ii1)[i] - box->y;
		z_b = access_tz(ii1)[i] - box->z;
                eval_L_base(FMMV, box, x_b, y_b, z_b, &pot_b, &dx_b, &dy_b, &dz_b, GRAD);
		access_pot(ii1)[i] += pot_b;
       	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		access_gradx(ii1)[i] += dx_b;
		access_grady(ii1)[i] += dy_b;
		access_gradz(ii1)[i] += dz_b;
        #endif
	}
}	

void eval_M(FmmvHandle *FMMV, Box *target, Box *source)
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	V4_TYPE x, y, z, pot;
	V4_BASETYPE x_b, y_b, z_b, pot_b;
	V4_TYPE sx, sy, sz;
	int i, i0, i1, ii, ii0, ii1, iii;
	V4_TYPE dx, dy, dz;
	V4_BASETYPE dx_b, dy_b, dz_b;

	if ((!source)||(!source->M)) return;

	sx = V4_SET1(source->x);
	sy = V4_SET1(source->y);
	sz = V4_SET1(source->z);

	i0 = target->firstTarget;
	i1 = i0 + target->noOfTargets;
	ii0 = i0>>2;
	if (i0&3) ii0++;	
	ii1 = i1>>2;

	for (ii=ii0; ii<ii1; ii++ ) {
		x = V4_SUB(V4_LOAD(access_tx(ii)), sx);
		y = V4_SUB(V4_LOAD(access_ty(ii)), sy);
		z = V4_SUB(V4_LOAD(access_tz(ii)), sz);
                eval_M_base_simd4(FMMV, source, x, y, z, &pot, &dx, &dy, &dz, GRAD);
	        V4_STORE(access_pot(ii), V4_ADD(pot, V4_LOAD(access_pot(ii))));
       	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		V4_STORE(access_gradx(ii), V4_ADD(dx, V4_LOAD(access_gradx(ii))));
		V4_STORE(access_grady(ii), V4_ADD(dy, V4_LOAD(access_grady(ii))));
		V4_STORE(access_gradz(ii), V4_ADD(dz, V4_LOAD(access_gradz(ii))));
        #endif
	}	
	
	if (i0&3) {
	    ii0--; 	
	    iii = (ii0==ii1 ? i1&3 : 4);
	    for (i=i0&3; i<iii; i++) {
		x_b = access_tx(ii0)[i] - source->x;
		y_b = access_ty(ii0)[i] - source->y;
		z_b = access_tz(ii0)[i] - source->z;
                eval_M_base(FMMV, source, x_b, y_b, z_b, &pot_b, &dx_b, &dy_b, &dz_b, GRAD);
		access_pot(ii0)[i] += pot_b;
       	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		access_gradx(ii0)[i] += dx_b;
		access_grady(ii0)[i] += dy_b;
		access_gradz(ii0)[i] += dz_b;
        #endif
	    }
	    if (ii0==ii1) return;
	}    

	for (i=0; i<(i1&3); i++) {
		x_b = access_tx(ii1)[i] - source->x;
		y_b = access_ty(ii1)[i] - source->y;
		z_b = access_tz(ii1)[i] - source->z;
                eval_M_base(FMMV, source, x_b, y_b, z_b, &pot_b, &dx_b, &dy_b, &dz_b, GRAD);
		access_pot(ii1)[i] += pot_b;
       	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		access_gradx(ii1)[i] += dx_b;
		access_grady(ii1)[i] += dy_b;
		access_gradz(ii1)[i] += dz_b;
        #endif
	}
}	

#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_DIPOLE)  \
   ||(FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD))

void gen_L_eval_M(FmmvHandle *FMMV, Box *list3, Box *list4)
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	V4_TYPE x, y, z, q, pot;
	V4_BASETYPE x_b, y_b, z_b, q_b, pot_b;
	V4_TYPE tx, ty, tz;
	int i, i0, i1, ii, ii0, ii1, iii;
	V4_TYPE mx, my, mz;
	V4_BASETYPE mx_b, my_b, mz_b;
	V4_TYPE dx, dy, dz;
	V4_BASETYPE dx_b, dy_b, dz_b;
	
	if (!list3 || !list4) return;

	if (!list3->L) {
            int pL = FMMV->pL;
            int len = (pL+1)*(pL+2);
	    list3->L = (V4_BASETYPE*) FMMV_MALLOC(FMMV, len*sizeof(V4_BASETYPE));
	    memset(list3->L, 0, len*sizeof(V4_BASETYPE));
	}	

	tx = V4_SET1(list3->x);
	ty = V4_SET1(list3->y);
	tz = V4_SET1(list3->z);

	i0 = list4->firstParticle;
	i1 = i0 + list4->noOfParticles;
	ii0 = i0>>2;
	if (i0&3) ii0++;	
	ii1 = i1>>2;

	for (ii=ii0; ii<ii1; ii++ ) {
		x = V4_SUB(V4_LOAD(access_x(ii)), tx);
		y = V4_SUB(V4_LOAD(access_y(ii)), ty);
		z = V4_SUB(V4_LOAD(access_z(ii)), tz);
		q = V4_LOAD(access_q(ii));
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
		mx = V4_LOAD(access_mx(ii));
		my = V4_LOAD(access_my(ii));
		mz = V4_LOAD(access_mz(ii));
        #endif
                gen_L_eval_M_base_simd4(FMMV, list3, x, y, z, q, mx, my, mz, 
                                  &pot, &dx, &dy, &dz, DIPOLE, GRAD);
	        V4_STORE(access_pot(ii), V4_ADD(pot, V4_LOAD(access_pot(ii))));
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
		V4_STORE(access_gradx(ii), V4_ADD(dx, V4_LOAD(access_gradx(ii))));
		V4_STORE(access_grady(ii), V4_ADD(dy, V4_LOAD(access_grady(ii))));
		V4_STORE(access_gradz(ii), V4_ADD(dz, V4_LOAD(access_gradz(ii))));
        #endif
	}
	
	if (i0&3) {
	    ii0--; 	
	    iii = (ii0==ii1 ? i1&3 : 4);
	    for (i=i0&3; i<iii; i++) {
		x_b = access_x(ii0)[i] - list3->x;
		y_b = access_y(ii0)[i] - list3->y;
		z_b = access_z(ii0)[i] - list3->z;
		q_b = access_q(ii0)[i];
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
		mx_b = access_mx(ii0)[i];
		my_b = access_my(ii0)[i];
		mz_b = access_mz(ii0)[i];
        #endif
                gen_L_eval_M_base(FMMV, list3, x_b, y_b, z_b, q_b, mx_b, my_b, mz_b, 
                                  &pot_b, &dx_b, &dy_b, &dz_b, DIPOLE, GRAD);
		access_pot(ii0)[i] += pot_b;
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
		access_gradx(ii0)[i] += dx_b;
		access_grady(ii0)[i] += dy_b;
		access_gradz(ii0)[i] += dz_b;
        #endif
	    }
	    if (ii0==ii1) return;
	}    

	for (i=0; i<(i1&3); i++) {
		x_b = access_x(ii1)[i] - list3->x;
		y_b = access_y(ii1)[i] - list3->y;
		z_b = access_z(ii1)[i] - list3->z;
		q_b = access_q(ii1)[i];
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
		mx_b = access_mx(ii1)[i];
		my_b = access_my(ii1)[i];
		mz_b = access_mz(ii1)[i];
        #endif
                gen_L_eval_M_base(FMMV, list3, x_b, y_b, z_b, q_b, mx_b, my_b, mz_b, 
                                  &pot_b, &dx_b, &dy_b, &dz_b, DIPOLE, GRAD);
		access_pot(ii1)[i] += pot_b;
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
		access_gradx(ii1)[i] += dx_b;
		access_grady(ii1)[i] += dy_b;
		access_gradz(ii1)[i] += dz_b;
        #endif
	}
}

#endif



