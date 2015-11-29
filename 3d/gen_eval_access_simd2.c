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

#include "simd.h"
#if (FMM_PRECISION==0)
   #include"simd2s.h"
#else
   #include"simd2d.h"
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

void gen_M_base_simd2(FmmvHandle *FMMV, Box *box, V2_TYPE x, V2_TYPE y, V2_TYPE z, V2_TYPE q, V2_TYPE mx, V2_TYPE my, V2_TYPE mz, int dipole);
void gen_L_base_simd2(FmmvHandle *FMMV, Box *target, V2_TYPE x, V2_TYPE y, V2_TYPE z, V2_TYPE q, V2_TYPE mx, V2_TYPE my, V2_TYPE mz, int dipole);
void eval_L_base_simd2(FmmvHandle *FMMV, Box *box, V2_TYPE x, V2_TYPE y, V2_TYPE z, V2_TYPE *pot, V2_TYPE *dx, V2_TYPE *dy, V2_TYPE *dz, int grad);
void eval_M_base_simd2(FmmvHandle *FMMV, Box *source, V2_TYPE x, V2_TYPE y, V2_TYPE z, V2_TYPE *pot, V2_TYPE *dx, V2_TYPE *dy, V2_TYPE *dz, int grad);
void gen_L_eval_M_base_simd2(FmmvHandle *FMMV, Box *list3, 
                        V2_TYPE x, V2_TYPE y, V2_TYPE z, V2_TYPE q, V2_TYPE mx, V2_TYPE my, V2_TYPE mz,
                        V2_TYPE *pot, V2_TYPE *dx, V2_TYPE *dy, V2_TYPE *dz,
                        int dipole, int grad);

extern void gen_M_base(FmmvHandle *FMMV, Box *box, _FLOAT_ x, _FLOAT_ y, _FLOAT_ z, _FLOAT_ q, _FLOAT_ mx, _FLOAT_ my, _FLOAT_ mz, int dipole);
extern void gen_L_base(FmmvHandle *FMMV, Box *target, _FLOAT_ x, _FLOAT_ y, _FLOAT_ z, _FLOAT_ q, _FLOAT_ mx, _FLOAT_ my, _FLOAT_ mz, int dipole);
extern void eval_L_base(FmmvHandle *FMMV, Box *box,  _FLOAT_ x, _FLOAT_ y, _FLOAT_ z, _FLOAT_ *pot, _FLOAT_ *dx, _FLOAT_ *dy, _FLOAT_ *dz, int grad);
extern void eval_M_base(FmmvHandle *FMMV, Box *source, _FLOAT_ x, _FLOAT_ y, _FLOAT_ z, _FLOAT_ *pot, _FLOAT_ *dx, _FLOAT_ *dy, _FLOAT_ *dz, int grad);
extern void gen_L_eval_M_base(FmmvHandle *FMMV, Box *list3, 
                        _FLOAT_ x, _FLOAT_ y, _FLOAT_ z, 
                        _FLOAT_ q, _FLOAT_ mx, _FLOAT_ my, _FLOAT_ mz, 
                        _FLOAT_ *pot, _FLOAT_ *dx, _FLOAT_ *dy, _FLOAT_ *dz, int dipole, int grad);




void gen_M(FmmvHandle *FMMV, Box * box)
{
    DEFINE_IDA_LOCAL_ALIASES(FMMV)

    V2_TYPE x, y, z, q;
    V2_BASETYPE x_b, y_b, z_b, q_b;
    int i0, i1, ii, ii0, ii1;
    V2_TYPE mx, my, mz;
    V2_TYPE bx, by, bz;
    V2_BASETYPE mx_b, my_b, mz_b;

    if (!box)
	return;

    if (!box->M) {
        int p = FMMV->pM;
        int len = (p+1)*(p+2);
        box->M = (V2_BASETYPE*) FMMV_MALLOC(FMMV, len*sizeof(V2_BASETYPE));
        memset(box->M, 0, len*sizeof(V2_BASETYPE));
    }	

    bx = V2_SET1(box->x);
    by = V2_SET1(box->y);
    bz = V2_SET1(box->z);

    i0 = box->firstParticle;
    i1 = i0 + box->noOfParticles;
    ii0 = i0 >> 1;
    if (i0 & 1)
	ii0++;
    ii1 = i1 >> 1;

    for (ii = ii0; ii < ii1; ii++) {
	x = V2_SUB(V2_LOAD(access_x(ii)), bx);
	y = V2_SUB(V2_LOAD(access_y(ii)), by);
	z = V2_SUB(V2_LOAD(access_z(ii)), bz);
	q = V2_LOAD(access_q(ii));
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	mx = V2_LOAD(access_mx(ii));
	my = V2_LOAD(access_my(ii));
	mz = V2_LOAD(access_mz(ii));
        #endif
        gen_M_base_simd2(FMMV, box, x, y, z, q, mx, my, mz, DIPOLE);
    }
    if (i0 & 1) {
	if (i1 & 1) {
	    x = V2_SUB(V2_LOAD2(access_x(ii0 - 1) + 1, access_x(ii1)), bx);
	    y = V2_SUB(V2_LOAD2(access_y(ii0 - 1) + 1, access_y(ii1)), by);
	    z = V2_SUB(V2_LOAD2(access_z(ii0 - 1) + 1, access_z(ii1)), bz);
	    q = V2_LOAD2(access_q(ii0 - 1) + 1, access_q(ii1));
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
   	    mx = V2_LOAD2(access_mx(ii0 - 1) + 1, access_mx(ii1));
	    my = V2_LOAD2(access_my(ii0 - 1) + 1, access_my(ii1));
            mz = V2_LOAD2(access_mz(ii0 - 1) + 1, access_mz(ii1));
        #endif
            gen_M_base_simd2(FMMV, box, x, y, z, q, mx, my, mz, DIPOLE);
	} else {
	    x_b = access_x(ii0 - 1)[1] - box->x;
	    y_b = access_y(ii0 - 1)[1] - box->y;
	    z_b = access_z(ii0 - 1)[1] - box->z;
	    q_b = access_q(ii0 - 1)[1];
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	    mx_b = access_mx(ii0 - 1)[1];
	    my_b = access_my(ii0 - 1)[1];
	    mz_b = access_mz(ii0 - 1)[1];
        #endif
            gen_M_base(FMMV, box, x_b, y_b, z_b, q_b, mx_b, my_b, mz_b, DIPOLE);
	}
    } else if (i1 & 1) {
	x_b = access_x(ii1)[0] - box->x;
	y_b = access_y(ii1)[0] - box->y;
	z_b = access_z(ii1)[0] - box->z;
	q_b = access_q(ii1)[0];
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	mx_b = access_mx(ii1)[0];
	my_b = access_my(ii1)[0];
	mz_b = access_mz(ii1)[0];
        #endif
        gen_M_base(FMMV, box, x_b, y_b, z_b, q_b, mx_b, my_b, mz_b, DIPOLE);
    }
}


void gen_L(FmmvHandle *FMMV, Box *target, Box *source)
{
    DEFINE_IDA_LOCAL_ALIASES(FMMV)

    V2_TYPE x, y, z, q;
    V2_BASETYPE x_b, y_b, z_b, q_b;
    V2_TYPE tx, ty, tz;
    int i0, i1, ii, ii0, ii1;
    V2_TYPE mx, my, mz;
    V2_BASETYPE mx_b, my_b, mz_b;

    if (!target || !source) return;

    if (!target->L) {
        int p = FMMV->pL;
        int len = (p+1)*(p+2);
        target->L = (V2_BASETYPE*) FMMV_MALLOC(FMMV, len*sizeof(V2_BASETYPE));
        memset(target->L, 0, len*sizeof(V2_BASETYPE));
    }	

    tx = V2_SET1(target->x);
    ty = V2_SET1(target->y);
    tz = V2_SET1(target->z);

    i0 = source->firstParticle;
    i1 = i0 + source->noOfParticles;
    ii0 = i0 >> 1;
    if (i0 & 1)
	ii0++;
    ii1 = i1 >> 1;

    for (ii = ii0; ii < ii1; ii++) {
	x = V2_SUB(V2_LOAD(access_x(ii)), tx);
	y = V2_SUB(V2_LOAD(access_y(ii)), ty);
	z = V2_SUB(V2_LOAD(access_z(ii)), tz);
	q = V2_LOAD(access_q(ii));
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	mx = V2_LOAD(access_mx(ii));
	my = V2_LOAD(access_my(ii));
	mz = V2_LOAD(access_mz(ii));
        #endif
        gen_L_base_simd2(FMMV, target, x, y, z, q, mx, my, mz, DIPOLE);
    }

    if (i0 & 1) {
	if (i1 & 1) {
	    x = V2_SUB(V2_LOAD2(access_x(ii0 - 1) + 1, access_x(ii1)), tx);
	    y = V2_SUB(V2_LOAD2(access_y(ii0 - 1) + 1, access_y(ii1)), ty);
	    z = V2_SUB(V2_LOAD2(access_z(ii0 - 1) + 1, access_z(ii1)), tz);
	    q = V2_LOAD2(access_q(ii0 - 1) + 1, access_q(ii1));
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	    mx = V2_LOAD2(access_mx(ii0 - 1) + 1, access_mx(ii1));
	    my = V2_LOAD2(access_my(ii0 - 1) + 1, access_my(ii1));
	    mz = V2_LOAD2(access_mz(ii0 - 1) + 1, access_mz(ii1));
        #endif
            gen_L_base_simd2(FMMV, target, x, y, z, q, mx, my, mz, DIPOLE);
	} else {
	    x_b = access_x(ii0 - 1)[1] - target->x;
	    y_b = access_y(ii0 - 1)[1] - target->y;
	    z_b = access_z(ii0 - 1)[1] - target->z;
	    q_b = access_q(ii0 - 1)[1];
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	    mx_b = access_mx(ii0 - 1)[1];
	    my_b = access_my(ii0 - 1)[1];
	    mz_b = access_mz(ii0 - 1)[1];
        #endif
            gen_L_base(FMMV, target, x_b, y_b, z_b, q_b, mx_b, my_b, mz_b, DIPOLE);
	}
    } else if (i1 & 1) {
	x_b = access_x(ii1)[0] - target->x;
	y_b = access_y(ii1)[0] - target->y;
	z_b = access_z(ii1)[0] - target->z;
	q_b = access_q(ii1)[0];
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	mx_b = access_mx(ii1)[0];
	my_b = access_my(ii1)[0];
	mz_b = access_mz(ii1)[0];
        #endif
        gen_L_base(FMMV, target, x_b, y_b, z_b, q_b, mx_b, my_b, mz_b, DIPOLE);
    }
}

void eval_L(FmmvHandle *FMMV, Box * box)
{
    DEFINE_IDA_LOCAL_ALIASES(FMMV)

    V2_TYPE x, y, z;
    V2_TYPE pot;
    V2_BASETYPE x_b, y_b, z_b, pot_b;
    V2_TYPE bx, by, bz;
    int i0, i1, ii, ii0, ii1;
    V2_TYPE dx, dy, dz;
    V2_BASETYPE dx_b, dy_b, dz_b;
    #if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
        ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
    V2_TYPE gradx, grady, gradz;
    #endif

    if ((!box)||(!box->L)) return;

    bx = V2_SET1(box->x);
    by = V2_SET1(box->y);
    bz = V2_SET1(box->z);

    i0 = box->firstTarget;
    i1 = i0 + box->noOfTargets;
    ii0 = i0 >> 1;
    if (i0 & 1)
	ii0++;
    ii1 = i1 >> 1;

    for (ii = ii0; ii < ii1; ii++) {
	x = V2_SUB(V2_LOAD(access_tx(ii)), bx);
	y = V2_SUB(V2_LOAD(access_ty(ii)), by);
	z = V2_SUB(V2_LOAD(access_tz(ii)), bz);
        eval_L_base_simd2(FMMV, box, x, y, z, &pot, &dx, &dy, &dz, GRAD);
	V2_STORE(access_pot(ii), V2_ADD(V2_LOAD(access_pot(ii)), pot));
       	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	V2_STORE(access_gradx(ii), V2_ADD(dx, V2_LOAD(access_gradx(ii))));
	V2_STORE(access_grady(ii), V2_ADD(dy, V2_LOAD(access_grady(ii))));
	V2_STORE(access_gradz(ii), V2_ADD(dz, V2_LOAD(access_gradz(ii))));
        #endif
    }
    if (i0 & 1) {
	if (i1 & 1) {
	    x = V2_SUB(V2_LOAD2(access_tx(ii0 - 1) + 1, access_tx(ii1)), bx);
	    y = V2_SUB(V2_LOAD2(access_ty(ii0 - 1) + 1, access_ty(ii1)), by);
	    z = V2_SUB(V2_LOAD2(access_tz(ii0 - 1) + 1, access_tz(ii1)), bz);
            eval_L_base_simd2(FMMV, box, x, y, z, &pot, &dx, &dy, &dz, GRAD);
	    pot = V2_ADD(V2_LOAD2(access_pot(ii0 - 1) + 1, access_pot(ii1)), pot);
	    V2_STORE2(access_pot(ii0 - 1) + 1, access_pot(ii1), pot);
       	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	    gradx = V2_ADD(dx, V2_LOAD2(access_gradx(ii0 - 1) + 1, access_gradx(ii1)));
	    grady = V2_ADD(dy, V2_LOAD2(access_grady(ii0 - 1) + 1, access_gradx(ii1)));
	    gradz = V2_ADD(dz, V2_LOAD2(access_gradz(ii0 - 1) + 1, access_gradz(ii1)));
	    V2_STORE2(access_gradx(ii0 - 1) + 1, access_gradx(ii1), gradx);
	    V2_STORE2(access_grady(ii0 - 1) + 1, access_grady(ii1), grady);
	    V2_STORE2(access_gradz(ii0 - 1) + 1, access_gradz(ii1), gradz);
        #endif
	} else {
	    x_b = access_tx(ii0 - 1)[1] - box->x;
	    y_b = access_ty(ii0 - 1)[1] - box->y;
	    z_b = access_tz(ii0 - 1)[1] - box->z;
            eval_L_base(FMMV, box, x_b, y_b, z_b, &pot_b, &dx_b, &dy_b, &dz_b, GRAD);
	    access_pot(ii0 - 1)[1] += pot_b;
       	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	    access_gradx(ii0 - 1)[1] += dx_b;
	    access_grady(ii0 - 1)[1] += dy_b;
	    access_gradz(ii0 - 1)[1] += dz_b;
        #endif
	}
    } else if (i1 & 1) {
	x_b = access_tx(ii1)[0] - box->x;
	y_b = access_ty(ii1)[0] - box->y;
	z_b = access_tz(ii1)[0] - box->z;
        eval_L_base(FMMV, box, x_b, y_b, z_b, &pot_b, &dx_b, &dy_b, &dz_b, GRAD);
	access_pot(ii1)[0] += pot_b;
       	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	access_gradx(ii1)[0] += dx_b;
	access_grady(ii1)[0] += dy_b;
	access_gradz(ii1)[0] += dz_b;
        #endif
    }
}


void eval_M(FmmvHandle *FMMV, Box * target, Box * source)
{
    DEFINE_IDA_LOCAL_ALIASES(FMMV)

    V2_TYPE x, y, z;
    V2_TYPE pot;
    V2_BASETYPE x_b, y_b, z_b, pot_b;
    V2_TYPE sx, sy, sz;
    int i0, i1, ii, ii0, ii1;
    V2_TYPE dx, dy, dz;
    V2_BASETYPE dx_b, dy_b, dz_b;
    #if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
       ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
    V2_TYPE gradx, grady, gradz;
    #endif

    if ((!source)||(!source->M)) return;

    sx = V2_SET1(source->x);
    sy = V2_SET1(source->y);
    sz = V2_SET1(source->z);

    i0 = target->firstTarget;
    i1 = i0 + target->noOfTargets;
    ii0 = i0 >> 1;
    if (i0 & 1)
	ii0++;
    ii1 = i1 >> 1;

    for (ii = ii0; ii < ii1; ii++) {
	x = V2_SUB(V2_LOAD(access_tx(ii)), sx);
	y = V2_SUB(V2_LOAD(access_ty(ii)), sy);
	z = V2_SUB(V2_LOAD(access_tz(ii)), sz);
        eval_M_base_simd2(FMMV, source, x, y, z, &pot, &dx, &dy, &dz, GRAD);
	V2_STORE(access_pot(ii), V2_ADD(pot, V2_LOAD(access_pot(ii))));
       	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	V2_STORE(access_gradx(ii), V2_ADD(dx, V2_LOAD(access_gradx(ii))));
	V2_STORE(access_grady(ii), V2_ADD(dy, V2_LOAD(access_grady(ii))));
	V2_STORE(access_gradz(ii), V2_ADD(dz, V2_LOAD(access_gradz(ii))));
        #endif
    }
    if (i0 & 1) {
	if (i1 & 1) {
	    x = V2_SUB(V2_LOAD2(access_tx(ii0 - 1) + 1, access_tx(ii1)), sx);
	    y = V2_SUB(V2_LOAD2(access_ty(ii0 - 1) + 1, access_ty(ii1)), sy);
	    z = V2_SUB(V2_LOAD2(access_tz(ii0 - 1) + 1, access_tz(ii1)), sz);
            eval_M_base_simd2(FMMV, source, x, y, z, &pot, &dx, &dy, &dz, GRAD);
	    pot = V2_ADD(pot, V2_LOAD2(access_pot(ii0 - 1) + 1, access_pot(ii1)));
	    V2_STORE2(access_pot(ii0 - 1) + 1, access_pot(ii1), pot);
       	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	    gradx = V2_ADD(dx, V2_LOAD2(access_gradx(ii0 - 1) + 1, access_gradx(ii1)));
	    grady = V2_ADD(dy, V2_LOAD2(access_grady(ii0 - 1) + 1, access_grady(ii1)));
	    gradz = V2_ADD(dz, V2_LOAD2(access_gradz(ii0 - 1) + 1, access_gradz(ii1)));
	    V2_STORE2(access_gradx(ii0 - 1) + 1, access_gradx(ii1), gradx);
	    V2_STORE2(access_grady(ii0 - 1) + 1, access_grady(ii1), grady);
	    V2_STORE2(access_gradz(ii0 - 1) + 1, access_gradz(ii1), gradz);
        #endif
	} else {
	    x_b = access_tx(ii0 - 1)[1] - source->x;
	    y_b = access_ty(ii0 - 1)[1] - source->y;
	    z_b = access_tz(ii0 - 1)[1] - source->z;
            eval_M_base(FMMV, source, x_b, y_b, z_b, &pot_b, &dx_b, &dy_b, &dz_b, GRAD);
	    access_pot(ii0 - 1)[1] += pot_b;
       	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	    access_gradx(ii0 - 1)[1] += dx_b;
	    access_grady(ii0 - 1)[1] += dy_b;
	    access_gradz(ii0 - 1)[1] += dz_b;
        #endif
	}
    } else if (i1 & 1) {
	x_b = access_tx(ii1)[0] - source->x;
	y_b = access_ty(ii1)[0] - source->y;
	z_b = access_tz(ii1)[0] - source->z;
        eval_M_base(FMMV, source, x_b, y_b, z_b, &pot_b, &dx_b, &dy_b, &dz_b, GRAD);
	access_pot(ii1)[0] += pot_b;
       	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	access_gradx(ii1)[0] += dx_b;
	access_grady(ii1)[0] += dy_b;
	access_gradz(ii1)[0] += dz_b;
        #endif
    }
}

#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_DIPOLE)  \
   ||(FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD))


void gen_L_eval_M(FmmvHandle *FMMV, Box * list3, Box * list4)
{
    DEFINE_IDA_LOCAL_ALIASES(FMMV)

    V2_TYPE x, y, z, q, pot;
    V2_BASETYPE x_b, y_b, z_b, q_b, pot_b;
    V2_TYPE tx, ty, tz;
    int i0, i1, ii, ii0, ii1;
    V2_TYPE mx, my, mz;
    V2_BASETYPE mx_b, my_b, mz_b;
    V2_TYPE dx, dy, dz;
    V2_BASETYPE dx_b, dy_b, dz_b;
    #if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD))     
    V2_TYPE gradx, grady, gradz;
    #endif

    if (!list3 || !list4)
	return;

    if (!list3->L) {
        int pL = FMMV->pL;
        int len = (pL+1)*(pL+2);
        list3->L = (V2_BASETYPE*) FMMV_MALLOC(FMMV, len*sizeof(V2_BASETYPE));
        memset(list3->L, 0, len*sizeof(V2_BASETYPE));
    }	

    tx = V2_SET1(list3->x);
    ty = V2_SET1(list3->y);
    tz = V2_SET1(list3->z);

    i0 = list4->firstParticle;
    i1 = i0 + list4->noOfParticles;
    ii0 = i0 >> 1;
    if (i0 & 1)
	ii0++;
    ii1 = i1 >> 1;

    for (ii = ii0; ii < ii1; ii++) {
	x = V2_SUB(V2_LOAD(access_x(ii)), tx);
	y = V2_SUB(V2_LOAD(access_y(ii)), ty);
	z = V2_SUB(V2_LOAD(access_z(ii)), tz);
	q = V2_LOAD(access_q(ii));
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
	mx = V2_LOAD(access_mx(ii));
	my = V2_LOAD(access_my(ii));
	mz = V2_LOAD(access_mz(ii));
        #endif
        gen_L_eval_M_base_simd2(FMMV, list3, x, y, z, q, mx, my, mz, 
                                &pot, &dx, &dy, &dz, DIPOLE, GRAD);
	V2_STORE(access_pot(ii), V2_ADD(pot, V2_LOAD(access_pot(ii))));
        #if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD))     
	V2_STORE(access_gradx(ii), V2_ADD(dx, V2_LOAD(access_gradx(ii))));
	V2_STORE(access_grady(ii), V2_ADD(dy, V2_LOAD(access_grady(ii))));
	V2_STORE(access_gradz(ii), V2_ADD(dz, V2_LOAD(access_gradz(ii))));
        #endif
    }

    if (i0 & 1) {
	if (i1 & 1) {
	    x = V2_SUB(V2_LOAD2(access_x(ii0 - 1) + 1, access_x(ii1)), tx);
	    y = V2_SUB(V2_LOAD2(access_y(ii0 - 1) + 1, access_y(ii1)), ty);
	    z = V2_SUB(V2_LOAD2(access_z(ii0 - 1) + 1, access_z(ii1)), tz);
   	    q = V2_LOAD2(access_q(ii0 - 1) + 1, access_q(ii1));
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
	    mx = V2_LOAD2(access_mx(ii0 - 1) + 1, access_mx(ii1));
	    my = V2_LOAD2(access_my(ii0 - 1) + 1, access_my(ii1));
	    mz = V2_LOAD2(access_mz(ii0 - 1) + 1, access_mz(ii1));
        #endif
            gen_L_eval_M_base_simd2(FMMV, list3, x, y, z, q, mx, my, mz, 
                                &pot, &dx, &dy, &dz, DIPOLE, GRAD);
	    pot = V2_ADD(pot, V2_LOAD2(access_pot(ii0 - 1) + 1, access_pot(ii1)));
	    V2_STORE2(access_pot(ii0 - 1) + 1, access_pot(ii1), pot);
        #if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD))     
	    gradx = V2_ADD(dx, V2_LOAD2(access_gradx(ii0 - 1) + 1, access_gradx(ii1)));
	    grady = V2_ADD(dy, V2_LOAD2(access_grady(ii0 - 1) + 1, access_grady(ii1)));
	    gradz = V2_ADD(dz, V2_LOAD2(access_gradz(ii0 - 1) + 1, access_gradz(ii1)));
	    V2_STORE2(access_gradx(ii0 - 1) + 1, access_gradx(ii1), gradx);
	    V2_STORE2(access_grady(ii0 - 1) + 1, access_grady(ii1), grady);
	    V2_STORE2(access_gradz(ii0 - 1) + 1, access_gradz(ii1), gradz);
        #endif
	} else {
	    x_b = access_x(ii0 - 1)[1] - list3->x;
	    y_b = access_y(ii0 - 1)[1] - list3->y;
	    z_b = access_z(ii0 - 1)[1] - list3->z;
	    q_b = access_q(ii0 - 1)[1];
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
	    mx_b = access_mx(ii0 - 1)[1];
	    my_b = access_my(ii0 - 1)[1];
	    mz_b = access_mz(ii0 - 1)[1];
        #endif
            gen_L_eval_M_base(FMMV, list3, x_b, y_b, z_b, q_b, mx_b, my_b, mz_b, 
                              &pot_b, &dx_b, &dy_b, &dz_b, DIPOLE, GRAD);
	    access_pot(ii0 - 1)[1] += pot_b;
        #if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD))     
	    access_gradx(ii0 - 1)[1] += dx_b;
	    access_grady(ii0 - 1)[1] += dy_b;
	    access_gradz(ii0 - 1)[1] += dz_b;
        #endif
	}
    } else if (i1 & 1) {
	x_b = access_x(ii1)[0] - list3->x;
	y_b = access_y(ii1)[0] - list3->y;
	z_b = access_z(ii1)[0] - list3->z;
	q_b = access_q(ii1)[0];
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
	mx_b = access_mx(ii1)[0];
	my_b = access_my(ii1)[0];
	mz_b = access_mz(ii1)[0];
        #endif
        gen_L_eval_M_base(FMMV, list3, x_b, y_b, z_b, q_b, mx_b, my_b, mz_b, 
                          &pot_b, &dx_b, &dy_b, &dz_b, DIPOLE, GRAD);
	access_pot(ii1)[0] += pot_b;
        #if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD))     
	access_gradx(ii1)[0] += dx_b;
	access_grady(ii1)[0] += dy_b;
	access_gradz(ii1)[0] += dz_b;
        #endif
    }
}

#endif
