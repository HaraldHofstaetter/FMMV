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

void gen_M_base_simd4(FmmvHandle *FMMV, Box *box, V4_TYPE x, V4_TYPE y, V4_TYPE q, V4_TYPE mx, V4_TYPE my, int dipole);
void gen_L_base_simd4(FmmvHandle *FMMV, Box *target, V4_TYPE x, V4_TYPE y, V4_TYPE q, V4_TYPE mx, V4_TYPE my, int dipole);
void eval_L_base_simd4(FmmvHandle *FMMV, Box *box, V4_TYPE x, V4_TYPE y, V4_TYPE *pot, V4_TYPE *dx, V4_TYPE *dy, int grad);
void eval_M_base_simd4(FmmvHandle *FMMV, Box *source, V4_TYPE x, V4_TYPE y, V4_TYPE *pot, V4_TYPE *dx, V4_TYPE *dy, int grad);
void gen_L_eval_M_base_simd4(FmmvHandle *FMMV, Box *list3, 
                        V4_TYPE x, V4_TYPE y, V4_TYPE q, V4_TYPE mx, V4_TYPE my,
                        V4_TYPE *pot, V4_TYPE *dx, V4_TYPE *dy,
                        int dipole, int grad);

extern void gen_M_base(FmmvHandle *FMMV, Box *box, _FLOAT_ x, _FLOAT_ y, _FLOAT_ q, _FLOAT_ mx, _FLOAT_ my,  int dipole);
extern void gen_L_base(FmmvHandle *FMMV, Box *target, _FLOAT_ x, _FLOAT_ y,  _FLOAT_ q, _FLOAT_ mx, _FLOAT_ my, int dipole);
extern void eval_L_base(FmmvHandle *FMMV, Box *box,  _FLOAT_ x, _FLOAT_ y,  _FLOAT_ *pot, _FLOAT_ *dx, _FLOAT_ *dy, int grad);
extern void eval_M_base(FmmvHandle *FMMV, Box *source, _FLOAT_ x, _FLOAT_ y,  _FLOAT_ *pot, _FLOAT_ *dx, _FLOAT_ *dy, int grad);
extern void gen_L_eval_M_base(FmmvHandle *FMMV, Box *list3, 
                        _FLOAT_ x, _FLOAT_ y, 
                        _FLOAT_ q, _FLOAT_ mx, _FLOAT_ my,
                        _FLOAT_ *pot, _FLOAT_ *dx, _FLOAT_ *dy, int dipole, int grad);


void gen_M(FmmvHandle *FMMV, Box *box)
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	V4_TYPE x, y, q;
	V4_BASETYPE x_b, y_b, q_b;
	V4_TYPE bx, by;
	int i, i0, i1, ii, ii0, ii1, iii;
	V4_TYPE mx, my;
	V4_BASETYPE mx_b, my_b;

	if (!box) return;

	if (!box->M) {
            int p = FMMV->pM;
            int len = (p+1)*(p+2);
	    box->M = (V4_BASETYPE*) FMMV_MALLOC(FMMV, len*sizeof(V4_BASETYPE));
	    memset(box->M, 0, len*sizeof(V4_BASETYPE));
	}	

	bx = V4_SET1(box->x);
	by = V4_SET1(box->y);

	i0 = box->firstParticle;
	i1 = i0 + box->noOfParticles;
	ii0 = i0>>2;
	if (i0&3) ii0++;	
	ii1 = i1>>2;

	for (ii=ii0; ii<ii1; ii++ ) {
		x = V4_SUB(V4_LOAD(access_x(ii)), bx);
		y = V4_SUB(V4_LOAD(access_y(ii)), by);
		q = V4_LOAD(access_q(ii));
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		mx = V4_LOAD(access_mx(ii));
		my = V4_LOAD(access_my(ii));
        #endif
                gen_M_base_simd4(FMMV, box, x, y, q, mx, my, DIPOLE);
	}
	
	if (i0&3) {
	    ii0--; 	
	    iii = (ii0==ii1 ? i1&3 : 4);
	    for (i=i0&3; i<iii; i++) {
		x_b = access_x(ii0)[i] - box->x;
		y_b = access_y(ii0)[i] - box->y;
		q_b = access_q(ii0)[i];
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		mx_b = access_mx(ii0)[i];
		my_b = access_my(ii0)[i];
        #endif
                gen_M_base(FMMV, box, x_b, y_b, q_b, mx_b, my_b, DIPOLE);
	     }
	     if (ii0==ii1) return;
	}

	for (i=0; i<(i1&3); i++) {
		x_b = access_x(ii1)[i] - box->x;
		y_b = access_y(ii1)[i] - box->y;
		q_b = access_q(ii1)[i];
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		mx_b = access_mx(ii1)[i];
		my_b = access_my(ii1)[i];
        #endif
                gen_M_base(FMMV, box, x_b, y_b, q_b, mx_b, my_b, DIPOLE);
	}
}	


void gen_L(FmmvHandle *FMMV, Box *target, Box *source)
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	V4_TYPE x, y, q;
	V4_BASETYPE x_b, y_b, q_b;
	V4_TYPE tx, ty;
	int i, i0, i1, ii, ii0, ii1, iii;
	V4_TYPE mx, my;
	V4_BASETYPE mx_b, my_b;
	
	if (!target || !source) return;

	if (!target->L) {
            int p = FMMV->pL;
            int len = (p+1)*(p+2);
	    target->L = (V4_BASETYPE*) FMMV_MALLOC(FMMV, len*sizeof(V4_BASETYPE));
	    memset(target->L, 0, len*sizeof(V4_BASETYPE));
	}	

	tx = V4_SET1(target->x);
	ty = V4_SET1(target->y);

	i0 = source->firstParticle;
	i1 = i0 + source->noOfParticles;
	ii0 = i0>>2;
	if (i0&3) ii0++;	
	ii1 = i1>>2;

	for (ii=ii0; ii<ii1; ii++ ) {
		x = V4_SUB(V4_LOAD(access_x(ii)), tx);
		y = V4_SUB(V4_LOAD(access_y(ii)), ty);
		q = V4_LOAD(access_q(ii));
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		mx = V4_LOAD(access_mx(ii));
		my = V4_LOAD(access_my(ii));
        #endif
                gen_L_base_simd4(FMMV, target, x, y, q, mx, my, DIPOLE);
	}
	
	if (i0&3) {
	    ii0--; 	
	    iii = (ii0==ii1 ? i1&3 : 4);
	    for (i=i0&3; i<iii; i++) {
		x_b = access_x(ii0)[i] - target->x;
		y_b = access_y(ii0)[i] - target->y;
		q_b = access_q(ii0)[i];
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		mx_b = access_mx(ii0)[i];
		my_b = access_my(ii0)[i];
        #endif
                gen_L_base(FMMV, target, x_b, y_b, q_b, mx_b, my_b, DIPOLE);
	    }
	    if (ii0==ii1) return;
	}    

	for (i=0; i<(i1&3); i++) {
		x_b = access_x(ii1)[i] - target->x;
		y_b = access_y(ii1)[i] - target->y;
		q_b = access_q(ii1)[i];
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		mx_b = access_mx(ii1)[i];
		my_b = access_my(ii1)[i];
        #endif
                gen_L_base(FMMV, target, x_b, y_b, q_b, mx_b, my_b, DIPOLE);
	}
}	

void eval_L(FmmvHandle *FMMV, Box *box)
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	V4_TYPE x, y, pot;
	V4_BASETYPE x_b, y_b, pot_b;
	V4_TYPE bx, by;
	int i, i0, i1, ii, ii0, ii1, iii;
	V4_TYPE dx, dy;
	V4_BASETYPE dx_b, dy_b;

	if ((!box)||(!box->L)) return;

	bx = V4_SET1(box->x);
	by = V4_SET1(box->y);

	i0 = box->firstTarget;
	i1 = i0 + box->noOfTargets;
	ii0 = i0>>2;
	if (i0&3) ii0++;	
	ii1 = i1>>2;

	for (ii=ii0; ii<ii1; ii++ ) {
		x = V4_SUB(V4_LOAD(access_tx(ii)), bx);
		y = V4_SUB(V4_LOAD(access_ty(ii)), by);
                eval_L_base_simd4(FMMV, box, x, y, &pot, &dx, &dy, GRAD);
		V4_STORE(access_pot(ii), V4_ADD(V4_LOAD(access_pot(ii)),pot));
       	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		V4_STORE(access_gradx(ii), V4_ADD(dx, V4_LOAD(access_gradx(ii))));
		V4_STORE(access_grady(ii), V4_ADD(dy, V4_LOAD(access_grady(ii))));
        #endif
	}	
	
	if (i0&3) {
	    ii0--; 	
	    iii = (ii0==ii1 ? i1&3 : 4);

	    for (i=i0&3; i<iii; i++) {
		x_b = access_tx(ii0)[i] - box->x;
		y_b = access_ty(ii0)[i] - box->y;
                eval_L_base(FMMV, box, x_b, y_b, &pot_b, &dx_b, &dy_b, GRAD);
		access_pot(ii0)[i] += pot_b;
       	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		access_gradx(ii0)[i] += dx_b;
		access_grady(ii0)[i] += dy_b;
        #endif
	    }
	    if (ii0==ii1) return;
	}    

	for (i=0; i<(i1&3); i++) {
		x_b = access_tx(ii1)[i] - box->x;
		y_b = access_ty(ii1)[i] - box->y;
                eval_L_base(FMMV, box, x_b, y_b, &pot_b, &dx_b, &dy_b, GRAD);
		access_pot(ii1)[i] += pot_b;
       	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		access_gradx(ii1)[i] += dx_b;
		access_grady(ii1)[i] += dy_b;
        #endif
	}
}	

void eval_M(FmmvHandle *FMMV, Box *target, Box *source)
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	V4_TYPE x, y, pot;
	V4_BASETYPE x_b, y_b, pot_b;
	V4_TYPE sx, sy;
	int i, i0, i1, ii, ii0, ii1, iii;
	V4_TYPE dx, dy;
	V4_BASETYPE dx_b, dy_b;

	if ((!source)||(!source->M)) return;

	sx = V4_SET1(source->x);
	sy = V4_SET1(source->y);

	i0 = target->firstTarget;
	i1 = i0 + target->noOfTargets;
	ii0 = i0>>2;
	if (i0&3) ii0++;	
	ii1 = i1>>2;

	for (ii=ii0; ii<ii1; ii++ ) {
		x = V4_SUB(V4_LOAD(access_tx(ii)), sx);
		y = V4_SUB(V4_LOAD(access_ty(ii)), sy);
                eval_M_base_simd4(FMMV, source, x, y, &pot, &dx, &dy, GRAD);
	        V4_STORE(access_pot(ii), V4_ADD(pot, V4_LOAD(access_pot(ii))));
       	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		V4_STORE(access_gradx(ii), V4_ADD(dx, V4_LOAD(access_gradx(ii))));
		V4_STORE(access_grady(ii), V4_ADD(dy, V4_LOAD(access_grady(ii))));
        #endif
	}	
	
	if (i0&3) {
	    ii0--; 	
	    iii = (ii0==ii1 ? i1&3 : 4);
	    for (i=i0&3; i<iii; i++) {
		x_b = access_tx(ii0)[i] - source->x;
		y_b = access_ty(ii0)[i] - source->y;
                eval_M_base(FMMV, source, x_b, y_b, &pot_b, &dx_b, &dy_b, GRAD);
		access_pot(ii0)[i] += pot_b;
       	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		access_gradx(ii0)[i] += dx_b;
		access_grady(ii0)[i] += dy_b;
        #endif
	    }
	    if (ii0==ii1) return;
	}    

	for (i=0; i<(i1&3); i++) {
		x_b = access_tx(ii1)[i] - source->x;
		y_b = access_ty(ii1)[i] - source->y;
                eval_M_base(FMMV, source, x_b, y_b, &pot_b, &dx_b, &dy_b, GRAD);
		access_pot(ii1)[i] += pot_b;
       	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		access_gradx(ii1)[i] += dx_b;
		access_grady(ii1)[i] += dy_b;
        #endif
	}
}	

#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_DIPOLE)  \
   ||(FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD))

void gen_L_eval_M(FmmvHandle *FMMV, Box *list3, Box *list4)
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	V4_TYPE x, y, q, pot;
	V4_BASETYPE x_b, y_b, q_b, pot_b;
	V4_TYPE tx, ty;
	int i, i0, i1, ii, ii0, ii1, iii;
	V4_TYPE mx, my;
	V4_BASETYPE mx_b, my_b;
	V4_TYPE dx, dy;
	V4_BASETYPE dx_b, dy_b;
	
	if (!list3 || !list4) return;

	if (!list3->L) {
            int pL = FMMV->pL;
            int len = (pL+1)*(pL+2);
	    list3->L = (V4_BASETYPE*) FMMV_MALLOC(FMMV, len*sizeof(V4_BASETYPE));
	    memset(list3->L, 0, len*sizeof(V4_BASETYPE));
	}	

	tx = V4_SET1(list3->x);
	ty = V4_SET1(list3->y);

	i0 = list4->firstParticle;
	i1 = i0 + list4->noOfParticles;
	ii0 = i0>>2;
	if (i0&3) ii0++;	
	ii1 = i1>>2;

	for (ii=ii0; ii<ii1; ii++ ) {
		x = V4_SUB(V4_LOAD(access_x(ii)), tx);
		y = V4_SUB(V4_LOAD(access_y(ii)), ty);
		q = V4_LOAD(access_q(ii));
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
		mx = V4_LOAD(access_mx(ii));
		my = V4_LOAD(access_my(ii));
        #endif
                gen_L_eval_M_base_simd4(FMMV, list3, x, y, q, mx, my, 
                                  &pot, &dx, &dy, DIPOLE, GRAD);
	        V4_STORE(access_pot(ii), V4_ADD(pot, V4_LOAD(access_pot(ii))));
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
		V4_STORE(access_gradx(ii), V4_ADD(dx, V4_LOAD(access_gradx(ii))));
		V4_STORE(access_grady(ii), V4_ADD(dy, V4_LOAD(access_grady(ii))));
        #endif
	}
	
	if (i0&3) {
	    ii0--; 	
	    iii = (ii0==ii1 ? i1&3 : 4);
	    for (i=i0&3; i<iii; i++) {
		x_b = access_x(ii0)[i] - list3->x;
		y_b = access_y(ii0)[i] - list3->y;
		q_b = access_q(ii0)[i];
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
		mx_b = access_mx(ii0)[i];
		my_b = access_my(ii0)[i];
        #endif
                gen_L_eval_M_base(FMMV, list3, x_b, y_b, q_b, mx_b, my_b, 
                                  &pot_b, &dx_b, &dy_b, DIPOLE, GRAD);
		access_pot(ii0)[i] += pot_b;
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
		access_gradx(ii0)[i] += dx_b;
		access_grady(ii0)[i] += dy_b;
        #endif
	    }
	    if (ii0==ii1) return;
	}    

	for (i=0; i<(i1&3); i++) {
		x_b = access_x(ii1)[i] - list3->x;
		y_b = access_y(ii1)[i] - list3->y;
		q_b = access_q(ii1)[i];
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
		mx_b = access_mx(ii1)[i];
		my_b = access_my(ii1)[i];
        #endif
                gen_L_eval_M_base(FMMV, list3, x_b, y_b, q_b, mx_b, my_b, 
                                  &pot_b, &dx_b, &dy_b, DIPOLE, GRAD);
		access_pot(ii1)[i] += pot_b;
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
		access_gradx(ii1)[i] += dx_b;
		access_grady(ii1)[i] += dy_b;
        #endif
	}
}

#endif



