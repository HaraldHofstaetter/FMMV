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

#include"_fmmv.h"


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

#include"fmmv_access.h"


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

	_FLOAT_ x, y, q;
	_FLOAT_ mx, my;
	int i;
	
	if (!box) return;

	if (!box->M) {
	    int p = FMMV->pM;
            int len = 2*(p+1); /* 2*p complex numbers */
	    box->M = (_FLOAT_*) FMMV_MALLOC(FMMV, len*sizeof(_FLOAT_));
	    memset(box->M, 0, len*sizeof(_FLOAT_));
	}	
        
	for (i=box->firstParticle; 
	     i<box->firstParticle+box->noOfParticles; i++) {
		x = access_x(i) - box->x;
		y = access_y(i) - box->y;
		q = access_q(i);
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		mx = access_mx(i);
		my = access_my(i);
        #endif
                gen_M_base(FMMV, box, x, y, q, mx, my, DIPOLE);
	}
}	

void gen_L(FmmvHandle *FMMV, Box *target, Box *source)
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	_FLOAT_ x, y, q;
	_FLOAT_ mx, my;
	int i;
	
	if (!target || !source) return;

	if (!target->L) {
	    int p = FMMV->pL;
            int len = 2*(p+1); /* 2*p complex numbers */
	    target->L = (_FLOAT_*) FMMV_MALLOC(FMMV, len*sizeof(_FLOAT_));
	    memset(target->L, 0, len*sizeof(_FLOAT_));
	}	
	
	for (i=source->firstParticle; 
	     i<source->firstParticle+source->noOfParticles; i++) {
		x = access_x(i) - target->x;
		y = access_y(i) - target->y;
		q = access_q(i);
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	        mx = access_mx(i);
	        my = access_my(i);
        #endif
                gen_L_base(FMMV, target, x, y, q, mx, my, DIPOLE);
	}	
}


void eval_L(FmmvHandle *FMMV, Box *box)
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)
	_FLOAT_ x, y, pot;
	int i;
	_FLOAT_ dx, dy;

	if ((!box)||(!box->L)) return;
	
	for (i=box->firstTarget; 
	     i<box->firstTarget+box->noOfTargets; i++) {
		x = access_tx(i) - box->x;
		y = access_ty(i) - box->y;
                eval_L_base(FMMV, box,x, y, &pot, &dx, &dy, GRAD);
		access_pot(i) += pot;
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		access_gradx(i) += dx;
		access_grady(i) += dy;
	#endif
	}	
}

void eval_M(FmmvHandle *FMMV, Box *target, Box *source)
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	_FLOAT_ x, y, pot;
	_FLOAT_ dx, dy;
	int i;

	if ((!source)||(!source->M)) return;

	for (i=target->firstTarget; 
	     i<target->firstTarget+target->noOfTargets; i++) { 
		x = access_tx(i) - source->x;
		y = access_ty(i) - source->y;
                eval_M_base(FMMV, source, x, y, &pot, &dx, &dy, GRAD);
		access_pot(i) += pot;
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		access_gradx(i) += dx;
		access_grady(i) += dy;
	#endif
	}	
}


#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
   ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
void gen_L_eval_M(FmmvHandle *FMMV, Box *list3, Box *list4)
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	_FLOAT_ x, y, q, pot;
	_FLOAT_ mx, my;
	_FLOAT_ dx, dy;
	int i;
	
	if (!list3 || !list4) return;

	if (!list3->L) {
	    int p = FMMV->pL;
            int len = 2*(p+1); /* 2*p complex numbers */
	    list3->L = (_FLOAT_*) FMMV_MALLOC(FMMV, len*sizeof(_FLOAT_));
	    memset(list3->L, 0, len*sizeof(_FLOAT_));
	}	

	for (i=list4->firstParticle; 
	     i<list4->firstParticle+list4->noOfParticles; i++) {
		x = access_x(i) - list3->x;
		y = access_y(i) - list3->y;
		q = access_q(i);
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
		mx = access_mx(i);
		my = access_my(i);
        #endif
                gen_L_eval_M_base(FMMV, list3, x, y, q, mx, my, 
                                  &pot, &dx, &dy, DIPOLE, GRAD);
		access_pot(i) += pot;
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
		access_gradx(i) += dx;
		access_grady(i) += dy;
        #endif
	}	
}
#endif

