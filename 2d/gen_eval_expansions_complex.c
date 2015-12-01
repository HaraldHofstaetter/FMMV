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

#include"_fmmv.h"

#if (defined(C_CHARGE))
	#define FMM_KIND FMM_C_CHARGE
	#define gen_M gen_M_c_charge
	#define eval_L eval_L_c_charge
	#define gen_L gen_L_c_charge
	#define eval_M eval_M_c_charge
	#define gen_L_eval_M gen_L_eval_M_c_charge
        #define ST 0
        #define CHARGE 1	
        #define DIPOLE 0	
        #define GRAD 0	
#elif (defined(C_DIPOLE))
	#define FMM_KIND FMM_C_DIPOLE
	#define gen_M gen_M_c_dipole
	#define eval_L eval_L_c_dipole
	#define gen_L gen_L_c_dipole
	#define eval_M eval_M_c_dipole
	#define gen_L_eval_M gen_L_eval_M_c_dipole
        #define ST 0
        #define CHARGE 0	
        #define DIPOLE 1	
        #define GRAD 0	
#elif (defined(C_CHARGE_DIPOLE))
	#define FMM_KIND FMM_C_CHARGE_DIPOLE
	#define gen_M gen_M_c_charge_dipole
	#define eval_L eval_L_c_charge_dipole
	#define gen_L gen_L_c_charge_dipole
	#define eval_M eval_M_c_charge_dipole
	#define gen_L_eval_M gen_L_eval_M_c_charge_dipole
        #define ST 0
        #define CHARGE 1	
        #define DIPOLE 1	
        #define GRAD 0	
#elif (defined(C_CHARGE_GRAD))
	#define FMM_KIND FMM_C_CHARGE_GRAD
	#define gen_M gen_M_c_charge_grad
	#define eval_L eval_L_c_charge_grad
	#define gen_L gen_L_c_charge_grad
	#define eval_M eval_M_c_charge_grad
	#define gen_L_eval_M gen_L_eval_M_c_charge_grad
        #define ST 0
        #define CHARGE 1	
        #define DIPOLE 0	
        #define GRAD 1	
#elif (defined(C_DIPOLE_GRAD))
	#define FMM_KIND FMM_C_DIPOLE_GRAD
	#define gen_M gen_M_c_dipole_grad
	#define eval_L eval_L_c_dipole_grad
	#define gen_L gen_L_c_dipole_grad
	#define eval_M eval_M_c_dipole_grad
	#define gen_L_eval_M gen_L_eval_M_c_dipole_grad
        #define ST 0
        #define CHARGE 0	
        #define DIPOLE 1	
        #define GRAD 1	
#elif (defined(C_CHARGE_DIPOLE_GRAD))
	#define FMM_KIND FMM_C_CHARGE_DIPOLE_GRAD
	#define gen_M gen_M_c_charge_dipole_grad
	#define eval_L eval_L_c_charge_dipole_grad
	#define gen_L gen_L_c_charge_dipole_grad
	#define eval_M eval_M_c_charge_dipole_grad
	#define gen_L_eval_M gen_L_eval_M_c_charge_dipole_grad
        #define ST 0
        #define CHARGE 1	
        #define DIPOLE 1	
        #define GRAD 1	
#elif (defined(ST_C_CHARGE))
	#define FMM_KIND FMM_ST_C_CHARGE
	#define gen_M gen_M_ST_c_charge
	#define eval_L eval_L_ST_c_charge
	#define gen_L gen_L_ST_c_charge
	#define eval_M eval_M_ST_c_charge
	#define gen_L_eval_M gen_L_eval_M_ST_c_charge
        #define ST 1
        #define CHARGE 1	
        #define DIPOLE 0	
        #define GRAD 0	
#elif (defined(ST_C_DIPOLE))
	#define FMM_KIND FMM_ST_C_DIPOLE
	#define gen_M gen_M_ST_c_dipole
	#define eval_L eval_L_ST_c_dipole
	#define gen_L gen_L_ST_c_dipole
	#define eval_M eval_M_ST_c_dipole
	#define gen_L_eval_M gen_L_eval_M_ST_c_dipole
        #define ST 1
        #define CHARGE 0	
        #define DIPOLE 1	
        #define GRAD 0	
#elif (defined(ST_C_CHARGE_DIPOLE))
	#define FMM_KIND FMM_ST_C_CHARGE_DIPOLE
	#define gen_M gen_M_ST_c_charge_dipole
	#define eval_L eval_L_ST_c_charge_dipole
	#define gen_L gen_L_ST_c_charge_dipole
	#define eval_M eval_M_ST_c_charge_dipole
	#define gen_L_eval_M gen_L_eval_M_ST_c_charge_dipole
        #define ST 1
        #define CHARGE 1	
        #define DIPOLE 1	
        #define GRAD 0	
#elif (defined(ST_C_CHARGE_GRAD))
	#define FMM_KIND FMM_ST_C_CHARGE_GRAD
	#define gen_M gen_M_ST_c_charge_grad
	#define eval_L eval_L_ST_c_charge_grad
	#define gen_L gen_L_ST_c_charge_grad
	#define eval_M eval_M_ST_c_charge_grad
	#define gen_L_eval_M gen_L_eval_M_ST_c_charge_grad
        #define ST 1
        #define CHARGE 1	
        #define DIPOLE 0	
        #define GRAD 1	
#elif (defined(ST_C_DIPOLE_GRAD))
	#define FMM_KIND FMM_ST_C_DIPOLE_GRAD
	#define gen_M gen_M_ST_c_dipole_grad
	#define eval_L eval_L_ST_c_dipole_grad
	#define gen_L gen_L_ST_c_dipole_grad
	#define eval_M eval_M_ST_c_dipole_grad
	#define gen_L_eval_M gen_L_eval_M_ST_c_dipole_grad
        #define ST 1
        #define CHARGE 0	
        #define DIPOLE 1	
        #define GRAD 1	
#elif (defined(ST_C_CHARGE_DIPOLE_GRAD))
	#define FMM_KIND FMM_ST_C_CHARGE_DIPOLE_GRAD
	#define gen_M gen_M_ST_c_charge_dipole_grad
	#define eval_L eval_L_ST_c_charge_dipole_grad
	#define gen_L gen_L_ST_c_charge_dipole_grad
	#define eval_M eval_M_ST_c_charge_dipole_grad
	#define gen_L_eval_M gen_L_eval_M_ST_c_charge_dipole_grad
        #define ST 1
        #define CHARGE 1	
        #define DIPOLE 1	
        #define GRAD 1	
#endif

#include"fmmv_access.h"


extern void gen_M_base_complex(FmmvHandle *FMMV, Box *box, _FLOAT_ x, _FLOAT_ y, _FLOAT_ qx, _FLOAT_ qy,
                         _FLOAT_ mx, _FLOAT_ my, int charge,  int dipole);
extern void gen_L_base_complex(FmmvHandle *FMMV, Box *target, _FLOAT_ x, _FLOAT_ y,  _FLOAT_ qx, _FLOAT_ qy,
                        _FLOAT_ mx, _FLOAT_ my, int charge, int dipole);
extern void eval_L_base_complex(FmmvHandle *FMMV, Box *box,  _FLOAT_ x, _FLOAT_ y,  _FLOAT_ *potx, _FLOAT_ *poty,
                        _FLOAT_ *dx, _FLOAT_ *dy, int grad);
extern void eval_M_base_complex(FmmvHandle *FMMV, Box *source, _FLOAT_ x, _FLOAT_ y,  _FLOAT_ *potx, _FLOAT_ *poty,
                        _FLOAT_ *dx, _FLOAT_ *dy, int grad);
extern void gen_L_eval_M_base_complex(FmmvHandle *FMMV, Box *list3, 
                        _FLOAT_ x, _FLOAT_ y, 
                        _FLOAT_ qx, _FLOAT_ qy, _FLOAT_ mx, _FLOAT_ my,
                        _FLOAT_ *potx, _FLOAT_ *poty, _FLOAT_ *dx, _FLOAT_ *dy, int charge, int dipole, int grad);

void gen_M(FmmvHandle *FMMV, Box *box)
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	_FLOAT_ x, y;
	_FLOAT_ qx, qy;
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
#if(CHARGE)        
		qx = access_qx(i);
		qy = access_qy(i);
#endif
#if(DIPOLE)
		mx = access_mx(i);
		my = access_my(i);
#endif
                gen_M_base_complex(FMMV, box, x, y, qx, qy, mx, my, CHARGE, DIPOLE);
	}
}	

void gen_L(FmmvHandle *FMMV, Box *target, Box *source)
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	_FLOAT_ x, y;
	_FLOAT_ qx, qy;
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
#if(CHARGE)        
		qx = access_qx(i);
		qy = access_qy(i);
#endif
#if(DIPOLE)
		mx = access_mx(i);
		my = access_my(i);
#endif
                gen_L_base_complex(FMMV, target, x, y, qx, qy, mx, my, CHARGE, DIPOLE);
	}	
}


void eval_L(FmmvHandle *FMMV, Box *box)
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)
	_FLOAT_ x, y;
        _FLOAT_ potx, poty;
	_FLOAT_ dx, dy;
	int i;

	if ((!box)||(!box->L)) return;
	
	for (i=box->firstTarget; 
	     i<box->firstTarget+box->noOfTargets; i++) {
		x = access_tx(i) - box->x;
		y = access_ty(i) - box->y;
                eval_L_base_complex(FMMV, box,x, y, &potx, &poty, &dx, &dy, GRAD);
		access_potx(i) += potx;
		access_poty(i) += poty;
#if(GRAD)           
		access_gradx(i) += dx;
		access_grady(i) += dy;
#endif
	}	
}

void eval_M(FmmvHandle *FMMV, Box *target, Box *source)
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	_FLOAT_ x, y;
        _FLOAT_ potx, poty;
	_FLOAT_ dx, dy;
	int i;

	if ((!source)||(!source->M)) return;

	for (i=target->firstTarget; 
	     i<target->firstTarget+target->noOfTargets; i++) { 
		x = access_tx(i) - source->x;
		y = access_ty(i) - source->y;
                eval_M_base_complex(FMMV, source, x, y, &potx, &poty, &dx, &dy, GRAD);
		access_potx(i) += potx;
		access_poty(i) += poty;
#if(GRAD)           
		access_gradx(i) += dx;
		access_grady(i) += dy;
#endif
	}	
}


#if(!ST)
void gen_L_eval_M(FmmvHandle *FMMV, Box *list3, Box *list4)
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	_FLOAT_ x, y; 
        _FLOAT_ qx, qy;
        _FLOAT_ potx, poty;
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
#if(CHARGE)        
		qx = access_qx(i);
		qy = access_qy(i);
#endif
#if(DIPOLE)
		mx = access_mx(i);
		my = access_my(i);
#endif
                gen_L_eval_M_base_complex(FMMV, list3, x, y, qx, qy, mx, my, 
                                  &potx, &poty, &dx, &dy, CHARGE, DIPOLE, GRAD);
		access_potx(i) += potx;
		access_poty(i) += poty;
#if(GRAD)           
		access_gradx(i) += dx;
		access_grady(i) += dy;
#endif
	}	
}
#endif

