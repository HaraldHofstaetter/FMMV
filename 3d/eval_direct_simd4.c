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

#include "simd.h"
#if (FMM_PRECISION==0)
   #include"simd4s.h"
#else
   #include"simd4d.h"
#endif
#include"math.h"

#if (defined(STANDARD))
	#define FMM_KIND FMM_STANDARD	
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_standard_acc0
	#define eval_direct_yukawa eval_direct_yukawa_standard_acc0
	#define eval_direct_periodic eval_direct_periodic_standard_acc0
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_standard_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_standard_acc1
	#define eval_direct_yukawa eval_direct_yukawa_standard_acc1
	#define eval_direct_periodic eval_direct_periodic_standard_acc1
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_standard_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_standard_acc2
	#define eval_direct_yukawa eval_direct_yukawa_standard_acc2
	#define eval_direct_periodic eval_direct_periodic_standard_acc2
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_standard_acc2
    #endif	
#elif (defined(GRAD))
	#define FMM_KIND FMM_GRAD	
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_grad_acc0
	#define eval_direct_yukawa eval_direct_yukawa_grad_acc0
	#define eval_direct_periodic eval_direct_periodic_grad_acc0
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_grad_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_grad_acc1
	#define eval_direct_yukawa eval_direct_yukawa_grad_acc1
	#define eval_direct_periodic eval_direct_periodic_grad_acc1
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_grad_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_grad_acc2
	#define eval_direct_yukawa eval_direct_yukawa_grad_acc2
	#define eval_direct_periodic eval_direct_periodic_grad_acc2
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_grad_acc2
    #endif	
#elif (defined(DIPOLE))
	#define FMM_KIND FMM_DIPOLE	
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_dipole_acc0
	#define eval_direct_yukawa eval_direct_yukawa_dipole_acc0
	#define eval_direct_periodic eval_direct_periodic_dipole_acc0
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_dipole_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_dipole_acc1
	#define eval_direct_yukawa eval_direct_yukawa_dipole_acc1
	#define eval_direct_periodic eval_direct_periodic_dipole_acc1
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_dipole_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_dipole_acc2
	#define eval_direct_yukawa eval_direct_yukawa_dipole_acc2
	#define eval_direct_periodic eval_direct_periodic_dipole_acc2
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_dipole_acc2
    #endif	
#elif (defined(DIPOLE_GRAD))
	#define FMM_KIND FMM_DIPOLE_GRAD	
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_dipole_grad_acc0
	#define eval_direct_yukawa eval_direct_yukawa_dipole_grad_acc0
	#define eval_direct_periodic eval_direct_periodic_dipole_grad_acc0
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_dipole_grad_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_dipole_grad_acc1
	#define eval_direct_yukawa eval_direct_yukawa_dipole_grad_acc1
	#define eval_direct_periodic eval_direct_periodic_dipole_grad_acc1
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_dipole_grad_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_dipole_grad_acc2
	#define eval_direct_yukawa eval_direct_yukawa_dipole_grad_acc2
	#define eval_direct_periodic eval_direct_periodic_dipole_grad_acc2
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_dipole_grad_acc2
    #endif	
#elif (defined(ST_STANDARD))
	#define FMM_KIND FMM_ST_STANDARD	
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_ST_standard_acc0
	#define eval_direct_yukawa eval_direct_yukawa_ST_standard_acc0
	#define eval_direct_periodic eval_direct_periodic_ST_standard_acc0
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_ST_standard_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_ST_standard_acc1
	#define eval_direct_yukawa eval_direct_yukawa_ST_standard_acc1
	#define eval_direct_periodic eval_direct_periodic_ST_standard_acc1
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_ST_standard_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_ST_standard_acc2
	#define eval_direct_yukawa eval_direct_yukawa_ST_standard_acc2
	#define eval_direct_periodic eval_direct_periodic_ST_standard_acc2
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_ST_standard_acc2
    #endif	
#elif (defined(ST_GRAD))
	#define FMM_KIND FMM_ST_GRAD	
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_ST_grad_acc0
	#define eval_direct_yukawa eval_direct_yukawa_ST_grad_acc0
	#define eval_direct_periodic eval_direct_periodic_ST_grad_acc0
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_ST_grad_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_ST_grad_acc1
	#define eval_direct_yukawa eval_direct_yukawa_ST_grad_acc1
	#define eval_direct_periodic eval_direct_periodic_ST_grad_acc1
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_ST_grad_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_ST_grad_acc2
	#define eval_direct_yukawa eval_direct_yukawa_ST_grad_acc2
	#define eval_direct_periodic eval_direct_periodic_ST_grad_acc2
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_ST_grad_acc2
    #endif	
#elif (defined(ST_DIPOLE))
	#define FMM_KIND FMM_ST_DIPOLE	
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_ST_dipole_acc0
	#define eval_direct_yukawa eval_direct_yukawa_ST_dipole_acc0
	#define eval_direct_periodic eval_direct_periodic_ST_dipole_acc0
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_ST_dipole_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_ST_dipole_acc1
	#define eval_direct_yukawa eval_direct_yukawa_ST_dipole_acc1
	#define eval_direct_periodic eval_direct_periodic_ST_dipole_acc1
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_ST_dipole_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_ST_dipole_acc2
	#define eval_direct_yukawa eval_direct_yukawa_ST_dipole_acc2
	#define eval_direct_periodic eval_direct_periodic_ST_dipole_acc2
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_ST_dipole_acc2
    #endif	
#elif (defined(ST_DIPOLE_GRAD))
	#define FMM_KIND FMM_ST_DIPOLE_GRAD	
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_ST_dipole_grad_acc0
	#define eval_direct_yukawa eval_direct_yukawa_ST_dipole_grad_acc0
	#define eval_direct_periodic eval_direct_periodic_ST_dipole_grad_acc0
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_ST_dipole_grad_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_ST_dipole_grad_acc1
	#define eval_direct_yukawa eval_direct_yukawa_ST_dipole_grad_acc1
	#define eval_direct_periodic eval_direct_periodic_ST_dipole_grad_acc1
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_ST_dipole_grad_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_ST_dipole_grad_acc2
	#define eval_direct_yukawa eval_direct_yukawa_ST_dipole_grad_acc2
	#define eval_direct_periodic eval_direct_periodic_ST_dipole_grad_acc2
	#define eval_direct_yukawa_periodic eval_direct_yukawa_periodic_ST_dipole_grad_acc2
    #endif	
#endif

#ifdef PERIODIC
void eval_direct_yukawa_periodic(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
#else
void eval_direct_yukawa(FmmvHandle *FMMV, Box *target, Box *source);
#endif


#define EVAL_DIRECT_ACCURACY ACCURACY

#include"fmmv_access.h"


#ifdef PERIODIC
void eval_direct_periodic(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx0, _FLOAT_ dy0, _FLOAT_ dz0)
#else
void eval_direct(FmmvHandle *FMMV, Box *target, Box *source)
#endif
{
    #if ((FMM_KIND==FMM_ST_STANDARD)||(FMM_KIND==FMM_ST_GRAD) \
       ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
    if (!(isTarget(target)&&isSource(source))) return;
    #else
    if (!(target&&source)||(target->firstTarget>source->firstParticle)) return;
    #endif

    if (FMMV->beta!=0){
#ifdef PERIODIC
           eval_direct_yukawa_periodic(FMMV, target, source, dx0, dy0, dz0);
#else
           eval_direct_yukawa(FMMV, target, source);
#endif
           return;
    }
    else {

	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	V4_TYPE x,y,z, xi,yi,zi,qj, one_over_r, poti;
	#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
	   ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
	V4_TYPE qi;
	#endif
	int i,j,ni,nj,i0,j0,j00,j1,k;
	int ii, jj, jj0, jj1, jjj;
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))	
	V4_TYPE one_o_r_3;
	#endif
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	V4_TYPE qj_o_r_3;
	V4_TYPE gradxi, gradyi, gradzi;
	#endif
	#if (FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD)
	V4_TYPE qi_o_r_3;
	#endif	
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))	
	V4_TYPE mxj, myj, mzj;
	V4_TYPE m_times_rj;	
	#endif
	#if (FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)
	V4_TYPE mxi, myi, mzi;
	V4_TYPE m_times_ri;	
	#endif
	#if ((FMM_KIND==FMM_DIPOLE_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	V4_TYPE three = V4_SET1(3.0);
	V4_TYPE one_o_r_5;
	#endif	
	
#ifdef PERIODIC
	V4_TYPE dx, dy, dz;

	dx = V4_SET1(dx0);
	dy = V4_SET1(dy0);
	dz = V4_SET1(dz0);
#endif	

	i0 = target->firstTarget;
	ni = target->noOfTargets;
	j0 = source->firstParticle;
	nj = source->noOfParticles;

        #if ((FMM_KIND==FMM_ST_STANDARD)||(FMM_KIND==FMM_ST_GRAD) \
           ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
        FMMV->noOfDirectInteractions += ni*nj;
	j00 = j0;
        #else
        if (i0==j0) {
                FMMV->noOfDirectInteractions += (ni*(ni-1))/2;
        }
        else { 
                FMMV->noOfDirectInteractions += ni*nj;
        }
        #endif
	j1 = j0+nj;
	jj1 = j1>>2;

	for (i=i0; i<i0+ni; i++) {
		ii = i>>2;
		k = i&3;
		#ifdef PERIODIC
		xi = V4_SUB(V4_LOAD1(access_tx(ii)+k), dx);
		yi = V4_SUB(V4_LOAD1(access_ty(ii)+k), dy);
		zi = V4_SUB(V4_LOAD1(access_tz(ii)+k), dz);
		#else
		xi = V4_LOAD1(access_tx(ii)+k);
		yi = V4_LOAD1(access_ty(ii)+k);
		zi = V4_LOAD1(access_tz(ii)+k);
		#endif
		#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
		   ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
		qi = V4_LOAD1(access_q(ii)+k);
		#endif
		#if (FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)
		mxi = V4_LOAD1(access_mx(ii)+k);
		myi = V4_LOAD1(access_my(ii)+k);
		mzi = V4_LOAD1(access_mz(ii)+k);		
		#endif			
		poti = V4_LOAD0(access_pot(ii)+k);
		#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
		   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		gradxi = V4_LOAD0(access_gradx(ii)+k);
		gradyi = V4_LOAD0(access_grady(ii)+k);
		gradzi = V4_LOAD0(access_gradz(ii)+k);
		#endif

		#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
		   ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
		j00 = (i<j0 ? j0 : i+1);
		#endif
		jj0 = j00>>2;
		if (j00&3) jj0++;
		
		for (jj=jj0; jj<jj1; jj++) {
			x = V4_SUB(xi, V4_LOAD(access_x(jj)));
			y = V4_SUB(yi, V4_LOAD(access_y(jj)));
			z = V4_SUB(zi, V4_LOAD(access_z(jj)));
			one_over_r = 
			#if (EVAL_DIRECT_ACCURACY==0)
			V4_RECIP_SQRT0(
			#elif (EVAL_DIRECT_ACCURACY==1)
			V4_RECIP_SQRT1(
			#elif (EVAL_DIRECT_ACCURACY==2)
			V4_RECIP_SQRT2(
			#endif
			  V4_FMA(x, x, V4_FMA(y, y, V4_MUL(z, z))));
			qj = V4_LOAD(access_q(jj));
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
			  ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			mxj = V4_LOAD(access_mx(jj));			
			myj = V4_LOAD(access_my(jj));			
			mzj = V4_LOAD(access_mz(jj));			
			one_o_r_3 = V4_MUL(V4_MUL(one_over_r, one_over_r), one_over_r);
			m_times_rj = V4_FMA(x, mxj, V4_FMA(y, myj, V4_MUL(z, mzj)));
			poti = V4_ADD(poti, V4_FMA(qj, one_over_r, V4_MUL(m_times_rj,one_o_r_3)));
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
			m_times_ri = V4_FMA(x, mxi, V4_FMA(y, myi, V4_MUL(z, mzi)));
			V4_STORE(access_pot(jj), V4_ADD(V4_LOAD(access_pot(jj)), V4_FNMS(m_times_ri,one_o_r_3, V4_MUL(qi, one_over_r))));
			#endif
			#else
			poti = V4_FMA(qj, one_over_r, poti);
			#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD)) 
			V4_STORE(access_pot(jj), V4_FMA(qi, one_over_r, V4_LOAD(access_pot(jj))));
			#endif
			#endif
			
			#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_ST_GRAD))
			one_o_r_3 = V4_MUL(V4_MUL(one_over_r, one_over_r), one_over_r);
			qj_o_r_3 = V4_MUL(qj, one_o_r_3);
			gradxi = V4_FNMS(x, qj_o_r_3, gradxi);
			gradyi = V4_FNMS(y, qj_o_r_3, gradyi);
			gradzi = V4_FNMS(z, qj_o_r_3, gradzi);
			#endif
			#if (FMM_KIND==FMM_GRAD)
			qi_o_r_3 = V4_MUL(qi, one_o_r_3);
			V4_STORE(access_gradx(jj), V4_FMA(x, qi_o_r_3, V4_LOAD(access_gradx(jj))));
			V4_STORE(access_grady(jj), V4_FMA(y, qi_o_r_3, V4_LOAD(access_grady(jj))));
			V4_STORE(access_gradz(jj), V4_FMA(z, qi_o_r_3, V4_LOAD(access_gradz(jj))));
			#endif
			#if ((FMM_KIND==FMM_DIPOLE_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			one_o_r_5 = V4_MUL(one_o_r_3, V4_MUL(one_over_r, one_over_r));
			qj_o_r_3 = V4_FMA(qj, one_o_r_3, V4_MUL(three, V4_MUL(m_times_rj, one_o_r_5)));
			gradxi = V4_FNMS(x, qj_o_r_3, V4_FMA(mxj, one_o_r_3, gradxi));
			gradyi = V4_FNMS(y, qj_o_r_3, V4_FMA(myj, one_o_r_3, gradyi));
			gradzi = V4_FNMS(z, qj_o_r_3, V4_FMA(mzj, one_o_r_3, gradzi));
			#endif
			#if (FMM_KIND==FMM_DIPOLE_GRAD)
			qi_o_r_3 = V4_FNMS(three, V4_MUL(m_times_ri, one_o_r_5), V4_MUL(qi, one_o_r_3));
			V4_STORE(access_gradx(jj), V4_FMA(x, qi_o_r_3, V4_FMA(mxi, one_o_r_3, V4_LOAD(access_gradx(jj)))));
			V4_STORE(access_grady(jj), V4_FMA(y, qi_o_r_3, V4_FMA(myi, one_o_r_3, V4_LOAD(access_grady(jj)))));
			V4_STORE(access_gradz(jj), V4_FMA(z, qi_o_r_3, V4_FMA(mzi, one_o_r_3, V4_LOAD(access_gradz(jj)))));
			#endif
		}

		V4_STORE0(access_pot(ii)+k, V4_HORIZADD1(poti));
		#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
		   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		V4_STORE0(access_gradx(ii)+k, V4_HORIZADD1(gradxi));
		V4_STORE0(access_grady(ii)+k, V4_HORIZADD1(gradyi));
		V4_STORE0(access_gradz(ii)+k, V4_HORIZADD1(gradzi));
		#endif

	    if ((j1>j00)&&((j00&3)||(j1&3))) {
		V4_BASETYPE x,y,z, xi,yi,zi,qj, one_over_r;
		#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
	   	||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
		V4_BASETYPE qi;
		#endif
		#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	   	||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))	
		V4_BASETYPE one_o_r_3;
		#endif
		#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	   	||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		V4_BASETYPE qj_o_r_3;
		#endif
		#if (FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD)
		V4_BASETYPE qi_o_r_3;
		#endif	
		#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   	||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))	
		V4_BASETYPE mxj, myj, mzj;
		V4_BASETYPE m_times_rj;	
		#endif
		#if (FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)
		V4_BASETYPE mxi, myi, mzi;
		V4_BASETYPE m_times_ri;	
		#endif
		#if ((FMM_KIND==FMM_DIPOLE_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		V4_BASETYPE one_o_r_5;
		#endif	
		
		#ifdef PERIODIC
		xi = access_tx(ii)[k] - dx0;
		yi = access_ty(ii)[k] - dy0;
		zi = access_tz(ii)[k] - dz0;
		#else
		xi = access_tx(ii)[k];
		yi = access_ty(ii)[k];
		zi = access_tz(ii)[k];
		#endif
		#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
		   ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
		qi = access_q(ii)[k];
		#endif
		#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
		mxi = access_mx(ii)[k];
		myi = access_my(ii)[k];
		mzi = access_mz(ii)[k];		
		#endif
		
		if (j00&3) {
		    jj0--;
		    jjj = (jj0==jj1 ? j1&3 : 4);
		    for (j=j00&3; j<jjj; j++) {
			x = xi - access_x(jj0)[j];
			y = yi - access_y(jj0)[j];
			z = zi - access_z(jj0)[j];

			#if (EVAL_DIRECT_ACCURACY==0)
			one_over_r = RECIP_SQRT0(x*x + y*y + z*z); 
			#elif (EVAL_DIRECT_ACCURACY==1)
			one_over_r = RECIP_SQRT1(x*x + y*y + z*z); 
			#elif (EVAL_DIRECT_ACCURACY==2)
			one_over_r = RECIP_SQRT2(x*x + y*y + z*z); 
			#endif
			
			qj = access_q(jj0)[j];
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
			  ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			mxj = access_mx(jj0)[j];
			myj = access_my(jj0)[j];
			mzj = access_mz(jj0)[j];
			one_o_r_3 = one_over_r*one_over_r*one_over_r;
			m_times_rj = x*mxj + y*myj + z*mzj;
			access_pot(ii)[k] += qj*one_over_r + m_times_rj*one_o_r_3;
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
			m_times_ri = x*mxi + y*myi + z*mzi;
			access_pot(jj0)[j] += qi*one_over_r - m_times_ri*one_o_r_3;
			#endif
			#else
			access_pot(ii)[k] += qj*one_over_r;
			#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD)) 
			access_pot(jj0)[j] += qi*one_over_r;
			#endif
			#endif
			
			#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_ST_GRAD))
			one_o_r_3 = one_over_r*one_over_r*one_over_r;
			qj_o_r_3 = qj*one_o_r_3;
			access_gradx(ii)[k] -= x*qj_o_r_3;
			access_grady(ii)[k] -= y*qj_o_r_3;
			access_gradz(ii)[k] -= z*qj_o_r_3;
			#endif
			#if (FMM_KIND==FMM_GRAD)
			qi_o_r_3 = qi*one_o_r_3;
			access_gradx(jj0)[j] += x*qi_o_r_3;
			access_grady(jj0)[j] += y*qi_o_r_3;
			access_gradz(jj0)[j] += z*qi_o_r_3;
			#endif
			#if ((FMM_KIND==FMM_DIPOLE_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			one_o_r_5 = one_o_r_3*one_over_r*one_over_r;
			qj_o_r_3 = qj*one_o_r_3 + 3.0*m_times_rj*one_o_r_5;
			access_gradx(ii)[k] -= x*qj_o_r_3 - mxj*one_o_r_3;
			access_grady(ii)[k] -= y*qj_o_r_3 - myj*one_o_r_3;
			access_gradz(ii)[k] -= z*qj_o_r_3 - mzj*one_o_r_3;
			#endif
			#if (FMM_KIND==FMM_DIPOLE_GRAD)
			qi_o_r_3 = qi*one_o_r_3 - 3.0*m_times_ri*one_o_r_5;	
			access_gradx(jj0)[j] += x*qi_o_r_3 + mxi*one_o_r_3;
			access_grady(jj0)[j] += y*qi_o_r_3 + myi*one_o_r_3;
			access_gradz(jj0)[j] += z*qi_o_r_3 + mzi*one_o_r_3;	
			#endif
		    }
		    if (jj0==jj1) continue;

	    	}

    		for (j=0; j<(j1&3); j++) {
			x = xi - access_x(jj1)[j];
			y = yi - access_y(jj1)[j];
			z = zi - access_z(jj1)[j];

			#if (EVAL_DIRECT_ACCURACY==0)
			one_over_r = RECIP_SQRT0(x*x + y*y + z*z); 
			#elif (EVAL_DIRECT_ACCURACY==1)
			one_over_r = RECIP_SQRT1(x*x + y*y + z*z); 
			#elif (EVAL_DIRECT_ACCURACY==2)
			one_over_r = RECIP_SQRT2(x*x + y*y + z*z); 
			#endif
			
			qj = access_q(jj1)[j];
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
			  ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			mxj = access_mx(jj1)[j];
			myj = access_my(jj1)[j];
			mzj = access_mz(jj1)[j];
			one_o_r_3 = one_over_r*one_over_r*one_over_r;
			m_times_rj = x*mxj + y*myj + z*mzj;
			access_pot(ii)[k] += qj*one_over_r + m_times_rj*one_o_r_3;
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
			m_times_ri = x*mxi + y*myi + z*mzi;
			access_pot(jj1)[j] += qi*one_over_r - m_times_ri*one_o_r_3;
			#endif
			#else
			access_pot(ii)[k] += qj*one_over_r;
			#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD)) 
			access_pot(jj1)[j] += qi*one_over_r;
			#endif
			#endif
			
			#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_ST_GRAD))
			one_o_r_3 = one_over_r*one_over_r*one_over_r;
			qj_o_r_3 = qj*one_o_r_3;
			access_gradx(ii)[k] -= x*qj_o_r_3;
			access_grady(ii)[k] -= y*qj_o_r_3;
			access_gradz(ii)[k] -= z*qj_o_r_3;
			#endif
			#if (FMM_KIND==FMM_GRAD)
			qi_o_r_3 = qi*one_o_r_3;
			access_gradx(jj1)[j] += x*qi_o_r_3;
			access_grady(jj1)[j] += y*qi_o_r_3;
			access_gradz(jj1)[j] += z*qi_o_r_3;
			#endif
			#if ((FMM_KIND==FMM_DIPOLE_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			one_o_r_5 = one_o_r_3*one_over_r*one_over_r;
			qj_o_r_3 = qj*one_o_r_3 + 3.0*m_times_rj*one_o_r_5;
			access_gradx(ii)[k] -= x*qj_o_r_3 - mxj*one_o_r_3;
			access_grady(ii)[k] -= y*qj_o_r_3 - myj*one_o_r_3;
			access_gradz(ii)[k] -= z*qj_o_r_3 - mzj*one_o_r_3;
			#endif
			#if (FMM_KIND==FMM_DIPOLE_GRAD)
			qi_o_r_3 = qi*one_o_r_3 - 3.0*m_times_ri*one_o_r_5;	
			access_gradx(jj1)[j] += x*qi_o_r_3 + mxi*one_o_r_3;
			access_grady(jj1)[j] += y*qi_o_r_3 + myi*one_o_r_3;
			access_gradz(jj1)[j] += z*qi_o_r_3 + mzi*one_o_r_3;	
			#endif
	       }
	    }
	}
    }
}


//*******************************************************************************
// TODO: better handling of simd-exponential function
#if (FMM_PRECISION==0)
__m128 exp_ps(__m128 x); 
__m128 exp0_ps(__m128 x); 
#define V4_EXP2(x) exp_ps(x)

#define V4_EXP1(x) V4_EXP2(x)
#define V4_EXP0(x) exp0_ps(x)

#else

#define V4_EXP2(x) V4_SET(exp(((V4_BASETYPE*) &(x))[3]), \
                          exp(((V4_BASETYPE*) &(x))[2]), \
                          exp(((V4_BASETYPE*) &(x))[1]), \
                          exp(((V4_BASETYPE*) &(x))[0]))

#define V4_EXP1(x) V4_EXP2(x)
#define V4_EXP0(x) V4_EXP2(x)

#endif

#ifdef PERIODIC
void eval_direct_yukawa_periodic(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx0, _FLOAT_ dy0, _FLOAT_ dz0)
#else
void eval_direct_yukawa(FmmvHandle *FMMV, Box *target, Box *source)
#endif
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)
        
        V4_TYPE beta = V4_SET1(FMMV->beta);

        V4_TYPE beta_r, minus_beta_r, exp_minus_beta_r_over_r;

	V4_TYPE x, y, z, one_over_r, r2, poti;
	V4_TYPE xi,yi,zi, qj;
	#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
	   ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
	V4_TYPE qi;
	#endif
	int i,j,ni,nj,i0,j0,j00,j1,k;
	int ii, jj, jj0, jj1, jjj;
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	V4_TYPE hh0;
	#endif
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	V4_TYPE gradxi, gradyi, gradzi;
	V4_TYPE hh1;
	#endif
	#if ((FMM_KIND==FMM_DIPOLE_GRAD) \
	   ||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	V4_TYPE  hh2;
	V4_TYPE beta2 = V4_MUL(beta, beta);	
	#endif
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	V4_TYPE mxj, myj, mzj;
	V4_TYPE m_times_rj;	
	#endif
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
	V4_TYPE mxi, myi, mzi;
	V4_TYPE m_times_ri;	
	#endif

#ifdef PERIODIC
	V4_TYPE dx, dy, dz;

	dx = V4_SET1(dx0);
	dy = V4_SET1(dy0);
	dz = V4_SET1(dz0);
#endif	

	i0 = target->firstTarget;
	ni = target->noOfTargets;
	j0 = source->firstParticle;
	nj = source->noOfParticles;

	#if ((FMM_KIND==FMM_ST_STANDARD)||(FMM_KIND==FMM_ST_GRAD) \
   	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	j00 =j0;
	FMMV->noOfDirectInteractions += ni*nj;
	#else
	if (i0==j0) {
		FMMV->noOfDirectInteractions += (ni*(ni-1))/2;
	}
	else { 
		FMMV->noOfDirectInteractions += (ni*(ni-1))/2;
	}
	#endif


	j1 = j0+nj;
	jj1 = j1>>2;

	for (i=i0; i<i0+ni; i++) {
		ii = i>>2;
		k = i&3;
		#ifdef PERIODIC
		xi = V4_SUB(V4_LOAD1(access_tx(ii)+k), dx);
		yi = V4_SUB(V4_LOAD1(access_ty(ii)+k), dy);
		zi = V4_SUB(V4_LOAD1(access_tz(ii)+k), dz);
		#else
		xi = V4_LOAD1(access_tx(ii)+k);
		yi = V4_LOAD1(access_ty(ii)+k);
		zi = V4_LOAD1(access_tz(ii)+k);
		#endif
		#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
		   ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
		qi = V4_LOAD1(access_q(ii)+k);
		#endif
		#if (FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)
		mxi = V4_LOAD1(access_mx(ii)+k);
		myi = V4_LOAD1(access_my(ii)+k);
		mzi = V4_LOAD1(access_mz(ii)+k);		
		#endif			
		poti = V4_LOAD0(access_pot(ii)+k);
		#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
		   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		gradxi = V4_LOAD0(access_gradx(ii)+k);
		gradyi = V4_LOAD0(access_grady(ii)+k);
		gradzi = V4_LOAD0(access_gradz(ii)+k);
		#endif

		#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
		   ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
		j00 = (i<j0 ? j0 : i+1);
		#endif
		jj0 = j00>>2;
		if (j00&3) jj0++;
		
		for (jj=jj0; jj<jj1; jj++) {
			x = V4_SUB(xi, V4_LOAD(access_x(jj)));
			y = V4_SUB(yi, V4_LOAD(access_y(jj)));
			z = V4_SUB(zi, V4_LOAD(access_z(jj)));

                        r2 = V4_FMA(x, x, V4_FMA(y, y, V4_MUL(z, z)));
			#if (EVAL_DIRECT_ACCURACY==0)
			one_over_r = V4_RECIP_SQRT0(r2); 
			#elif (EVAL_DIRECT_ACCURACY==1)
			one_over_r = V4_RECIP_SQRT1(r2); 
			#elif (EVAL_DIRECT_ACCURACY==2)
			one_over_r = V4_RECIP_SQRT2(r2); 
			#endif
                        beta_r = V4_MUL(V4_MUL(beta, r2), one_over_r);
                        minus_beta_r = V4_NEG(beta_r);
			#if (EVAL_DIRECT_ACCURACY==0)
                        exp_minus_beta_r_over_r = V4_MUL(V4_EXP0(minus_beta_r), one_over_r);
			#elif (EVAL_DIRECT_ACCURACY==1)
                        exp_minus_beta_r_over_r = V4_MUL(V4_EXP1(minus_beta_r), one_over_r);
			#elif (EVAL_DIRECT_ACCURACY==2)
                        exp_minus_beta_r_over_r = V4_MUL(V4_EXP2(minus_beta_r), one_over_r);
			#endif

			qj = V4_LOAD(access_q(jj));
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
			  ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			mxj = V4_LOAD(access_mx(jj));			
			myj = V4_LOAD(access_my(jj));			
			mzj = V4_LOAD(access_mz(jj));			
			m_times_rj = V4_FMA(x, mxj, V4_FMA(y, myj, V4_MUL(z, mzj)));

                        hh0 = V4_MUL(V4_MUL(V4_MUL(V4_ADD(V4_SET1(1.0), beta_r), one_over_r), one_over_r), exp_minus_beta_r_over_r);
			poti = V4_ADD(poti, V4_ADD(V4_MUL(qj,exp_minus_beta_r_over_r), V4_MUL(m_times_rj,hh0)));
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
			m_times_ri = V4_FMA(x, mxi, V4_FMA(y, myi, V4_MUL(z, mzi)));
                        V4_STORE(access_pot(jj), V4_ADD(V4_LOAD(access_pot(jj)), V4_SUB(V4_MUL(qi,exp_minus_beta_r_over_r), V4_MUL(m_times_ri,hh0))));
			#endif
			#else
			poti = V4_FMA(qj, exp_minus_beta_r_over_r, poti);
			#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD)) 
			V4_STORE(access_pot(jj), V4_FMA(qi, exp_minus_beta_r_over_r, V4_LOAD(access_pot(jj))));
			#endif
			#endif
			
			#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_ST_GRAD))
                        hh0 = V4_MUL(V4_MUL(V4_MUL(V4_ADD(V4_SET1(1.0), beta_r), one_over_r), one_over_r), exp_minus_beta_r_over_r);
                        hh1 = V4_MUL(qj, hh0);
			gradxi = V4_FNMS(x, hh1, gradxi);
			gradyi = V4_FNMS(y, hh1, gradyi);
			gradzi = V4_FNMS(z, hh1, gradzi);
			#endif
			#if (FMM_KIND==FMM_GRAD)
                        hh1 = V4_MUL(qi, hh0);
			V4_STORE(access_gradx(jj), V4_FMA(x, hh1, V4_LOAD(access_gradx(jj))));
			V4_STORE(access_grady(jj), V4_FMA(y, hh1, V4_LOAD(access_grady(jj))));
			V4_STORE(access_gradz(jj), V4_FMA(z, hh1, V4_LOAD(access_gradz(jj))));
			#endif
			#if ((FMM_KIND==FMM_DIPOLE_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
                        hh1 = V4_MUL(V4_MUL(one_over_r, one_over_r), m_times_rj);
                        hh2 = V4_ADD(V4_MUL(V4_ADD(qj, V4_MUL(V4_SET1(3.0),hh1)),hh0), V4_MUL(V4_MUL(beta2, hh1), exp_minus_beta_r_over_r));
			gradxi = V4_FNMS(x, hh2, V4_FMA(mxj, hh0, gradxi));
			gradyi = V4_FNMS(y, hh2, V4_FMA(myj, hh0, gradyi));
			gradzi = V4_FNMS(z, hh2, V4_FMA(mzj, hh0, gradzi));
			#endif
			#if (FMM_KIND==FMM_DIPOLE_GRAD)
                        hh1 = V4_MUL(V4_MUL(one_over_r, one_over_r),m_times_ri);
                        hh2 = V4_SUB(V4_MUL(V4_SUB(qi, V4_MUL(V4_SET1(3.0),hh1)),hh0), V4_MUL(V4_MUL(beta2, hh1), exp_minus_beta_r_over_r));
			V4_STORE(access_gradx(jj), V4_FMA(x, hh2, V4_FMA(mxi, hh0, V4_LOAD(access_gradx(jj)))));
			V4_STORE(access_grady(jj), V4_FMA(y, hh2, V4_FMA(myi, hh0, V4_LOAD(access_grady(jj)))));
			V4_STORE(access_gradz(jj), V4_FMA(z, hh2, V4_FMA(mzi, hh0, V4_LOAD(access_gradz(jj)))));
			#endif
                }

		V4_STORE0(access_pot(ii)+k, V4_HORIZADD1(poti));
		#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
		   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		V4_STORE0(access_gradx(ii)+k, V4_HORIZADD1(gradxi));
		V4_STORE0(access_grady(ii)+k, V4_HORIZADD1(gradyi));
		V4_STORE0(access_gradz(ii)+k, V4_HORIZADD1(gradzi));
		#endif

	    if ((j1>j00)&&((j00&3)||(j1&3))) {
                V4_BASETYPE beta = FMMV->beta;
                V4_BASETYPE beta_r, exp_minus_beta_r_over_r;
	        V4_BASETYPE x, y, z, one_over_r, r2;
	        V4_BASETYPE xi,yi,zi, qj;
	        #if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
	           ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
	        V4_BASETYPE qi;
	        #endif
	        #if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	           ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	        V4_BASETYPE hh0;
	        #endif
	        #if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	           ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	        V4_BASETYPE hh1;
	        #endif
	        #if ((FMM_KIND==FMM_DIPOLE_GRAD) \
	           ||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	        V4_BASETYPE  hh2;
	        V4_BASETYPE beta2 = beta*beta;	
	        #endif
	        #if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	           ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	        V4_BASETYPE mxj, myj, mzj;
	        V4_BASETYPE m_times_rj;	
	        #endif
	        #if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
	        V4_BASETYPE mxi, myi, mzi;
	        V4_BASETYPE m_times_ri;	
	        #endif

		#ifdef PERIODIC
		xi = access_tx(ii)[k] - dx0;
		yi = access_ty(ii)[k] - dy0;
		zi = access_tz(ii)[k] - dz0;
		#else
		xi = access_tx(ii)[k];
		yi = access_ty(ii)[k];
		zi = access_tz(ii)[k];
		#endif
		#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
		   ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
		qi = access_q(ii)[k];
		#endif
		#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
		mxi = access_mx(ii)[k];
		myi = access_my(ii)[k];
		mzi = access_mz(ii)[k];		
		#endif
		
		if (j00&3) {
		    jj0--;
		    jjj = (jj0==jj1 ? j1&3 : 4);
		    for (j=j00&3; j<jjj; j++) {
			x = xi - access_x(jj0)[j];
			y = yi - access_y(jj0)[j];
			z = zi - access_z(jj0)[j];

                        r2 = x*x + y*y + z*z;
			#if (EVAL_DIRECT_ACCURACY==0)
			one_over_r = RECIP_SQRT0(r2); 
			#elif (EVAL_DIRECT_ACCURACY==1)
			one_over_r = RECIP_SQRT1(r2); 
			#elif (EVAL_DIRECT_ACCURACY==2)
			one_over_r = RECIP_SQRT2(r2); 
			#endif
                        beta_r = beta*r2*one_over_r;
			#if (EVAL_DIRECT_ACCURACY==0)
                        exp_minus_beta_r_over_r = EXP0(-beta_r)*one_over_r;
			#elif (EVAL_DIRECT_ACCURACY==1)
                        exp_minus_beta_r_over_r = EXP1(-beta_r)*one_over_r;
			#elif (EVAL_DIRECT_ACCURACY==2)
                        exp_minus_beta_r_over_r = EXP2(-beta_r)*one_over_r;
			#endif


			qj = access_q(jj0)[j];
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
			  ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			mxj = access_mx(jj0)[j];
			myj = access_my(jj0)[j];
			mzj = access_mz(jj0)[j];
			m_times_rj = x*mxj + y*myj + z*mzj;
                        hh0 = (1.0+beta_r)*one_over_r*one_over_r*exp_minus_beta_r_over_r;
			access_pot(ii)[k] += qj*exp_minus_beta_r_over_r + m_times_rj*hh0;
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
			m_times_ri = x*mxi + y*myi + z*mzi;
			access_pot(jj0)[j] += qi*exp_minus_beta_r_over_r - m_times_ri*hh0;
			#endif
			#else
			access_pot(ii)[k] += qj*exp_minus_beta_r_over_r;
			#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD)) 
			access_pot(jj0)[j] += qi*exp_minus_beta_r_over_r;
			#endif
			#endif
			
			#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_ST_GRAD))
                        hh0 = (1.0+beta_r)*one_over_r*one_over_r*exp_minus_beta_r_over_r;
                        hh1 = qj*hh0;
			access_gradx(ii)[k] -= x*hh1;
			access_grady(ii)[k] -= y*hh1;
			access_gradz(ii)[k] -= z*hh1;
			#endif
			#if (FMM_KIND==FMM_GRAD)
                        hh1 = qi*hh0;
			access_gradx(jj0)[j] += x*hh1;
			access_grady(jj0)[j]+= y*hh1;
			access_gradz(jj0)[j]+= z*hh1;
			#endif
			#if ((FMM_KIND==FMM_DIPOLE_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
                        hh1 = one_over_r*one_over_r*m_times_rj;
                        hh2 = (qj + 3.0*hh1)*hh0 + beta2*hh1*exp_minus_beta_r_over_r;
			access_gradx(ii)[k] -= x*hh2 - mxj*hh0;
			access_grady(ii)[k]-= y*hh2 - myj*hh0;
			access_gradz(ii)[k] -= z*hh2 - mzj*hh0;
			#endif
			#if (FMM_KIND==FMM_DIPOLE_GRAD)
                        hh1 = one_over_r*one_over_r*m_times_ri;
                        hh2 = (qi - 3.0*hh1)*hh0 - beta2*hh1*exp_minus_beta_r_over_r;
			access_gradx(jj0)[j] += x*hh2 + mxi*hh0;
			access_grady(jj0)[j] += y*hh2 + myi*hh0;
			access_gradz(jj0)[j] += z*hh2 + mzi*hh0;	
			#endif
		    }
		    if (jj0==jj1) continue;

	    	}

    		for (j=0; j<(j1&3); j++) {
			x = xi - access_x(jj1)[j];
			y = yi - access_y(jj1)[j];
			z = zi - access_z(jj1)[j];

                        r2 = x*x + y*y + z*z;
			#if (EVAL_DIRECT_ACCURACY==0)
			one_over_r = RECIP_SQRT0(r2); 
			#elif (EVAL_DIRECT_ACCURACY==1)
			one_over_r = RECIP_SQRT1(r2); 
			#elif (EVAL_DIRECT_ACCURACY==2)
			one_over_r = RECIP_SQRT2(r2); 
			#endif
                        beta_r = beta*r2*one_over_r;
			#if (EVAL_DIRECT_ACCURACY==0)
                        exp_minus_beta_r_over_r = EXP0(-beta_r)*one_over_r;
			#elif (EVAL_DIRECT_ACCURACY==1)
                        exp_minus_beta_r_over_r = EXP1(-beta_r)*one_over_r;
			#elif (EVAL_DIRECT_ACCURACY==2)
                        exp_minus_beta_r_over_r = EXP2(-beta_r)*one_over_r;
			#endif


			qj = access_q(jj1)[j];
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
			  ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			mxj = access_mx(jj1)[j];
			myj = access_my(jj1)[j];
			mzj = access_mz(jj1)[j];
			m_times_rj = x*mxj + y*myj + z*mzj;
                        hh0 = (1.0+beta_r)*one_over_r*one_over_r*exp_minus_beta_r_over_r;
			access_pot(ii)[k] += qj*exp_minus_beta_r_over_r + m_times_rj*hh0;
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
			m_times_ri = x*mxi + y*myi + z*mzi;
			access_pot(jj1)[j] += qi*exp_minus_beta_r_over_r - m_times_ri*hh0;
			#endif
			#else
			access_pot(ii)[k] += qj*exp_minus_beta_r_over_r;
			#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD)) 
			access_pot(jj1)[j] += qi*exp_minus_beta_r_over_r;
			#endif
			#endif
			
			#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_ST_GRAD))
                        hh0 = (1.0+beta_r)*one_over_r*one_over_r*exp_minus_beta_r_over_r;
                        hh1 = qj*hh0;
			access_gradx(ii)[k] -= x*hh1;
			access_grady(ii)[k] -= y*hh1;
			access_gradz(ii)[k] -= z*hh1;
			#endif
			#if (FMM_KIND==FMM_GRAD)
                        hh1 = qi*hh0;
			access_gradx(jj1)[j] += x*hh1;
			access_grady(jj1)[j] += y*hh1;
			access_gradz(jj1)[j] += z*hh1;
			#endif
			#if ((FMM_KIND==FMM_DIPOLE_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
                        hh1 = one_over_r*one_over_r*m_times_rj;
                        hh2 = (qj + 3.0*hh1)*hh0 + beta2*hh1*exp_minus_beta_r_over_r;
			access_gradx(ii)[k] -= x*hh2 - mxj*hh0;
			access_grady(ii)[k] -= y*hh2 - myj*hh0;
			access_gradz(ii)[k] -= z*hh2 - mzj*hh0;
			#endif
			#if (FMM_KIND==FMM_DIPOLE_GRAD)
                        hh1 = one_over_r*one_over_r*m_times_ri;
                        hh2 = (qi - 3.0*hh1)*hh0 - beta2*hh1*exp_minus_beta_r_over_r;
			access_gradx(jj1)[j] += x*hh2 + mxi*hh0;
			access_grady(jj1)[j] += y*hh2 + myi*hh0;
			access_gradz(jj1)[j] += z*hh2 + mzi*hh0;	
			#endif
	    	}

            }
	}   
}














































