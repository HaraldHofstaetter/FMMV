/*
 * FMMV - the Fastest Multipole Method of Vienna
 * Copyright (c) 2006-2010 Harald Hofstaetter
 * Institute of Mathematics
 * University of Vienna
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

#if (FMM_PRECISION==0)
   #include"simd2s.h"
#else
   #include"simd2d.h"
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

	V2_TYPE x,y,z, xi,yi,zi,qj, one_over_r, poti;
	#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
	   ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
	V2_TYPE qi, potj;
	#endif
	int i,j,k,ni,nj,i0,j0, j00,j1;
	int ii, jj, jj0, jj1, left, right;
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))	
	V2_TYPE one_o_r_3;
	#endif
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	V2_TYPE qj_o_r_3;
	V2_TYPE gradxi, gradyi, gradzi;
	#endif
	#if (FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD)
	V2_TYPE qi_o_r_3;
	V2_TYPE gradxj, gradyj, gradzj;
	#endif	
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))	
	V2_TYPE mxj, myj, mzj;
	V2_TYPE m_times_rj;	
	#endif
	#if (FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)
	V2_TYPE mxi, myi, mzi;
	V2_TYPE m_times_ri;	
	#endif
	#if ((FMM_KIND==FMM_DIPOLE_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	V2_TYPE three = V2_SET1(3.0);
	V2_TYPE one_o_r_5;
	#endif	
	
#ifdef PERIODIC
	V2_TYPE dx, dy, dz;

	dx = V2_SET1(dx0);
	dy = V2_SET1(dy0);
	dz = V2_SET1(dz0);
#endif	
	i0 = target->firstTarget;
	ni = target->noOfTargets;
	j0 = source->firstParticle;
	nj = source->noOfParticles;

        #if ((FMM_KIND==FMM_ST_STANDARD)||(FMM_KIND==FMM_ST_GRAD) \
           ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
        FMMV->noOfDirectInteractions += ni*nj;
	j00 =j0;
	jj0 = j00>>1;
	if (j00&1) jj0++;
        #else
        if (i0==j0) {
                FMMV->noOfDirectInteractions += (ni*(ni-1))/2;
        }
        else { 
                FMMV->noOfDirectInteractions += ni*nj;
        }
        #endif

	j1 = j0+nj;
	jj1 = j1>>1;
	
	for (i=i0; i<i0+ni; i++) {
		ii = i>>1;
		k = i&1;
		#ifdef PERIODIC
		xi = V2_SUB(V2_LOAD1(access_tx(ii)+k), dx);
		yi = V2_SUB(V2_LOAD1(access_ty(ii)+k), dy);
		zi = V2_SUB(V2_LOAD1(access_tz(ii)+k), dz);
		#else
		xi = V2_LOAD1(access_tx(ii)+k);
		yi = V2_LOAD1(access_ty(ii)+k);
		zi = V2_LOAD1(access_tz(ii)+k);
		#endif
		#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
		   ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
		qi = V2_LOAD1(access_q(ii)+k);
		#endif
		#if (FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)
		mxi = V2_LOAD1(access_mx(ii)+k);
		myi = V2_LOAD1(access_my(ii)+k);
		mzi = V2_LOAD1(access_mz(ii)+k);		
		#endif			
		poti = V2_LOAD0(access_pot(ii)+k);
		#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
		   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		gradxi = V2_LOAD0(access_gradx(ii)+k);
		gradyi = V2_LOAD0(access_grady(ii)+k);
		gradzi = V2_LOAD0(access_gradz(ii)+k);
		#endif
		
		#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
		   ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
		j00 = (i<j0 ? j0 : i+1);
		jj0 = j00>>1;
		if (j00&1) jj0++;
		#endif

		for (jj=jj0; jj<jj1; jj++) {
			x = V2_SUB(xi, V2_LOAD(access_x(jj)));
			y = V2_SUB(yi, V2_LOAD(access_y(jj)));
			z = V2_SUB(zi, V2_LOAD(access_z(jj)));
			
			one_over_r = 
			#if (EVAL_DIRECT_ACCURACY==0)
			V2_RECIP_SQRT0(V2_FMA(x, x, V2_FMA(y, y, V2_MUL(z, z))));
			#elif (EVAL_DIRECT_ACCURACY==1)
			V2_RECIP_SQRT1(V2_FMA(x, x, V2_FMA(y, y, V2_MUL(z, z))));
			#elif (EVAL_DIRECT_ACCURACY==2)
			V2_RECIP_SQRT2(V2_FMA(x, x, V2_FMA(y, y, V2_MUL(z, z))));
			#endif
			
			qj = V2_LOAD(access_q(jj));
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
			  ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			mxj = V2_LOAD(access_mx(jj));			
			myj = V2_LOAD(access_my(jj));			
			mzj = V2_LOAD(access_mz(jj));			
			one_o_r_3 = V2_MUL(V2_MUL(one_over_r, one_over_r), one_over_r);
			m_times_rj = V2_FMA(x, mxj, V2_FMA(y, myj, V2_MUL(z, mzj)));
			poti = V2_ADD(poti, V2_FMA(qj, one_over_r, V2_MUL(m_times_rj,one_o_r_3)));
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
			m_times_ri = V2_FMA(x, mxi, V2_FMA(y, myi, V2_MUL(z, mzi)));
			V2_STORE(access_pot(jj), V2_ADD(V2_LOAD(access_pot(jj)), V2_FNMS(m_times_ri,one_o_r_3, V2_MUL(qi, one_over_r))));
			#endif
			#else
			poti = V2_FMA(qj, one_over_r, poti);
			#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD)) 
			V2_STORE(access_pot(jj), V2_FMA(qi, one_over_r, V2_LOAD(access_pot(jj))));
			#endif
			#endif
			
			#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_ST_GRAD))
			one_o_r_3 = V2_MUL(V2_MUL(one_over_r, one_over_r), one_over_r);
			qj_o_r_3 = V2_MUL(qj, one_o_r_3);
			gradxi = V2_FNMS(x, qj_o_r_3, gradxi);
			gradyi = V2_FNMS(y, qj_o_r_3, gradyi);
			gradzi = V2_FNMS(z, qj_o_r_3, gradzi);
			#endif
			#if (FMM_KIND==FMM_GRAD)
			qi_o_r_3 = V2_MUL(qi, one_o_r_3);
			V2_STORE(access_gradx(jj), V2_FMA(x, qi_o_r_3, V2_LOAD(access_gradx(jj))));
			V2_STORE(access_grady(jj), V2_FMA(y, qi_o_r_3, V2_LOAD(access_grady(jj))));
			V2_STORE(access_gradz(jj), V2_FMA(z, qi_o_r_3, V2_LOAD(access_gradz(jj))));		
			#endif
			#if ((FMM_KIND==FMM_DIPOLE_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			one_o_r_5 = V2_MUL(one_o_r_3, V2_MUL(one_over_r, one_over_r));
			qj_o_r_3 = V2_FMA(qj, one_o_r_3, V2_MUL(three, V2_MUL(m_times_rj, one_o_r_5)));
			gradxi = V2_FNMS(x, qj_o_r_3, V2_FMA(mxj, one_o_r_3, gradxi));
			gradyi = V2_FNMS(y, qj_o_r_3, V2_FMA(myj, one_o_r_3, gradyi));
			gradzi = V2_FNMS(z, qj_o_r_3, V2_FMA(mzj, one_o_r_3, gradzi));
			#endif
			#if (FMM_KIND==FMM_DIPOLE_GRAD)
			qi_o_r_3 = V2_FNMS(three, V2_MUL(m_times_ri, one_o_r_5), V2_MUL(qi, one_o_r_3));
			V2_STORE(access_gradx(jj), V2_FMA(x, qi_o_r_3, V2_FMA(mxi, one_o_r_3, V2_LOAD(access_gradx(jj)))));
			V2_STORE(access_grady(jj), V2_FMA(y, qi_o_r_3, V2_FMA(myi, one_o_r_3, V2_LOAD(access_grady(jj)))));
			V2_STORE(access_gradz(jj), V2_FMA(z, qi_o_r_3, V2_FMA(mzi, one_o_r_3, V2_LOAD(access_gradz(jj)))));
			#endif
		}
		
		#if ((FMM_KIND==FMM_ST_STANDARD)||(FMM_KIND==FMM_ST_GRAD) \
		  ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		right = (j00&1)&&1;  /* &&1 to get standard boolean value */
		left = (j1&1)&&(j1>j00); 
	        #else	
		right = (j00&1)&&(j1>j00)&&(j00>i);  
		left = (j1&1)&&(j1>j00)&&(j1-1>i); 
		#endif

		if (right&&left) {
			x = V2_SUB(xi, V2_LOAD2(access_x(jj0-1)+1, access_x(jj1)));
			y = V2_SUB(yi, V2_LOAD2(access_y(jj0-1)+1, access_y(jj1)));
			z = V2_SUB(zi, V2_LOAD2(access_z(jj0-1)+1, access_z(jj1)));
			
			one_over_r = 
			#if (EVAL_DIRECT_ACCURACY==0)
			V2_RECIP_SQRT0(V2_FMA(x, x, V2_FMA(y, y, V2_MUL(z, z))));
			#elif (EVAL_DIRECT_ACCURACY==1)
			V2_RECIP_SQRT1(V2_FMA(x, x, V2_FMA(y, y, V2_MUL(z, z))));
			#elif (EVAL_DIRECT_ACCURACY==2)
			V2_RECIP_SQRT2(V2_FMA(x, x, V2_FMA(y, y, V2_MUL(z, z))));
			#endif

			qj = V2_LOAD2(access_q(jj0-1)+1, access_q(jj1));
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
			  ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			mxj = V2_LOAD2(access_mx(jj0-1)+1, access_mx(jj1));
			myj = V2_LOAD2(access_my(jj0-1)+1, access_my(jj1));
			mzj = V2_LOAD2(access_mz(jj0-1)+1, access_mz(jj1));
			one_o_r_3 = V2_MUL(V2_MUL(one_over_r, one_over_r), one_over_r);
			m_times_rj = V2_FMA(x, mxj, V2_FMA(y, myj, V2_MUL(z, mzj)));
			poti = V2_ADD(poti, V2_FMA(qj, one_over_r, V2_MUL(m_times_rj,one_o_r_3)));
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
			m_times_ri = V2_FMA(x, mxi, V2_FMA(y, myi, V2_MUL(z, mzi)));
			potj = V2_ADD(V2_LOAD2(access_pot(jj0-1)+1, access_pot(jj1)), V2_FNMS(m_times_ri,one_o_r_3, V2_MUL(qi, one_over_r)));
			V2_STORE2(access_pot(jj0-1)+1, access_pot(jj1), potj);
			#endif
			#else
			poti = V2_FMA(qj, one_over_r, poti);
			#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD)) 
			potj= V2_FMA(qi, one_over_r, V2_LOAD2(access_pot(jj0-1)+1, access_pot(jj1)));
			V2_STORE2(access_pot(jj0-1)+1, access_pot(jj1), potj);
			#endif
			#endif
			
			#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_ST_GRAD))
			one_o_r_3 = V2_MUL(V2_MUL(one_over_r, one_over_r), one_over_r);
			qj_o_r_3 = V2_MUL(qj, one_o_r_3);
			gradxi = V2_FNMS(x, qj_o_r_3, gradxi);
			gradyi = V2_FNMS(y, qj_o_r_3, gradyi);
			gradzi = V2_FNMS(z, qj_o_r_3, gradzi);
			#endif
			#if (FMM_KIND==FMM_GRAD)
			qi_o_r_3 = V2_MUL(qi, one_o_r_3);	
			gradxj = V2_FMA(x, qi_o_r_3, V2_LOAD2(access_gradx(jj0-1)+1, access_gradx(jj1)));
			V2_STORE2(access_gradx(jj0-1)+1, access_gradx(jj1), gradxj);
			gradyj = V2_FMA(y, qi_o_r_3, V2_LOAD2(access_grady(jj0-1)+1, access_grady(jj1)));
			V2_STORE2(access_grady(jj0-1)+1, access_grady(jj1), gradyj);
			gradzj = V2_FMA(z, qi_o_r_3, V2_LOAD2(access_gradz(jj0-1)+1, access_gradz(jj1)));
			V2_STORE2(access_gradz(jj0-1)+1, access_gradz(jj1), gradzj);
			#endif
			#if ((FMM_KIND==FMM_DIPOLE_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			one_o_r_5 = V2_MUL(one_o_r_3, V2_MUL(one_over_r, one_over_r));
			qj_o_r_3 = V2_FMA(qj, one_o_r_3, V2_MUL(three, V2_MUL(m_times_rj, one_o_r_5)));
			gradxi = V2_FNMS(x, qj_o_r_3, V2_FMA(mxj, one_o_r_3, gradxi));
			gradyi = V2_FNMS(y, qj_o_r_3, V2_FMA(myj, one_o_r_3, gradyi));
			gradzi = V2_FNMS(z, qj_o_r_3, V2_FMA(mzj, one_o_r_3, gradzi));
			#endif
			#if (FMM_KIND==FMM_DIPOLE_GRAD)
			qi_o_r_3 = V2_FNMS(three, V2_MUL(m_times_ri, one_o_r_5), V2_MUL(qi, one_o_r_3));
			gradxj = V2_FMA(x, qi_o_r_3, V2_FMA(mxi, one_o_r_3, V2_LOAD2(access_gradx(jj0-1)+1, access_gradx(jj1))));
			V2_STORE2(access_gradx(jj0-1)+1, access_gradx(jj1), gradxj);
			gradyj = V2_FMA(y, qi_o_r_3, V2_FMA(myi, one_o_r_3, V2_LOAD2(access_grady(jj0-1)+1, access_grady(jj1))));
			V2_STORE2(access_grady(jj0-1)+1, access_grady(jj1), gradyj);
			gradzj = V2_FMA(z, qi_o_r_3, V2_FMA(mzi, one_o_r_3, V2_LOAD2(access_gradz(jj0-1)+1, access_gradz(jj1))));
			V2_STORE2(access_gradz(jj0-1)+1, access_gradz(jj1), gradzj);
			#endif
		}

		V2_STORE0(access_pot(ii)+k, V2_HORIZADD1(poti));
		#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
		   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		V2_STORE0(access_gradx(ii)+k, V2_HORIZADD1(gradxi));
		V2_STORE0(access_grady(ii)+k, V2_HORIZADD1(gradyi));
		V2_STORE0(access_gradz(ii)+k, V2_HORIZADD1(gradzi));
		#endif
		
		if (left!=right) { /* != ... xor */
		        V2_BASETYPE x,y,z, xi,yi,zi,qj, one_over_r;
		        #if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
	   	        ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
		        V2_BASETYPE qi;
		        #endif
		        #if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	   	        ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))	
		        V2_BASETYPE one_o_r_3;
		        #endif
		        #if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	   	        ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		        V2_BASETYPE qj_o_r_3;
		        #endif
		        #if (FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD)
		        V2_BASETYPE qi_o_r_3;
		        #endif	
		        #if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   	        ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))	
		        V2_BASETYPE mxj, myj, mzj;
		        V2_BASETYPE m_times_rj;	
		        #endif
		        #if (FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)
		        V2_BASETYPE mxi, myi, mzi;
		        V2_BASETYPE m_times_ri;	
		        #endif
		        #if ((FMM_KIND==FMM_DIPOLE_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		        V2_BASETYPE one_o_r_5;
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
    
		        if (!left) { 	
			        jj = jj0-1;
			        j = 1;
		        }
		        else {
			        jj = jj1;
			        j = 0;
		        }	    

			x = xi - access_x(jj)[j];
			y = yi - access_y(jj)[j];
			z = zi - access_z(jj)[j];

			#if (EVAL_DIRECT_ACCURACY==0)
			one_over_r = RECIP_SQRT0(x*x + y*y + z*z); 
			#elif (EVAL_DIRECT_ACCURACY==1)
			one_over_r = RECIP_SQRT1(x*x + y*y + z*z); 
			#elif (EVAL_DIRECT_ACCURACY==2)
			one_over_r = RECIP_SQRT2(x*x + y*y + z*z); 
			#endif
			
			qj = access_q(jj)[j];
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
			  ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			mxj = access_mx(jj)[j];
			myj = access_my(jj)[j];
			mzj = access_mz(jj)[j];
			one_o_r_3 = one_over_r*one_over_r*one_over_r;
			m_times_rj = x*mxj + y*myj + z*mzj;
			access_pot(ii)[k] += qj*one_over_r + m_times_rj*one_o_r_3;
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
			m_times_ri = x*mxi + y*myi + z*mzi;
			access_pot(jj)[j] += qi*one_over_r - m_times_ri*one_o_r_3;
			#endif
			#else
			access_pot(ii)[k] += qj*one_over_r;
			#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD)) 
			access_pot(jj)[j] += qi*one_over_r;
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
			access_gradx(jj)[j] += x*qi_o_r_3;
			access_grady(jj)[j] += y*qi_o_r_3;
			access_gradz(jj)[j] += z*qi_o_r_3;
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
			access_gradx(jj)[j] += x*qi_o_r_3 + mxi*one_o_r_3;
			access_grady(jj)[j] += y*qi_o_r_3 + myi*one_o_r_3;
			access_gradz(jj)[j] += z*qi_o_r_3 + mzi*one_o_r_3;	
			#endif
		}
	}
    }
}


// TODO: other accuracies for exp ...
// NOTE: order of vector components!
#if (FMM_PRECISION==0)
#define V2_EXP2(x) V2_SET(expf(((V2_BASETYPE*) &(x))[1]), \
                         expf(((V2_BASETYPE*) &(x))[0]))
#else
#define V2_EXP2(x) V2_SET(exp(((V2_BASETYPE*) &(x))[1]), \
                         exp(((V2_BASETYPE*) &(x))[0]))
#endif


#define V2_EXP1(x) V2_EXP2(x)
#define V2_EXP0(x) V2_EXP2(x)




#ifdef PERIODIC
void eval_direct_yukawa_periodic(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx0, _FLOAT_ dy0, _FLOAT_ dz0)
#else
void eval_direct_yukawa(FmmvHandle *FMMV, Box *target, Box *source)
#endif
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)
        
        V2_TYPE beta = V2_SET1(FMMV->beta);

        V2_TYPE beta_r, minus_beta_r, exp_minus_beta_r_over_r;

	V2_TYPE x, y, z, one_over_r, r2, poti;
	V2_TYPE xi,yi,zi, qj;
	#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
	   ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
	V2_TYPE qi;
	#endif
	int i,j,ni,nj,i0,j0,j00,j1,k;
	int ii, jj, jj0, jj1, left, right;
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	V2_TYPE hh0;
	#endif
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	V2_TYPE gradxi, gradyi, gradzi;
	V2_TYPE hh1;
	#endif
	#if ((FMM_KIND==FMM_DIPOLE_GRAD) \
	   ||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	V2_TYPE  hh2;
	V2_TYPE beta2 = V2_MUL(beta, beta);	
	#endif
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	V2_TYPE mxj, myj, mzj;
	V2_TYPE m_times_rj;	
	#endif
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
	V2_TYPE mxi, myi, mzi;
	V2_TYPE m_times_ri;	
	#endif

#ifdef PERIODIC
	V2_TYPE dx, dy, dz;

	dx = V2_SET1(dx0);
	dy = V2_SET1(dy0);
	dz = V2_SET1(dz0);
#endif	

	i0 = target->firstTarget;
	ni = target->noOfTargets;
	j0 = source->firstParticle;
	nj = source->noOfParticles;

        #if ((FMM_KIND==FMM_ST_STANDARD)||(FMM_KIND==FMM_ST_GRAD) \
           ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
        FMMV->noOfDirectInteractions += ni*nj;
	j00 =j0;
	jj0 = j00>>1;
	if (j00&1) jj0++;
        #else
        if (i0==j0) {
                FMMV->noOfDirectInteractions += (ni*(ni-1))/2;
        }
        else { 
                FMMV->noOfDirectInteractions += ni*nj;
        }
        #endif

	j1 = j0+nj;
	jj1 = j1>>1;
	
	for (i=i0; i<i0+ni; i++) {
		ii = i>>1;
		k = i&1;

		#ifdef PERIODIC
		xi = V2_SUB(V2_LOAD1(access_tx(ii)+k), dx);
		yi = V2_SUB(V2_LOAD1(access_ty(ii)+k), dy);
		zi = V2_SUB(V2_LOAD1(access_tz(ii)+k), dz);
		#else
		xi = V2_LOAD1(access_tx(ii)+k);
		yi = V2_LOAD1(access_ty(ii)+k);
		zi = V2_LOAD1(access_tz(ii)+k);
		#endif
		#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
		   ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
		qi = V2_LOAD1(access_q(ii)+k);
		#endif
		#if (FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)
		mxi = V2_LOAD1(access_mx(ii)+k);
		myi = V2_LOAD1(access_my(ii)+k);
		mzi = V2_LOAD1(access_mz(ii)+k);		
		#endif			
		poti = V2_LOAD0(access_pot(ii)+k);
		#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
		   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		gradxi = V2_LOAD0(access_gradx(ii)+k);
		gradyi = V2_LOAD0(access_grady(ii)+k);
		gradzi = V2_LOAD0(access_gradz(ii)+k);
		#endif

   	        #if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
		   ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
		j00 = (i<j0 ? j0 : i+1);
		jj0 = j00>>1;
		if (j00&1) jj0++;
		#endif
	
		for (jj=jj0; jj<jj1; jj++) {
			x = V2_SUB(xi, V2_LOAD(access_x(jj)));
			y = V2_SUB(yi, V2_LOAD(access_y(jj)));
			z = V2_SUB(zi, V2_LOAD(access_z(jj)));

                        r2 = V2_FMA(x, x, V2_FMA(y, y, V2_MUL(z, z)));
			#if (EVAL_DIRECT_ACCURACY==0)
			one_over_r = V2_RECIP_SQRT0(r2); 
			#elif (EVAL_DIRECT_ACCURACY==1)
			one_over_r = V2_RECIP_SQRT1(r2); 
			#elif (EVAL_DIRECT_ACCURACY==2)
			one_over_r = V2_RECIP_SQRT2(r2); 
			#endif
                        beta_r = V2_MUL(V2_MUL(beta, r2), one_over_r);
                        minus_beta_r = V2_NEG(beta_r);
			#if (EVAL_DIRECT_ACCURACY==0)
                        exp_minus_beta_r_over_r = V2_MUL(V2_EXP0(minus_beta_r), one_over_r);
			#elif (EVAL_DIRECT_ACCURACY==1)
                        exp_minus_beta_r_over_r = V2_MUL(V2_EXP1(minus_beta_r), one_over_r);
			#elif (EVAL_DIRECT_ACCURACY==2)
                        exp_minus_beta_r_over_r = V2_EXP2(minus_beta_r);
                        exp_minus_beta_r_over_r = V2_MUL(exp_minus_beta_r_over_r, one_over_r);
			#endif

			qj = V2_LOAD(access_q(jj));
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
			  ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			mxj = V2_LOAD(access_mx(jj));			
			myj = V2_LOAD(access_my(jj));			
			mzj = V2_LOAD(access_mz(jj));			
			m_times_rj = V2_FMA(x, mxj, V2_FMA(y, myj, V2_MUL(z, mzj)));

                        hh0 = V2_MUL(V2_MUL(V2_MUL(V2_ADD(V2_SET1(1.0), beta_r), one_over_r), one_over_r), exp_minus_beta_r_over_r);
			poti = V2_ADD(poti, V2_ADD(V2_MUL(qj,exp_minus_beta_r_over_r), V2_MUL(m_times_rj,hh0)));
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
			m_times_ri = V2_FMA(x, mxi, V2_FMA(y, myi, V2_MUL(z, mzi)));
                        V2_STORE(access_pot(jj), V2_ADD(V2_LOAD(access_pot(jj)), V2_SUB(V2_MUL(qi,exp_minus_beta_r_over_r), V2_MUL(m_times_ri,hh0))));
			#endif
			#else
			poti = V2_FMA(qj, exp_minus_beta_r_over_r, poti);
			#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD)) 
			V2_STORE(access_pot(jj), V2_FMA(qi, exp_minus_beta_r_over_r, V2_LOAD(access_pot(jj))));
			#endif
			#endif
			
			#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_ST_GRAD))
                        hh0 = V2_MUL(V2_MUL(V2_MUL(V2_ADD(V2_SET1(1.0), beta_r), one_over_r), one_over_r), exp_minus_beta_r_over_r);
                        hh1 = V2_MUL(qj, hh0);
			gradxi = V2_FNMS(x, hh1, gradxi);
			gradyi = V2_FNMS(y, hh1, gradyi);
			gradzi = V2_FNMS(z, hh1, gradzi);
			#endif
			#if (FMM_KIND==FMM_GRAD)
                        hh1 = V2_MUL(qi, hh0);
			V2_STORE(access_gradx(jj), V2_FMA(x, hh1, V2_LOAD(access_gradx(jj))));
			V2_STORE(access_grady(jj), V2_FMA(y, hh1, V2_LOAD(access_grady(jj))));
			V2_STORE(access_gradz(jj), V2_FMA(z, hh1, V2_LOAD(access_gradz(jj))));
			#endif
			#if ((FMM_KIND==FMM_DIPOLE_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
                        hh1 = V2_MUL(V2_MUL(one_over_r, one_over_r), m_times_rj);
                        hh2 = V2_ADD(V2_MUL(V2_ADD(qj, V2_MUL(V2_SET1(3.0),hh1)),hh0), V2_MUL(V2_MUL(beta2, hh1), exp_minus_beta_r_over_r));
			gradxi = V2_FNMS(x, hh2, V2_FMA(mxj, hh0, gradxi));
			gradyi = V2_FNMS(y, hh2, V2_FMA(myj, hh0, gradyi));
			gradzi = V2_FNMS(z, hh2, V2_FMA(mzj, hh0, gradzi));
			#endif
			#if (FMM_KIND==FMM_DIPOLE_GRAD)
                        hh1 = V2_MUL(V2_MUL(one_over_r, one_over_r),m_times_ri);
                        hh2 = V2_SUB(V2_MUL(V2_SUB(qi, V2_MUL(V2_SET1(3.0),hh1)),hh0), V2_MUL(V2_MUL(beta2, hh1), exp_minus_beta_r_over_r));
			V2_STORE(access_gradx(jj), V2_FMA(x, hh2, V2_FMA(mxi, hh0, V2_LOAD(access_gradx(jj)))));
			V2_STORE(access_grady(jj), V2_FMA(y, hh2, V2_FMA(myi, hh0, V2_LOAD(access_grady(jj)))));
			V2_STORE(access_gradz(jj), V2_FMA(z, hh2, V2_FMA(mzi, hh0, V2_LOAD(access_gradz(jj)))));
			#endif
                }

		#if ((FMM_KIND==FMM_ST_STANDARD)||(FMM_KIND==FMM_ST_GRAD) \
		  ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		right = (j00&1)&&1;  /* &&1 to get standard boolean value */
		left = (j1&1)&&(j1>j00); 
	        #else	
		right = (j00&1)&&(j1>j00)&&(j00>i);  
		left = (j1&1)&&(j1>j00)&&(j1-1>i); 
		#endif

		if (right&&left) {
			x = V2_SUB(xi, V2_LOAD2(access_x(jj0-1)+1, access_x(jj1)));
			y = V2_SUB(yi, V2_LOAD2(access_y(jj0-1)+1, access_y(jj1)));
			z = V2_SUB(zi, V2_LOAD2(access_z(jj0-1)+1, access_z(jj1)));

                        r2 = V2_FMA(x, x, V2_FMA(y, y, V2_MUL(z, z)));
			#if (EVAL_DIRECT_ACCURACY==0)
			one_over_r = V2_RECIP_SQRT0(r2); 
			#elif (EVAL_DIRECT_ACCURACY==1)
			one_over_r = V2_RECIP_SQRT1(r2); 
			#elif (EVAL_DIRECT_ACCURACY==2)
			one_over_r = V2_RECIP_SQRT2(r2); 
			#endif
                        beta_r = V2_MUL(V2_MUL(beta, r2), one_over_r);
                        minus_beta_r = V2_NEG(beta_r);
			#if (EVAL_DIRECT_ACCURACY==0)
                        exp_minus_beta_r_over_r = V2_MUL(V2_EXP0(minus_beta_r), one_over_r);
			#elif (EVAL_DIRECT_ACCURACY==1)
                        exp_minus_beta_r_over_r = V2_MUL(V2_EXP1(minus_beta_r), one_over_r);
			#elif (EVAL_DIRECT_ACCURACY==2)
                        exp_minus_beta_r_over_r = V2_MUL(V2_EXP2(minus_beta_r), one_over_r);
			#endif

			qj = V2_LOAD2(access_q(jj0-1)+1, access_q(jj1));
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
			  ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			mxj = V2_LOAD2(access_mx(jj0-1)+1, access_mx(jj1));
			myj = V2_LOAD2(access_my(jj0-1)+1, access_my(jj1));
			mzj = V2_LOAD2(access_mz(jj0-1)+1, access_mz(jj1));
			m_times_rj = V2_FMA(x, mxj, V2_FMA(y, myj, V2_MUL(z, mzj)));

                        hh0 = V2_MUL(V2_MUL(V2_MUL(V2_ADD(V2_SET1(1.0), beta_r), one_over_r), one_over_r), exp_minus_beta_r_over_r);
			poti = V2_ADD(poti, V2_ADD(V2_MUL(qj,exp_minus_beta_r_over_r), V2_MUL(m_times_rj,hh0)));
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
			m_times_ri = V2_FMA(x, mxi, V2_FMA(y, myi, V2_MUL(z, mzi)));
                        V2_STORE2(access_pot(jj0-1)+1, access_pot(jj1), V2_ADD(V2_LOAD2(access_pot(jj0-1)+1, access_pot(jj1)), 
                                  V2_SUB(V2_MUL(qi,exp_minus_beta_r_over_r), V2_MUL(m_times_ri,hh0))));
			#endif
			#else
			poti = V2_FMA(qj, exp_minus_beta_r_over_r, poti);
			#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD)) 
			V2_STORE2(access_pot(jj0-1)+1, access_pot(jj1), 
                                  V2_FMA(qi, exp_minus_beta_r_over_r, V2_LOAD2(access_pot(jj0-1)+1, access_pot(jj1))));
			#endif
			#endif
			
			#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_ST_GRAD))
                        hh0 = V2_MUL(V2_MUL(V2_MUL(V2_ADD(V2_SET1(1.0), beta_r), one_over_r), one_over_r), exp_minus_beta_r_over_r);
                        hh1 = V2_MUL(qj, hh0);
			gradxi = V2_FNMS(x, hh1, gradxi);
			gradyi = V2_FNMS(y, hh1, gradyi);
			gradzi = V2_FNMS(z, hh1, gradzi);
			#endif
			#if (FMM_KIND==FMM_GRAD)
                        hh1 = V2_MUL(qi, hh0);
			V2_STORE2(access_gradx(jj0-1)+1, access_gradx(jj1), V2_FMA(x, hh1, V2_LOAD2(access_gradx(jj0-1)+1, access_gradx(jj1))));
			V2_STORE2(access_grady(jj0-1)+1, access_grady(jj1), V2_FMA(y, hh1, V2_LOAD2(access_grady(jj0-1)+1, access_grady(jj1))));
			V2_STORE2(access_gradz(jj0-1)+1, access_gradz(jj1), V2_FMA(z, hh1, V2_LOAD2(access_gradz(jj0-1)+1, access_gradz(jj1))));
			#endif
			#if ((FMM_KIND==FMM_DIPOLE_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
                        hh1 = V2_MUL(V2_MUL(one_over_r, one_over_r), m_times_rj);
                        hh2 = V2_ADD(V2_MUL(V2_ADD(qj, V2_MUL(V2_SET1(3.0),hh1)),hh0), V2_MUL(V2_MUL(beta2, hh1), exp_minus_beta_r_over_r));
			gradxi = V2_FNMS(x, hh2, V2_FMA(mxj, hh0, gradxi));
			gradyi = V2_FNMS(y, hh2, V2_FMA(myj, hh0, gradyi));
			gradzi = V2_FNMS(z, hh2, V2_FMA(mzj, hh0, gradzi));
			#endif
			#if (FMM_KIND==FMM_DIPOLE_GRAD)
                        hh1 = V2_MUL(V2_MUL(one_over_r, one_over_r),m_times_ri);
                        hh2 = V2_SUB(V2_MUL(V2_SUB(qi, V2_MUL(V2_SET1(3.0),hh1)),hh0), V2_MUL(V2_MUL(beta2, hh1), exp_minus_beta_r_over_r));
			V2_STORE2(access_gradx(jj0-1)+1, access_gradx(jj1), V2_FMA(x, hh2, V2_FMA(mxi, hh0, V2_LOAD2(access_gradx(jj0-1)+1, access_gradx(jj1)))));
			V2_STORE2(access_grady(jj0-1)+1, access_grady(jj1), V2_FMA(y, hh2, V2_FMA(myi, hh0, V2_LOAD2(access_grady(jj0-1)+1, access_grady(jj1)))));
			V2_STORE2(access_gradz(jj0-1)+1, access_gradz(jj1), V2_FMA(z, hh2, V2_FMA(mzi, hh0, V2_LOAD2(access_gradz(jj0-1)+1, access_gradz(jj1)))));
			#endif
		}

		V2_STORE0(access_pot(ii)+k, V2_HORIZADD1(poti));
		#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
		   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		V2_STORE0(access_gradx(ii)+k, V2_HORIZADD1(gradxi));
		V2_STORE0(access_grady(ii)+k, V2_HORIZADD1(gradyi));
		V2_STORE0(access_gradz(ii)+k, V2_HORIZADD1(gradzi));
		#endif
		
		if (left!=right) { /* != ... xor */
                        V2_BASETYPE beta = FMMV->beta;
                        V2_BASETYPE beta_r, exp_minus_beta_r_over_r;
	                V2_BASETYPE x, y, z, one_over_r, r2;
	                V2_BASETYPE xi,yi,zi, qj;
	                #if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
	                   ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
	                V2_BASETYPE qi;
	                #endif
	                #if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	                   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	                V2_BASETYPE hh0;
	                #endif
	                #if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	                   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	                V2_BASETYPE hh1;
	                #endif
	                #if ((FMM_KIND==FMM_DIPOLE_GRAD) \
	                   ||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	                V2_BASETYPE  hh2;
	                V2_BASETYPE beta2 = beta*beta;	
	                #endif
	                #if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	                   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	                V2_BASETYPE mxj, myj, mzj;
	                V2_BASETYPE m_times_rj;	
	                #endif
	                #if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
	                V2_BASETYPE mxi, myi, mzi;
	                V2_BASETYPE m_times_ri;	
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
    
		        if (!left) { 	
			        jj = jj0-1;
			        j = 1;
		        }
		        else {
			        jj = jj1;
			        j = 0;
		        }	    

			x = xi - access_x(jj)[j];
			y = yi - access_y(jj)[j];
			z = zi - access_z(jj)[j];

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


			qj = access_q(jj)[j];
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
			  ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			mxj = access_mx(jj)[j];
			myj = access_my(jj)[j];
			mzj = access_mz(jj)[j];
			m_times_rj = x*mxj + y*myj + z*mzj;
                        hh0 = (1.0+beta_r)*one_over_r*one_over_r*exp_minus_beta_r_over_r;
			access_pot(ii)[k] += qj*exp_minus_beta_r_over_r + m_times_rj*hh0;
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
			m_times_ri = x*mxi + y*myi + z*mzi;
			access_pot(jj)[j] += qi*exp_minus_beta_r_over_r - m_times_ri*hh0;
			#endif
			#else
			access_pot(ii)[k] += qj*exp_minus_beta_r_over_r;
			#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD)) 
			access_pot(jj)[j] += qi*exp_minus_beta_r_over_r;
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
			access_gradx(jj)[j] += x*hh1;
			access_grady(jj)[j]+= y*hh1;
			access_gradz(jj)[j]+= z*hh1;
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
			access_gradx(jj)[j] += x*hh2 + mxi*hh0;
			access_grady(jj)[j] += y*hh2 + myi*hh0;
			access_gradz(jj)[j] += z*hh2 + mzi*hh0;	
			#endif

		}
	}
}


