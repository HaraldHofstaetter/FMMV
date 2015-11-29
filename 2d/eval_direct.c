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

#undef ST

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
	#define ST
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
	#define ST
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
	#define ST
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
	#define ST
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
void eval_direct_yukawa_periodic(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy);
#else
void eval_direct_yukawa(FmmvHandle *FMMV, Box *target, Box *source);
#endif


#define EVAL_DIRECT_ACCURACY ACCURACY

#include"fmmv_access.h"


#ifdef PERIODIC
void eval_direct_periodic(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy)
#else
void eval_direct(FmmvHandle *FMMV, Box *target, Box *source)
#endif
{
    #ifdef ST
    if (!(isTarget(target)&&isSource(source))) return;
    #else
    if (!(target&&source)||(target->firstTarget>source->firstParticle)) return;
    #endif

    if (FMMV->beta!=0){
#ifdef PERIODIC
           eval_direct_yukawa_periodic(FMMV, target, source, dx, dy);
#else
           eval_direct_yukawa(FMMV, target, source);
#endif
           return;
    }
    else {

	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	_FLOAT_ x, y, r2, log_r; 
	#if(defined(access_q))
        _FLOAT_ qj;
        #endif
	#if(defined(access_pot))
        _FLOAT_ poti;
        #endif
	#if(defined(access_gradx))
        _FLOAT_ gradxi, gradyi;
        _FLOAT_ x_over_r2, y_over_r2;
        #endif 
	#ifndef ST
	#if(defined(access_q))
        _FLOAT_ qi;
        #endif
	#if(defined(access_pot))
        _FLOAT_ potj;
        #endif
	#if(defined(access_gradx))
        _FLOAT_ gradxj, gradyj;
        #endif 
	#endif
	#if(defined(access_mx)||defined(access_gradx))
	_FLOAT_ one_over_r2;
	#endif
	_FLOAT_ xi,yi;
	#if(defined(access_mx))
	_FLOAT_ mxj, myj;
        _FLOAT_ m_x_over_r2j;
	#ifndef ST
	_FLOAT_ mxi, myi;
        _FLOAT_ m_x_over_r2i;
	#endif
	#endif

	int i,j,ni,nj,i0,j0,j00,j1;

		
        #ifdef ST
        if (!(isTarget(target)&&isSource(source))) return;
        #else
        if (!(target&&source)||(target->firstTarget>source->firstParticle)) return;
        #endif

	i0 = target->firstTarget;
	ni = target->noOfTargets;
	j0 = source->firstParticle;
	nj = source->noOfParticles;

        #ifdef ST
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

	for (i=i0; i<i0+ni; i++) {
		#ifdef PERIODIC
		xi = access_tx(i) - dx;
		yi = access_ty(i) - dy;
		#else
		xi = access_tx(i);
		yi = access_ty(i);
		#endif

        	#ifndef ST
	        #if(defined(access_q))
		qi = access_q(i);
                #endif
		#if(defined(access_mx))
		mxi = access_mx(i);
		myi = access_my(i);
		#endif
		#endif

	        #if(defined(access_pot))
		poti = access_pot(i);
                #endif
	        #if(defined(access_gradx))
                gradxi = access_gradx(i);
                gradyi = access_grady(i);
                #endif
		
        	#ifndef ST 
		j00 = (i<j0 ? j0 : i+1);
		#endif

		for (j=j00; j<j1; j++) {
	                #if(defined(access_q))
			qj = access_q(j);
                        #endif
	                #if(defined(access_pot)&&(!defined(ST)))
			potj = access_pot(j);
                        #endif
			x = xi - access_x(j);
			y = yi - access_y(j);
			r2 = x*x+y*y;

	                #if(defined(access_pot)&&defined(access_q))

			#if (EVAL_DIRECT_ACCURACY==0)
			log_r = 0.5*LOG0(r2);
			#elif (EVAL_DIRECT_ACCURACY==1)
			log_r = 0.5*LOG1(r2);
			#elif (EVAL_DIRECT_ACCURACY==2)
			log_r = 0.5*LOG2(r2);
			#endif

			poti += qj*log_r;
			#ifndef ST
			potj += qi*log_r;
			#endif

                        #endif

			#if(defined(access_mx)||defined(access_gradx))

			#if (EVAL_DIRECT_ACCURACY==0)
			one_over_r2 = RECIP0(r2); 
			#elif (EVAL_DIRECT_ACCURACY==1)
			one_over_r2 = RECIP1(r2); 
			#elif (EVAL_DIRECT_ACCURACY==2)
			one_over_r2 = RECIP2(r2); 
			#endif

                        #endif
	                
                        #if(defined(access_gradx))
			#ifndef ST
                        gradxj = access_gradx(j);
                        gradyj = access_grady(j);
                        #endif
	                #if(defined(access_q))
                        x_over_r2 = x*one_over_r2;
                        y_over_r2 = y*one_over_r2;
                        gradxi += qj*x_over_r2;
                        gradyi += qj*y_over_r2;
			#ifndef ST
                        gradxj -= qi*x_over_r2;
                        gradyj -= qi*y_over_r2;
                        #endif
                        #endif
                        #endif


			#if(defined(access_mx))

                        mxj = access_mx(j);
                        myj = access_my(j);
			m_x_over_r2j = (x*mxj + y*myj)*one_over_r2;
	                #if(defined(access_pot))
                        poti += m_x_over_r2j;
                        #endif
                        #if(defined(access_gradx))
                        gradxi += (mxj - 2.0*x*m_x_over_r2j)*one_over_r2;
                        gradyi += (myj - 2.0*y*m_x_over_r2j)*one_over_r2;
                        #endif

			#ifndef ST
			m_x_over_r2i = (x*mxi + y*myi)*one_over_r2;
	                #if(defined(access_pot))
                        potj -= m_x_over_r2i;
                        #endif
                        #if(defined(access_gradx)) 
                        gradxj += (mxi - 2.0*x*m_x_over_r2i)*one_over_r2;
                        gradyj += (myi - 2.0*y*m_x_over_r2i)*one_over_r2;
                        #endif
			#endif

			#endif

			#ifndef ST
	                #if(defined(access_pot))
                        access_pot(j) = potj;
			#endif
                        #if(defined(access_gradx))
                        access_gradx(j) = gradxj;
                        access_grady(j) = gradyj;
                        #endif
                        #endif
		}

	        #if(defined(access_pot))
                access_pot(i) = poti;
                #endif
                #if(defined(access_gradx))
                access_gradx(i) = gradxi;
                access_grady(i) = gradyi;
                #endif
	}
    }
}




#ifdef PERIODIC
void eval_direct_yukawa_periodic(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy)
#else
void eval_direct_yukawa(FmmvHandle *FMMV, Box *target, Box *source)
#endif
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)
        
        _FLOAT_ beta = FMMV->beta;
        
	_FLOAT_ x, y, r, k0_beta_r;
	#if(defined(access_q))
        _FLOAT_ qj;
        _FLOAT_ k0_correction = FMMV->k0_correction;
        #endif
	#if(defined(access_pot))
        _FLOAT_ poti;
        #endif
	#if(defined(access_gradx))
        _FLOAT_ gradxi, gradyi;
	#if(defined(access_q))
        _FLOAT_ qj_beta_k1_over_r;        
        #endif
        #endif 
	#ifndef ST
	#if(defined(access_q))
        _FLOAT_ qi;
        #endif
	#if(defined(access_pot))
        _FLOAT_ potj;
        #endif
	#if(defined(access_gradx))
        _FLOAT_ gradxj, gradyj;
	#if(defined(access_q))
        _FLOAT_ qi_beta_k1_over_r;        
        #endif
        #endif 
	#endif
	#if(defined(access_mx)||defined(access_gradx))
        _FLOAT_ one_over_r, beta_k1_over_r;
	#endif
	_FLOAT_ xi,yi;
        
	#if(defined(access_mx))
	_FLOAT_ mxj, myj, mj_times_r;        
	#ifndef ST
	_FLOAT_ mxi, myi, mi_times_r;
	#endif
        #if(defined(access_gradx))
        _FLOAT_ one_over_r2, beta_k1_over_r3, beta2_k0_over_r2;
        _FLOAT_ xmyj_m_ymxj;
	#ifndef ST
        _FLOAT_ xmyi_m_ymxi;
	#endif
	#endif
	#endif

	int i,j,ni,nj,i0,j0,j00,j1;
		
        #ifdef ST
        if (!(isTarget(target)&&isSource(source))) return;
        #else
        if (!(target&&source)||(target->firstTarget>source->firstParticle)) return;
        #endif

	i0 = target->firstTarget;
	ni = target->noOfTargets;
	j0 = source->firstParticle;
	nj = source->noOfParticles;

        #ifdef ST
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

	for (i=i0; i<i0+ni; i++) {
		#ifdef PERIODIC
		xi = access_tx(i) - dx;
		yi = access_ty(i) - dy;
		#else
		xi = access_tx(i);
		yi = access_ty(i);
		#endif

        	#ifndef ST
	        #if(defined(access_q))
		qi = access_q(i);
                #endif
		#if(defined(access_mx))
		mxi = access_mx(i);
		myi = access_my(i);
		#endif
		#endif

	        #if(defined(access_pot))
		poti = access_pot(i);
                #endif
	        #if(defined(access_gradx))
                gradxi = access_gradx(i);
                gradyi = access_grady(i);
                #endif
		
        	#ifndef ST 
		j00 = (i<j0 ? j0 : i+1);
		#endif

		for (j=j00; j<j1; j++) {
	                #if(defined(access_q))
			qj = access_q(j);
                        #endif
	                #if(defined(access_pot)&&(!defined(ST)))
			potj = access_pot(j);
                        #endif
			x = xi - access_x(j);
			y = yi - access_y(j);

	                r = sqrt(x*x + y*y);

	                #if(defined(access_pot)&&defined(access_q))

			#if (EVAL_DIRECT_ACCURACY==0)
			k0_beta_r = K00(beta*r);
			#elif (EVAL_DIRECT_ACCURACY==1)
			k0_beta_r = K01(beta*r);
			#elif (EVAL_DIRECT_ACCURACY==2)
			k0_beta_r = K02(beta*r);
			#endif

			poti += qj*(k0_beta_r + k0_correction);
			#ifndef ST
			potj += qi*(k0_beta_r + k0_correction);
			#endif

                        #endif

			#if(defined(access_mx)||defined(access_gradx))

			#if (EVAL_DIRECT_ACCURACY==0)
                        one_over_r = RECIP0(r);
                        beta_k1_over_r = beta*K10(beta*r)*one_over_r;
			#elif (EVAL_DIRECT_ACCURACY==1)
                        one_over_r = RECIP1(r);
                        beta_k1_over_r = beta*K11(beta*r)*one_over_r;
			#elif (EVAL_DIRECT_ACCURACY==2)
                        one_over_r = RECIP2(r);
                        beta_k1_over_r = beta*K12(beta*r)*one_over_r;
			#endif

                        #endif

                        #if(defined(access_gradx))
			#ifndef ST
                        gradxj = access_gradx(j);
                        gradyj = access_grady(j);
                        #endif
	                #if(defined(access_q))
                        qj_beta_k1_over_r = qj*beta_k1_over_r;
                        gradxi -= x*qj_beta_k1_over_r;
                        gradyi -= y*qj_beta_k1_over_r;
			#ifndef ST
                        qi_beta_k1_over_r = qi*beta_k1_over_r;
                        gradxj += x*qi_beta_k1_over_r;
                        gradyj += y*qi_beta_k1_over_r;
                        #endif
                        #endif
                        #endif

			#if(defined(access_mx))

                        mxj = access_mx(j);
                        myj = access_my(j);
                        mj_times_r = x*mxj + y*myj;
	                #if(defined(access_pot))
                        poti -= mj_times_r*beta_k1_over_r;
                        #endif
                        #if(defined(access_gradx))
                        one_over_r2 = one_over_r*one_over_r;
                        beta_k1_over_r3 = beta_k1_over_r*one_over_r2;
                        beta2_k0_over_r2 = beta*beta*k0_beta_r*one_over_r2;
                        xmyj_m_ymxj =  myj*x - mxj*y;
                        gradxi += x*mj_times_r*beta2_k0_over_r2 + (x*mj_times_r+y*xmyj_m_ymxj)*beta_k1_over_r3; 
                        gradyi += y*mj_times_r*beta2_k0_over_r2 + (y*mj_times_r-x*xmyj_m_ymxj)*beta_k1_over_r3; 
                        #endif

			#ifndef ST
                        mj_times_r = x*mxi + y*myi;
	                #if(defined(access_pot))
                        potj += mi_times_r*beta_k1_over_r;
                        #endif
                        #if(defined(access_gradx)) 
                        xmyi_m_ymxi =  myi*x - mxi*y;
                        gradxj -= x*mi_times_r*beta2_k0_over_r2 + (x*mi_times_r+y*xmyi_m_ymxi)*beta_k1_over_r3; 
                        gradyj -= y*mi_times_r*beta2_k0_over_r2 + (y*mi_times_r-x*xmyi_m_ymxi)*beta_k1_over_r3; 

                        #endif
			#endif

			#endif

			#ifndef ST
	                #if(defined(access_pot))
                        access_pot(j) = potj;
			#endif
                        #if(defined(access_gradx))
                        access_gradx(j) = gradxj;
                        access_grady(j) = gradyj;
                        #endif
                        #endif
		}

	        #if(defined(access_pot))
                access_pot(i) = poti;
                #endif
                #if(defined(access_gradx))
                access_gradx(i) = gradxi;
                access_grady(i) = gradyi;
                #endif
	}
}

	
	

	

