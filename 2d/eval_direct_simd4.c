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
   #include"simd4s.h"
#else
   #include"simd4d.h"
#endif
#include"math.h"

#undef ST

#if (defined(STANDARD))
        #define FMM_KIND FMM_STANDARD	
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_standard_acc0
	#define eval_direct_periodic eval_direct_periodic_standard_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_standard_acc1
	#define eval_direct_periodic eval_direct_periodic_standard_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_standard_acc2
	#define eval_direct_periodic eval_direct_periodic_standard_acc2
    #endif	
#elif (defined(GRAD))
        #define FMM_KIND FMM_GRAD	
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_grad_acc0
	#define eval_direct_periodic eval_direct_periodic_grad_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_grad_acc1
	#define eval_direct_periodic eval_direct_periodic_grad_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_grad_acc2
	#define eval_direct_periodic eval_direct_periodic_grad_acc2
    #endif	
#elif (defined(DIPOLE))
        #define FMM_KIND FMM_DIPOLE	
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_dipole_acc0
	#define eval_direct_periodic eval_direct_periodic_dipole_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_dipole_acc1
	#define eval_direct_periodic eval_direct_periodic_dipole_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_dipole_acc2
	#define eval_direct_periodic eval_direct_periodic_dipole_acc2
    #endif	
#elif (defined(DIPOLE_GRAD))
        #define FMM_KIND FMM_DIPOLE_GRAD	
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_dipole_grad_acc0
	#define eval_direct_periodic eval_direct_periodic_dipole_grad_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_dipole_grad_acc1
	#define eval_direct_periodic eval_direct_periodic_dipole_grad_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_dipole_grad_acc2
	#define eval_direct_periodic eval_direct_periodic_dipole_grad_acc2
    #endif	
#elif (defined(ST_STANDARD))
        #define FMM_KIND FMM_ST_STANDARD	
	#define ST
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_ST_standard_acc0
	#define eval_direct_periodic eval_direct_periodic_ST_standard_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_ST_standard_acc1
	#define eval_direct_periodic eval_direct_periodic_ST_standard_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_ST_standard_acc2
	#define eval_direct_periodic eval_direct_periodic_ST_standard_acc2
    #endif	
#elif (defined(ST_GRAD))
        #define FMM_KIND FMM_ST_GRAD	
	#define ST
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_ST_grad_acc0
	#define eval_direct_periodic eval_direct_periodic_ST_grad_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_ST_grad_acc1
	#define eval_direct_periodic eval_direct_periodic_ST_grad_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_ST_grad_acc2
	#define eval_direct_periodic eval_direct_periodic_ST_grad_acc2
    #endif	
#elif (defined(ST_DIPOLE))
        #define FMM_KIND FMM_ST_DIPOLE	
	#define ST
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_ST_dipole_acc0
	#define eval_direct_periodic eval_direct_periodic_ST_dipole_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_ST_dipole_acc1
	#define eval_direct_periodic eval_direct_periodic_ST_dipole_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_ST_dipole_acc2
	#define eval_direct_periodic eval_direct_periodic_ST_dipole_acc2
    #endif	
#elif (defined(ST_DIPOLE_GRAD))
        #define FMM_KIND FMM_ST_DIPOLE_GRAD	
	#define ST
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_ST_dipole_grad_acc0
	#define eval_direct_periodic eval_direct_periodic_ST_dipole_grad_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_ST_dipole_grad_acc1
	#define eval_direct_periodic eval_direct_periodic_ST_dipole_grad_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_ST_dipole_grad_acc2
	#define eval_direct_periodic eval_direct_periodic_ST_dipole_grad_acc2
    #endif	
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
    {
	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	int i,j,ni,nj,i0,j0,j00,j1,k;
	int ii, jj, jj0, jj1, jjj;

	V4_TYPE x, y, r2, log_r; 
	#if(defined(access_q))
        V4_TYPE qj;
        #endif
	#if(defined(access_pot))
        V4_TYPE poti;
        #endif
	#if(defined(access_gradx))
        V4_TYPE gradxi, gradyi;
        V4_TYPE x_over_r2, y_over_r2;
        #endif 
	#ifndef ST
	#if(defined(access_q))
        V4_TYPE qi;
        #endif
	#if(defined(access_pot))
        V4_TYPE potj;
        #endif
	#if(defined(access_gradx))
        V4_TYPE gradxj, gradyj;
        #endif 
	#endif
	#if(defined(access_mx)||defined(access_gradx))
	V4_TYPE one_over_r2;
	#endif
	V4_TYPE xi,yi;
	#if(defined(access_mx))
	V4_TYPE mxj, myj;
        V4_TYPE m_x_over_r2j;
	#ifndef ST
	V4_TYPE mxi, myi;
        V4_TYPE m_x_over_r2i;
	#endif
	#endif

#ifdef PERIODIC
	V4_TYPE dx, dy;

	dx = V4_SET1(dx0);
	dy = V4_SET1(dy0);
#endif	

	i0 = target->firstTarget;
	ni = target->noOfTargets;
	j0 = source->firstParticle;
	nj = source->noOfParticles;

        #ifdef ST
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
		#else
		xi = V4_LOAD1(access_tx(ii)+k);
		yi = V4_LOAD1(access_ty(ii)+k);
		#endif

        	#ifndef ST
		#if(defined(access_q))
		qi = V4_LOAD1(access_q(ii)+k);
		#endif
		#if(defined(access_mx))
		mxi = V4_LOAD1(access_mx(ii)+k);
		myi = V4_LOAD1(access_my(ii)+k);
		#endif			
		#endif
		
		#if(defined(access_pot))
		poti = V4_LOAD0(access_pot(ii)+k);
		#endif
		#if(defined(access_gradx))
		gradxi = V4_LOAD0(access_gradx(ii)+k);
		gradyi = V4_LOAD0(access_grady(ii)+k);
		#endif

	        #ifndef ST 
		j00 = (i<j0 ? j0 : i+1);
		#endif
		jj0 = j00>>2;
		if (j00&3) jj0++;
		
		for (jj=jj0; jj<jj1; jj++) {
		        #if(defined(access_q))
			qj = V4_LOAD(access_q(jj));
		        #endif
                        #if(defined(access_pot)&&(!defined(ST)))
		        potj = V4_LOAD0(access_pot(jj)+k);
		        #endif
			x = V4_SUB(xi, V4_LOAD(access_x(jj)));
			y = V4_SUB(yi, V4_LOAD(access_y(jj)));

                        r2 = V4_FMA(x, x, V4_MUL(y, y))
	                #if(defined(access_pot)&&defined(access_q))
			#if (EVAL_DIRECT_ACCURACY==0)
			log_r = V4_MUL(0.5, V4_LOG0(r2));
			#elif (EVAL_DIRECT_ACCURACY==1)
			log_r = V4_MUL(0.5, V4_LOG1(r2));
			#elif (EVAL_DIRECT_ACCURACY==2)
			log_r = V4_MUL(0.5, V4_LOG2(r2));
			#endif

                        poti = V4_ADD(poti, V4_MUL(qj, log_r));
			#ifndef ST
                        potj = V4_ADD(potj, V4_MUL(qi, log_r));
			#endif
                        #endif

			#if(defined(access_mx)||defined(access_gradx))

			#if (EVAL_DIRECT_ACCURACY==0)
			one_over_r2 = V4_RECIP0(r2); 
			#elif (EVAL_DIRECT_ACCURACY==1)
			one_over_r2 = V4_RECIP1(r2); 
			#elif (EVAL_DIRECT_ACCURACY==2)
			one_over_r2 = V4_RECIP2(r2); 
			#endif

                        #endif
	                
                        #if(defined(access_gradx))
			#ifndef ST
                        gradxj = access_gradx(jj);
                        gradyj = access_grady(jj);
                        #endif
	                #if(defined(access_q))
                        x_over_r2 = V4_MUL(x, one_over_r2);
                        y_over_r2 = V4_MUL(y, one_over_r2);
                        gradxi = V4_FMA(qj, x_over_r2, gradxi);
                        gradyi = V4_FMA(qj, y_over_r2, gradyi);
			#ifndef ST
                        gradxj = V4_FNMS(qi, x_over_r2, gradxj);
                        gradyj = V4_FNMS(qi, x_over_r2, gradyj);
                        #endif
                        #endif
                        #endif


			#if(defined(access_mx))

			mxj = V4_LOAD(access_mx(jj));			
			myj = V4_LOAD(access_my(jj));	

			m_x_over_r2j = V4_MUL((V4_FMA(x,mxj), V4_MUL(y,myj)),one_over_r2);
	                #if(defined(access_pot))
                        poti = V4_ADD(poti, m_x_over_r2j);
                        #endif
                        #if(defined(access_gradx))
                        m_x_over_r2j = V4_MUL(V4_SET1(2.0), m_x_over_r2j);
                        gradxi = V4_FMA(V4_FNMS(x, m_x_over_r2j, mxj),one_over_r2, gradxi)
                        gradyi = V4_FMA(V4_FNMS(y, m_x_over_r2j, myj),one_over_r2, gradyi)
                        #endif

			#ifndef ST
			m_x_over_r2j = V4_MUL((V4_FMA(x,mxi), V4_MUL(y,myi)),one_over_r2);
	                #if(defined(access_pot))
                        poti = V4_SUB(poti, m_x_over_r2j);
                        #endif
                        #if(defined(access_gradx)) 
                        m_x_over_r2j = V4_MUL(V4_SET1(2.0), m_x_over_r2i);
                        gradxj = V4_FMA(V4_FNMS(x, m_x_over_r2i, mxi),one_over_r2, gradxj)
                        gradyj = V4_FMA(V4_FNMS(y, m_x_over_r2i, myi),one_over_r2, gradyj)
                        #endif
			#endif

			#endif

			#ifndef ST
	                #if(defined(access_pot))
                        V4_STORE(access_pot(jj), potj);
			#endif
                        #if(defined(access_gradx))
                        V4_STORE(access_gradx(jj), gradxj);
                        V4_STORE(access_grady(jj), gradyj);
                        #endif
                        #endif
                }

		V4_STORE0(access_pot(ii)+k, V4_HORIZADD1(poti));
		#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
		   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		V4_STORE0(access_gradx(ii)+k, V4_HORIZADD1(gradxi));
		V4_STORE0(access_grady(ii)+k, V4_HORIZADD1(gradyi));
		#endif

	    if ((j1>j00)&&((j00&3)||(j1&3))) {
	        V4_BASETYPE x, y, r2, log_r; 
	        #if(defined(access_q))
                V4_BASETYPE qj;
                #endif
	        #if(defined(access_pot))
                V4_BASETYPE poti;
                #endif
	        #if(defined(access_gradx))
                _FLOAT_ gradxi, gradyi;
                V4_BASETYPE x_over_r2, y_over_r2;
                #endif 
	        #ifndef ST
	        #if(defined(access_q))
                V4_BASETYPE qi;
                #endif
	        #if(defined(access_pot))
                V4_BASETYPE potj;
                #endif
	        #if(defined(access_gradx))
                V4_BASETYPE gradxj, gradyj;
                #endif 
	        #endif
	        #if(defined(access_mx)||defined(access_gradx))
	        V4_BASETYPE one_over_r2;
	        #endif
	        V4_BASETYPE xi,yi;
	        #if(defined(access_mx))
	        V4_BASETYPE mxj, myj;
                V4_BASETYPE m_x_over_r2j;
	        #ifndef ST
	        V4_BASETYPE mxi, myi;
                V4_BASETYPE m_x_over_r2i;
	        #endif
	        #endif
           
		#ifdef PERIODIC
		xi = access_tx(ii)[k] - dx0;
		yi = access_ty(ii)[k] - dy0;
		#else
		xi = access_tx(ii)[k];
		yi = access_ty(ii)[k];
		#endif

                #ifndef ST
                #if(defined(access_q))
		qi = access_q(ii)[k];
		#endif
                #if(defined(access_mx))
		mxi = access_mx(ii)[k];
		myi = access_my(ii)[k];
		#endif
		#endif

	        #if(defined(access_pot))
		poti = access_pot(ii)[k];
                #endif
	        #if(defined(access_gradx))
                gradxi = access_gradx(ii)[k];
                gradyi = access_grady(ii)[k];
                #endif
                
		
		if (j00&3) {
		    jj0--;
		    jjj = (jj0==jj1 ? j1&3 : 4);
		    for (j=j00&3; j<jjj; j++) {
	                #if(defined(access_q))
			qj = access_q(jj0)[j];
                        #endif
	                #if(defined(access_pot)&&(!defined(ST)))
			potj = access_pot(jj0)[j];
                        #endif
			x = xi - access_x(jj0)[j];
			y = yi - access_y(jj0)[j];
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
                        gradxj = access_gradx(jj0)[j];
                        gradyj = access_grady(jj0)[j];
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

                        mxj = access_mx(jj0)[j];
                        myj = access_my(jj0)[j];
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
                        access_pot(jj0)[j] = potj;
			#endif
                        #if(defined(access_gradx))
                        access_gradx(jj0)[j] = gradxj;
                        access_grady(jj0)[j] = gradyj;
                        #endif
                        #endif

	                #if(defined(access_pot))
                        access_pot(ii)[k] = poti;
                        #endif
                        #if(defined(access_gradx))
                        access_gradx(ii)[k] = gradxi;
                        access_grady(ii)[k] = gradyi;
                        #endif
		    }
		    if (jj0==jj1) continue;

	    	}

    		for (j=0; j<(j1&3); j++) {
	                #if(defined(access_q))
			qj = access_q(jj1)[j];
                        #endif
	                #if(defined(access_pot)&&(!defined(ST)))
			potj = access_pot(jj1)[j];
                        #endif
			x = xi - access_x(jj1)[j];
			y = yi - access_y(jj1)[j];
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
                        gradxj = access_gradx(jj1)[j];
                        gradyj = access_grady(jj1)[j];
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

                        mxj = access_mx(jj1)[j];
                        myj = access_my(jj1)[j];
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
                        access_pot(jj1)[j] = potj;
			#endif
                        #if(defined(access_gradx))
                        access_gradx(jj1)[j] = gradxj;
                        access_grady(jj1)[j] = gradyj;
                        #endif
                        #endif

	                #if(defined(access_pot))
                        access_pot(ii)[k] = poti;
                        #endif
                        #if(defined(access_gradx))
                        access_gradx(ii)[k] = gradxi;
                        access_grady(ii)[k] = gradyi;
                        #endif

	       }
	    }
	}
    }
}



