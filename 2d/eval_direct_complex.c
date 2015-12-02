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

#if (defined(C_CHARGE))
    #define FMM_KIND FMM_C_CHARGE	
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_c_charge_acc0
	#define eval_direct_periodic eval_direct_periodic_c_charge_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_c_charge_acc1
	#define eval_direct_periodic eval_direct_periodic_c_charge_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_c_charge_acc2
	#define eval_direct_periodic eval_direct_periodic_c_charge_acc2
    #endif	
#elif (defined(C_DIPOLE))
    #define FMM_KIND FMM_C_DIPOLE	
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_c_dipole_acc0
	#define eval_direct_periodic eval_direct_periodic_c_dipole_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_c_dipole_acc1
	#define eval_direct_periodic eval_direct_periodic_c_dipole_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_c_dipole_acc2
	#define eval_direct_periodic eval_direct_periodic_c_dipole_acc2
    #endif	
#elif (defined(C_CHARGE_DIPOLE))
    #define FMM_KIND FMM_C_CHARGE_DIPOLE	
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_c_charge_dipole_acc0
	#define eval_direct_periodic eval_direct_periodic_c_charge_dipole_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_c_charge_dipole_acc1
	#define eval_direct_periodic eval_direct_periodic_c_charge_dipole_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_c_charge_dipole_acc2
	#define eval_direct_periodic eval_direct_periodic_c_charge_dipole_acc2
    #endif	
#elif (defined(C_CHARGE_GRAD))
    #define FMM_KIND FMM_C_CHARGE_GRAD	
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_c_charge_grad_acc0
	#define eval_direct_periodic eval_direct_periodic_c_charge_grad_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_c_charge_grad_acc1
	#define eval_direct_periodic eval_direct_periodic_c_charge_grad_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_c_charge_grad_acc2
	#define eval_direct_periodic eval_direct_periodic_c_charge_grad_acc2
    #endif	
#elif (defined(C_DIPOLE_GRAD))
    #define FMM_KIND FMM_C_DIPOLE_GRAD	
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_c_dipole_grad_acc0
	#define eval_direct_periodic eval_direct_periodic_c_dipole_grad_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_c_dipole_grad_acc1
	#define eval_direct_periodic eval_direct_periodic_c_dipole_grad_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_c_dipole_grad_acc2
	#define eval_direct_periodic eval_direct_periodic_c_dipole_grad_acc2
    #endif	
#elif (defined(C_CHARGE_DIPOLE_GRAD))
    #define FMM_KIND FMM_C_CHARGE_DIPOLE_GRAD	
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_c_charge_dipole_grad_acc0
	#define eval_direct_periodic eval_direct_periodic_c_charge_dipole_grad_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_c_charge_dipole_grad_acc1
	#define eval_direct_periodic eval_direct_periodic_c_charge_dipole_grad_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_c_charge_dipole_grad_acc2
	#define eval_direct_periodic eval_direct_periodic_c_charge_dipole_grad_acc2
    #endif	
#elif (defined(ST_C_CHARGE))
    #define FMM_KIND FMM_ST_C_CHARGE	
    #define ST
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_ST_c_charge_acc0
	#define eval_direct_periodic eval_direct_periodic_ST_c_charge_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_ST_c_charge_acc1
	#define eval_direct_periodic eval_direct_periodic_ST_c_charge_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_ST_c_charge_acc2
	#define eval_direct_periodic eval_direct_periodic_ST_c_charge_acc2
    #endif	
#elif (defined(ST_C_DIPOLE))
    #define FMM_KIND FMM_ST_C_DIPOLE	
    #define ST
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_ST_c_dipole_acc0
	#define eval_direct_periodic eval_direct_periodic_ST_c_dipole_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_ST_c_dipole_acc1
	#define eval_direct_periodic eval_direct_periodic_ST_c_dipole_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_ST_c_dipole_acc2
	#define eval_direct_periodic eval_direct_periodic_ST_c_dipole_acc2
    #endif	
#elif (defined(ST_C_CHARGE_DIPOLE))
    #define FMM_KIND FMM_ST_C_CHARGE_DIPOLE	
    #define ST
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_ST_c_charge_dipole_acc0
	#define eval_direct_periodic eval_direct_periodic_ST_c_charge_dipole_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_ST_c_charge_dipole_acc1
	#define eval_direct_periodic eval_direct_periodic_ST_c_charge_dipole_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_ST_c_charge_dipole_acc2
	#define eval_direct_periodic eval_direct_periodic_ST_c_charge_dipole_acc2
    #endif	
#elif (defined(ST_C_CHARGE_GRAD))
    #define FMM_KIND FMM_ST_C_CHARGE_GRAD	
    #define ST
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_ST_c_charge_grad_acc0
	#define eval_direct_periodic eval_direct_periodic_ST_c_charge_grad_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_ST_c_charge_grad_acc1
	#define eval_direct_periodic eval_direct_periodic_ST_c_charge_grad_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_ST_c_charge_grad_acc2
	#define eval_direct_periodic eval_direct_periodic_ST_c_charge_grad_acc2
    #endif	
#elif (defined(ST_C_DIPOLE_GRAD))
    #define FMM_KIND FMM_ST_C_DIPOLE_GRAD	
    #define ST
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_ST_c_dipole_grad_acc0
	#define eval_direct_periodic eval_direct_periodic_ST_c_dipole_grad_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_ST_c_dipole_grad_acc1
	#define eval_direct_periodic eval_direct_periodic_ST_c_dipole_grad_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_ST_c_dipole_grad_acc2
	#define eval_direct_periodic eval_direct_periodic_ST_c_dipole_grad_acc2
    #endif	
#elif (defined(ST_C_CHARGE_DIPOLE_GRAD))
    #define FMM_KIND FMM_ST_C_CHARGE_DIPOLE_GRAD	
    #define ST
    #if (ACCURACY==0)	
	#define eval_direct eval_direct_ST_c_charge_dipole_grad_acc0
	#define eval_direct_periodic eval_direct_periodic_ST_c_charge_dipole_grad_acc0
    #elif (ACCURACY==1)	
	#define eval_direct eval_direct_ST_c_charge_dipole_grad_acc1
	#define eval_direct_periodic eval_direct_periodic_ST_c_charge_dipole_grad_acc1
    #elif (ACCURACY==2)	
	#define eval_direct eval_direct_ST_c_charge_dipole_grad_acc2
	#define eval_direct_periodic eval_direct_periodic_ST_c_charge_dipole_grad_acc2
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

	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	_FLOAT_ x, y, r2, log_r, atan2_yx; 
	#if(defined(access_qx))
        _FLOAT_ qxj, qyj;
        #endif
	#if(defined(access_potx))
        _FLOAT_ potxi, potyi;
        #endif
	#if(defined(access_gradx))
        _FLOAT_ gradxi, gradyi;
        _FLOAT_ x_over_r2, y_over_r2;
        #endif 
	#ifndef ST
	#if(defined(access_qx))
        _FLOAT_ qxi, qyi;
        #endif
	#if(defined(access_potx))
        _FLOAT_ potxj, potyj;
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
        _FLOAT_ mr_re_over_r2j;
        _FLOAT_ mr_im_over_r2j;
	#ifndef ST
	_FLOAT_ mxi, myi;
        _FLOAT_ mr_re_over_r2i;
        _FLOAT_ mr_im_over_r2i;
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
	        #if(defined(access_qx))
		qxi = access_qx(i);
		qyi = access_qy(i);
                #endif
		#if(defined(access_mx))
		mxi = access_mx(i);
		myi = access_my(i);
		#endif
		#endif

	        #if(defined(access_potx))
		potxi = access_potx(i);
		potyi = access_poty(i);
                #endif
	        #if(defined(access_gradx))
                gradxi = access_gradx(i);
                gradyi = access_grady(i);
                #endif
		
        	#ifndef ST 
		j00 = (i<j0 ? j0 : i+1);
		#endif

		for (j=j00; j<j1; j++) {
	                #if(defined(access_qx))
			qxj = access_qx(j);
			qyj = access_qy(j);
                        #endif
	                #if(defined(access_potx)&&(!defined(ST)))
			potxj = access_potx(j);
			potyj = access_poty(j);
                        #endif
			x = xi - access_x(j);
			y = yi - access_y(j);
			r2 = x*x + y*y;

	                #if(defined(access_potx)&&defined(access_qx))

			#if (EVAL_DIRECT_ACCURACY==0)
			log_r = 0.5*LOG0(r2);
                        atan2_yx = ATAN20(y,x);
			#elif (EVAL_DIRECT_ACCURACY==1)
			log_r = 0.5*LOG1(r2);
                        atan2_yx = ATAN21(y,x);
			#elif (EVAL_DIRECT_ACCURACY==2)
			log_r = 0.5*LOG2(r2);
                        atan2_yx = ATAN22(y,x);
			#endif

			potxi += qxj*log_r - qyj*atan2_yx;
			potyi += qyj*log_r + qxj*atan2_yx;
			#ifndef ST
			potxj += qxi*log_r - qyi*atan2_yx;
			potyj += qyi*log_r + qxi*atan2_yx;
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
	                #if(defined(access_qx))
                        x_over_r2 = x*one_over_r2;
                        y_over_r2 = y*one_over_r2;
                        gradxi += qxj*x_over_r2 + qyj*y_over_r2;
                        gradyi += qyj*x_over_r2 - qxj*y_over_r2;
			#ifndef ST
                        gradxj -= qxi*x_over_r2 + qyi*y_over_r2;
                        gradyj -= qyi*x_over_r2 - qxi*y_over_r2;
                        #endif
                        #endif
                        #endif


			#if(defined(access_mx))

                        mxj = access_mx(j);
                        myj = access_my(j);
			mr_re_over_r2j = (x*mxj + y*myj)*one_over_r2;
			mr_im_over_r2j = (x*myj - y*mxj)*one_over_r2; /* TODO: CHECK SIGN!!!! */
	                #if(defined(access_potx))
                        potxi += mr_re_over_r2j;
                        potyi += mr_im_over_r2j;
                        #endif
                        #if(defined(access_gradx)) 
                        gradxi -= 2.0*(x*mr_re_over_r2j+y*mr_im_over_r2j)*one_over_r2;
                        gradyi -= 2.0*(x*mr_im_over_r2j-y*mr_re_over_r2j)*one_over_r2;
                        #endif

			#ifndef ST
			mr_re_over_r2i = (x*mxi + y*myi)*one_over_r2;
			mr_im_over_r2i = (x*myi - y*mxi)*one_over_r2; /* TODO: CHECK SIGN!!!! */
	                #if(defined(access_potx))
                        potxj -= mr_re_over_r2i;
                        potyj -= mr_im_over_r2i;
                        #endif
                        #if(defined(access_gradx))  
                        gradxj -= 2.0*(x*mr_re_over_r2i+y*mr_im_over_r2i)*one_over_r2;
                        gradyj -= 2.0*(x*mr_im_over_r2i-y*mr_re_over_r2i)*one_over_r2;
                        #endif
			#endif

			#endif

			#ifndef ST
	                #if(defined(access_potx))
                        access_potx(j) = potxj;
                        access_poty(j) = potyj;
			#endif
                        #if(defined(access_gradx))
                        access_gradx(j) = gradxj;
                        access_grady(j) = gradyj;
                        #endif
                        #endif
		}

	        #if(defined(access_potx))
                access_potx(i) = potxi;
                access_poty(i) = potyi;
                #endif
                #if(defined(access_gradx))
                access_gradx(i) = gradxi;
                access_grady(i) = gradyi;
                #endif
	}
}




