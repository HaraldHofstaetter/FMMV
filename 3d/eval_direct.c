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
void eval_direct_periodic(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz)
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
           eval_direct_yukawa_periodic(FMMV, target, source, dx, dy, dz);
#else
           eval_direct_yukawa(FMMV, target, source);
#endif
           return;
    }
    else {
	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	_FLOAT_ x, y, z, one_over_r;
	_FLOAT_ xi,yi,zi, qj;
	#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
	   ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
	_FLOAT_ qi;
	#endif
	int i,j,ni,nj,i0,j0,j00,j1;
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	_FLOAT_ one_o_r_3;
	#endif
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	_FLOAT_ qj_o_r_3;
	#endif
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
	_FLOAT_ qi_o_r_3;
	#endif
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	_FLOAT_ mxj, myj, mzj;
	_FLOAT_ m_times_rj;	
	#endif
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
	_FLOAT_ mxi, myi, mzi;
	_FLOAT_ m_times_ri;	
	#endif
	#if ((FMM_KIND==FMM_DIPOLE_GRAD) \
	   ||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	_FLOAT_ one_o_r_5;
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

	for (i=i0; i<i0+ni; i++) {
		#ifdef PERIODIC
		xi = access_tx(i) - dx;
		yi = access_ty(i) - dy;
		zi = access_tz(i) - dz;
		#else
		xi = access_tx(i);
		yi = access_ty(i);
		zi = access_tz(i);
		#endif
		#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
		   ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
		qi = access_q(i);
		#endif
		#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
		mxi = access_mx(i);
		myi = access_my(i);
		mzi = access_mz(i);		
		#endif
		
		#if ((FMM_KIND==FMM_ST_STANDARD)||(FMM_KIND==FMM_ST_GRAD) \
		   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		#else
		j00 = (i<j0 ? j0 : i+1);
		#endif

		for (j=j00; j<j1; j++) {
			x = xi - access_x(j);
			y = yi - access_y(j);
			z = zi - access_z(j);
			#if (EVAL_DIRECT_ACCURACY==0)
			one_over_r = RECIP_SQRT0(x*x + y*y + z*z); 
			#elif (EVAL_DIRECT_ACCURACY==1)
			one_over_r = RECIP_SQRT1(x*x + y*y + z*z); 
			#elif (EVAL_DIRECT_ACCURACY==2)
			one_over_r = RECIP_SQRT2(x*x + y*y + z*z); 
			#endif

			qj = access_q(j);
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
			  ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			mxj = access_mx(j);
			myj = access_my(j);
			mzj = access_mz(j);
			one_o_r_3 = one_over_r*one_over_r*one_over_r;
			m_times_rj = x*mxj + y*myj + z*mzj;
			access_pot(i) += qj*one_over_r + m_times_rj*one_o_r_3;
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
			m_times_ri = x*mxi + y*myi + z*mzi;
			access_pot(j) += qi*one_over_r - m_times_ri*one_o_r_3;
			#endif
			#else
			access_pot(i) += qj*one_over_r;
			#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD)) 
			access_pot(j) += qi*one_over_r;
			#endif
			#endif
			
			#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_ST_GRAD))
			one_o_r_3 = one_over_r*one_over_r*one_over_r;
			qj_o_r_3 = qj*one_o_r_3;
			access_gradx(i) -= x*qj_o_r_3;
			access_grady(i) -= y*qj_o_r_3;
			access_gradz(i) -= z*qj_o_r_3;
			#endif
			#if (FMM_KIND==FMM_GRAD)
			qi_o_r_3 = qi*one_o_r_3;
			access_gradx(j) += x*qi_o_r_3;
			access_grady(j) += y*qi_o_r_3;
			access_gradz(j) += z*qi_o_r_3;
			#endif
			#if ((FMM_KIND==FMM_DIPOLE_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			one_o_r_5 = one_o_r_3*one_over_r*one_over_r;
			qj_o_r_3 = qj*one_o_r_3 + 3.0*m_times_rj*one_o_r_5;
			access_gradx(i) -= x*qj_o_r_3 - mxj*one_o_r_3;
			access_grady(i) -= y*qj_o_r_3 - myj*one_o_r_3;
			access_gradz(i) -= z*qj_o_r_3 - mzj*one_o_r_3;
			#endif
			#if (FMM_KIND==FMM_DIPOLE_GRAD)
			qi_o_r_3 = qi*one_o_r_3 - 3.0*m_times_ri*one_o_r_5;	
			access_gradx(j) += x*qi_o_r_3 + mxi*one_o_r_3;
			access_grady(j) += y*qi_o_r_3 + myi*one_o_r_3;
			access_gradz(j) += z*qi_o_r_3 + mzi*one_o_r_3;	
			#endif
		}
	}
    }
}




#ifdef PERIODIC
void eval_direct_yukawa_periodic(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz)
#else
void eval_direct_yukawa(FmmvHandle *FMMV, Box *target, Box *source)
#endif
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)
        
        _FLOAT_ beta = FMMV->beta;

        _FLOAT_ beta_r, exp_minus_beta_r_over_r;

	_FLOAT_ x, y, z, one_over_r, r2;
	_FLOAT_ xi,yi,zi, qj;
	#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
	   ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
	_FLOAT_ qi;
	#endif
	int i,j,ni,nj,i0,j0,j00,j1;
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	_FLOAT_ hh0;
	#endif
	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	_FLOAT_ hh1;
	#endif
	#if ((FMM_KIND==FMM_DIPOLE_GRAD) \
	   ||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	_FLOAT_  hh2;
	_FLOAT_ beta2 = beta*beta;	
	#endif
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
 	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	_FLOAT_ mxj, myj, mzj;
	_FLOAT_ m_times_rj;	
	#endif
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
	_FLOAT_ mxi, myi, mzi;
	_FLOAT_ m_times_ri;	
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

	for (i=i0; i<i0+ni; i++) {
		#ifdef PERIODIC
		xi = access_tx(i) - dx;
		yi = access_ty(i) - dy;
		zi = access_tz(i) - dz;
		#else
		xi = access_tx(i);
		yi = access_ty(i);
		zi = access_tz(i);
		#endif
		#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD) \
		   ||(FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD))
		qi = access_q(i);
		#endif
		#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
		mxi = access_mx(i);
		myi = access_my(i);
		mzi = access_mz(i);		
		#endif
		
		#if ((FMM_KIND==FMM_ST_STANDARD)||(FMM_KIND==FMM_ST_GRAD) \
		   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		#else
		j00 = (i<j0 ? j0 : i+1);
		#endif

		for (j=j00; j<j1; j++) {
			x = xi - access_x(j);
			y = yi - access_y(j);
			z = zi - access_z(j);
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


			qj = access_q(j);
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
			  ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			mxj = access_mx(j);
			myj = access_my(j);
			mzj = access_mz(j);
			m_times_rj = x*mxj + y*myj + z*mzj;
                        hh0 = (1.0+beta_r)*one_over_r*one_over_r*exp_minus_beta_r_over_r;
			access_pot(i) += qj*exp_minus_beta_r_over_r + m_times_rj*hh0;
			#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD)) 
			m_times_ri = x*mxi + y*myi + z*mzi;
			access_pot(j) += qi*exp_minus_beta_r_over_r - m_times_ri*hh0;
			#endif
			#else
			access_pot(i) += qj*exp_minus_beta_r_over_r;
			#if ((FMM_KIND==FMM_STANDARD)||(FMM_KIND==FMM_GRAD)) 
			access_pot(j) += qi*exp_minus_beta_r_over_r;
			#endif
			#endif
			
			#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_ST_GRAD))
                        hh0 = (1.0+beta_r)*one_over_r*one_over_r*exp_minus_beta_r_over_r;
                        hh1 = qj*hh0;
			access_gradx(i) -= x*hh1;
			access_grady(i) -= y*hh1;
			access_gradz(i) -= z*hh1;
			#endif
			#if (FMM_KIND==FMM_GRAD)
                        hh1 = qi*hh0;
			access_gradx(j) += x*hh1;
			access_grady(j) += y*hh1;
			access_gradz(j) += z*hh1;
			#endif
			#if ((FMM_KIND==FMM_DIPOLE_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
                        hh1 = one_over_r*one_over_r*m_times_rj;
                        hh2 = (qj + 3.0*hh1)*hh0 + beta2*hh1*exp_minus_beta_r_over_r;
			access_gradx(i) -= x*hh2 - mxj*hh0;
			access_grady(i) -= y*hh2 - myj*hh0;
			access_gradz(i) -= z*hh2 - mzj*hh0;
			#endif
			#if (FMM_KIND==FMM_DIPOLE_GRAD)
                        hh1 = one_over_r*one_over_r*m_times_ri;
                        hh2 = (qi - 3.0*hh1)*hh0 - beta2*hh1*exp_minus_beta_r_over_r;
			access_gradx(j) += x*hh2 + mxi*hh0;
			access_grady(j) += y*hh2 + myi*hh0;
			access_gradz(j) += z*hh2 + mzi*hh0;	
			#endif
		}
	}
}

