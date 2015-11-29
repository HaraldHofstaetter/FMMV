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


//#include<assert.h>
//
//void extrinsic_correction(FmmvHandle *FMMV)
//{
////	assert(0);
//}	

#if (defined(STANDARD))
	#define FMM_KIND FMM_STANDARD	
	#define extrinsic_correction extrinsic_correction_standard
#elif (defined(GRAD))
	#define FMM_KIND FMM_GRAD	
	#define extrinsic_correction extrinsic_correction_grad
#elif (defined(DIPOLE))
	#define FMM_KIND FMM_DIPOLE	
	#define extrinsic_correction extrinsic_correction_dipole
#elif (defined(DIPOLE_GRAD))
	#define FMM_KIND FMM_DIPOLE_GRAD	
	#define extrinsic_correction extrinsic_correction_dipole_grad
#elif (defined(ST_STANDARD))
	#define FMM_KIND FMM_ST_STANDARD	
	#define extrinsic_correction extrinsic_correction_ST_standard
#elif (defined(ST_GRAD))
	#define FMM_KIND FMM_ST_GRAD	
	#define extrinsic_correction extrinsic_correction_ST_grad
#elif (defined(ST_DIPOLE))
	#define FMM_KIND FMM_ST_DIPOLE	
	#define extrinsic_correction extrinsic_correction_ST_dipole
#elif (defined(ST_DIPOLE_GRAD))
	#define FMM_KIND FMM_ST_DIPOLE_GRAD	
	#define extrinsic_correction extrinsic_correction_ST_dipole_grad
#endif

#include"fmmv_access.h"

void extrinsic_correction(FmmvHandle *FMMV)
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	_FLOAT_ /* dx, */ dy;
	_FLOAT_ Q;
	_FLOAT_  /* xi, */ yi, qi;
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	_FLOAT_ /* mxi, */ myi;
	#endif
	int i;
        _FLOAT_ PI = 3.14159265358979323846;
	
/*	dx = 0.0; */
	dy = 0.0;
	Q = 0.0;
	
	for (i=0; i<FMMV->NParticles; i++) {
		/* xi = access_x(i); */
		yi = access_y(i);
		qi = access_q(i);
		
		/*dx = dx + qi*xi; */
		dy = dy + qi*yi;

		Q = Q + qi*(/*xi*xi +*/ yi*yi);

		#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
		   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		/*mxi = access_mx(i); */
		myi = access_my(i);

		/* dx = dx + mxi; */
		dy = dy - myi; 
		
		Q = Q + 2*( /* mxi*xi */ -  myi*yi); 
		#endif
	}

	Q  *= -PI;
        /* dx *=2*PI; */
        dy *=2*PI;

	for (i=0; i<FMMV->NTargets; i++) {
		/* xi = access_tx(i); */
		yi = access_ty(i);
		
		access_pot(i) += /* dx*xi + */ dy*yi + Q;

		#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
		   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		/* access_gradx(i) += dx; */
		access_grady(i) += dy;
		#endif
	}
}	
		
	
 
