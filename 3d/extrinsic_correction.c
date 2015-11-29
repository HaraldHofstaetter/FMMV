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

	_FLOAT_ dx, dy, dz;
	_FLOAT_ Q;
	_FLOAT_ xi, yi, zi, qi;
	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	_FLOAT_ mxi, myi, mzi;
	#endif
	int i;
	
	dx = 0.0;
	dy = 0.0;
	dz = 0.0;
	Q = 0.0;
	
	for (i=0; i<FMMV->NParticles; i++) {
		xi = access_x(i);
		yi = access_y(i);
		zi = access_z(i);
		qi = access_q(i);
		
		dx = dx + qi*xi;
		dy = dy + qi*yi;
		dz = dz + qi*zi;

		Q = Q + qi*(xi*xi + yi*yi + zi*zi);

		#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
		   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		mxi = access_mx(i);
		myi = access_my(i);
		mzi = access_mz(i);

		dx = dx + mxi;
		dy = dy + myi;
		dz = dz + mzi;
		
		Q = Q + 2*(mxi*xi + myi*yi + mzi*zi);
		#endif
	}

	dx = -4.1887902047863909846*dx; /* -4*pi/3 */
	dy = -4.1887902047863909846*dy; /* -4*pi/3 */
	dz = -4.1887902047863909846*dz; /* -4*pi/3 */

	Q = 2.0943951023931954923*Q; /* 2*pi/3 */

	for (i=0; i<FMMV->NTargets; i++) {
		xi = access_tx(i); 
		yi = access_ty(i);
		zi = access_tz(i);
		
		access_pot(i) += dx*xi + dy*yi + dz*zi + Q;

		#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
		   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		access_gradx(i) += dx;
		access_grady(i) += dy;
		access_gradz(i) += dz;
		#endif
	}
}	
		
	
 
