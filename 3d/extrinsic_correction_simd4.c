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

#include"simd.h"
#if (FMM_PRECISION==0)
   #include"simd4s.h"
#else
   #include"simd4d.h"
#endif

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

	V4_TYPE dx, dy, dz;
	V4_TYPE Q;
	V4_TYPE xi, yi, zi, qi;
	_FLOAT_ dxf, dyf, dzf, Qf;
        #if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
           ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	V4_TYPE mxi, myi, mzi;
	V4_TYPE two = V4_SET1(2.0);
	#endif
	int i;
	
	dx = V4_ZERO;
	dy = V4_ZERO;
	dz = V4_ZERO;
	Q = V4_ZERO;
	
	for (i = 0; i < FMMV->NParticles>>2; i++) {
		xi = V4_LOAD(access_x(i));
		yi = V4_LOAD(access_y(i));
		zi = V4_LOAD(access_z(i));
		qi = V4_LOAD(access_q(i));

		dx = V4_ADD(dx, V4_MUL(qi, xi));
		dy = V4_ADD(dy, V4_MUL(qi, yi));
		dz = V4_ADD(dz, V4_MUL(qi, zi));

		Q = V4_ADD(Q, V4_MUL(qi, 
		       V4_ADD(V4_ADD(V4_MUL(xi, xi), V4_MUL(yi, yi)), V4_MUL(zi, zi))));

        	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	           ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		mxi = V4_LOAD(access_mx(i));
		myi = V4_LOAD(access_my(i));
		mzi = V4_LOAD(access_mz(i));

		dx = V4_ADD(dx, mxi);
		dy = V4_ADD(dy, myi);
		dz = V4_ADD(dz, mzi);
		
		Q = V4_ADD(Q, V4_MUL(two, V4_ADD(V4_ADD(V4_MUL(mxi,xi), V4_MUL(myi,yi)), V4_MUL(mzi,zi))));
		#endif
	}
	
	dxf = *((_FLOAT_*) &dx) + *(((_FLOAT_ *) &dx)+1) + *(((_FLOAT_ *) &dx)+2) + *(((_FLOAT_ *) &dx)+3);
	dyf = *((_FLOAT_*) &dy) + *(((_FLOAT_ *) &dy)+1) + *(((_FLOAT_ *) &dy)+2) + *(((_FLOAT_ *) &dy)+3);
	dzf = *((_FLOAT_*) &dz) + *(((_FLOAT_ *) &dz)+1) + *(((_FLOAT_ *) &dz)+2) + *(((_FLOAT_ *) &dz)+3);
	Qf = *((_FLOAT_*) &Q) + *(((_FLOAT_ *) &Q)+1) + *(((_FLOAT_ *) &Q)+2) + *(((_FLOAT_ *) &Q)+3);

	if (FMMV->NParticles&3) {
		int i = FMMV->NParticles>>2;
		int j;
		_FLOAT_ xif, yif, zif, qif;
        	#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
	           ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		_FLOAT_ mxif, myif, mzif;
		#endif

		for (j=0; j<(FMMV->NParticles&3); j++) {
			xif = access_x(i)[j];
			yif = access_y(i)[j];
			zif = access_z(i)[j];
			qif = access_q(i)[j];
		
			dxf = dxf + qif*xif;
			dyf = dyf + qif*yif;
			dzf = dzf + qif*zif;

			Qf = Qf + qif*(xif*xif + yif*yif + zif*zif);

        		#if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
		           ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			mxif = access_mx(i)[j];
			myif = access_my(i)[j];
			mzif = access_mz(i)[j];
			dxf = dxf + mxif;
			dyf = dyf + myif;
			dzf = dzf + mzif;
		
			Qf = Qf + 2*(mxif*xif + myif*yif + mzif*zif);
			#endif
		}	
	}

	dxf = -4.1887902047863909846*dxf; /* -4*pi/3 */
	dyf = -4.1887902047863909846*dyf; /* -4*pi/3 */
	dzf = -4.1887902047863909846*dzf; /* -4*pi/3 */
	Qf = 2.0943951023931954923*Qf; /* 2*pi/3 */

	dx = V4_LOAD1(&dxf);
	dy = V4_LOAD1(&dyf);
	dz = V4_LOAD1(&dzf);
	Q = V4_LOAD1(&Qf);
	
	for (i = 0; i < FMMV->NTargets>>2; i++) {
		xi = V4_LOAD(access_tx(i));
		yi = V4_LOAD(access_ty(i));
		zi = V4_LOAD(access_tz(i));

		V4_STORE(access_pot(i), V4_ADD(V4_LOAD(access_pot(i)), V4_ADD(V4_ADD(V4_ADD(V4_MUL(dx, xi), V4_MUL(dy, yi)), V4_MUL(dz, zi)),Q)));

		
        	#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
  		   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		V4_STORE(access_gradx(i), V4_ADD(V4_LOAD(access_gradx(i)), dx));
		V4_STORE(access_grady(i), V4_ADD(V4_LOAD(access_grady(i)), dy));
		V4_STORE(access_gradz(i), V4_ADD(V4_LOAD(access_gradz(i)), dz));
		#endif
	}	
	
	if (FMMV->NTargets&3) {
		int i = FMMV->NTargets>>2;
		int j;
		_FLOAT_ xif, yif, zif;

		for (j=0; j<(FMMV->NTargets&3); j++) {
			xif = access_tx(i)[j];
			yif = access_ty(i)[j];
			zif = access_tz(i)[j];

			access_pot(i)[j] += dxf*xif + dyf*yif + dzf*zif + Qf;
		
        		#if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
			   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
			access_gradx(i)[j] += dxf;
			access_grady(i)[j] += dyf;
			access_gradz(i)[j] += dzf;
			#endif
		}	
	}
}
