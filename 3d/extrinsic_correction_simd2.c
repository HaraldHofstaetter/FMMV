/*
 * FMMV - the Fastest Multipole Method of Vienna
 * Copyright (c) 2006 Harald Hofstaetter
 * Institute for Analysis and Scientific Computing
 * Vienna University of Technology
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
   #include"simd2s.h"
#else
   #include"simd2d.h"
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

	V2_TYPE dx, dy, dz;
	V2_TYPE Q;
	V2_TYPE xi, yi, zi, qi;
	_FLOAT_ dxf, dyf, dzf, Qf;
        #if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
           ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
	V2_TYPE mxi, myi, mzi;
	V2_TYPE two = V2_SET1(2.0);
	#endif
	int i;
	
	dx = V2_ZERO;
	dy = V2_ZERO;
	dz = V2_ZERO;
	Q = V2_ZERO;
	
	for (i = 0; i < FMMV->NParticles>>1; i++) {
		xi = V2_LOAD(access_x(i));
		yi = V2_LOAD(access_y(i));
		zi = V2_LOAD(access_z(i));
		qi = V2_LOAD(access_q(i));

		dx = V2_ADD(dx, V2_MUL(qi, xi));
		dy = V2_ADD(dy, V2_MUL(qi, yi));
		dz = V2_ADD(dz, V2_MUL(qi, zi));

		Q = V2_ADD(Q, V2_MUL(qi, 
		       V2_ADD(V2_ADD(V2_MUL(xi, xi), V2_MUL(yi, yi)), V2_MUL(zi, zi))));

                #if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
                   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		mxi = V2_LOAD(access_mx(i));
		myi = V2_LOAD(access_my(i));
		mzi = V2_LOAD(access_mz(i));

		dx = V2_ADD(dx, mxi);
		dy = V2_ADD(dy, myi);
		dz = V2_ADD(dz, mzi);
		
		Q = V2_ADD(Q, V2_MUL(two, V2_ADD(V2_ADD(V2_MUL(mxi, xi), V2_MUL(myi, yi)), V2_MUL(mzi, zi))));
		#endif
	}
	
	dxf = *((_FLOAT_*) &dx) + *(((_FLOAT_ *) &dx)+1);
	dyf = *((_FLOAT_*) &dy) + *(((_FLOAT_ *) &dy)+1);
	dzf = *((_FLOAT_*) &dz) + *(((_FLOAT_ *) &dz)+1);

	Qf = *((_FLOAT_*) &Q) + *(((_FLOAT_ *) &Q)+1);
	
	if (FMMV->NParticles&1) {
		int i = FMMV->NParticles>>1;
		_FLOAT_ xif = access_x(i)[0];
		_FLOAT_ yif = access_y(i)[0];
		_FLOAT_ zif = access_z(i)[0];
		_FLOAT_ qif = access_q(i)[0];
                #if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
                   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		_FLOAT_ mxif = access_mx(i)[0];
		_FLOAT_ myif = access_my(i)[0];
		_FLOAT_ mzif = access_mz(i)[0];
		#endif
		
		dxf = dxf + qif*xif;
		dyf = dyf + qif*yif;
		dzf = dzf + qif*zif;

		Qf = Qf + qif*(xif*xif + yif*yif + zif*zif);

                #if ((FMM_KIND==FMM_DIPOLE)||(FMM_KIND==FMM_DIPOLE_GRAD) \
                   ||(FMM_KIND==FMM_ST_DIPOLE)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		dxf = dxf + mxif;
		dyf = dyf + myif;
		dzf = dzf + mzif;
		
		Qf = Qf + 2.0*(mxif*xif + myif*yif + mzif*zif);
		#endif
	}
	
	dxf = -4.1887902047863909846*dxf; /* -4*pi/3 */
	dyf = -4.1887902047863909846*dyf; /* -4*pi/3 */
	dzf = -4.1887902047863909846*dzf; /* -4*pi/3 */

	Qf = 2.0943951023931954923*Qf; /* 2*pi/3 */

	dx = V2_LOAD1(&dxf);
	dy = V2_LOAD1(&dyf);
	dz = V2_LOAD1(&dzf);
	Q = V2_LOAD1(&Qf);
	
	for (i = 0; i < FMMV->NTargets>>1; i++) {
		xi = V2_LOAD(access_tx(i));
		yi = V2_LOAD(access_ty(i));
		zi = V2_LOAD(access_tz(i));

		V2_STORE(access_pot(i), V2_ADD(V2_LOAD(access_pot(i)), V2_ADD(V2_ADD(V2_ADD(V2_MUL(dx, xi), V2_MUL(dy, yi)), V2_MUL(dz, zi)),Q)));
		
                #if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
                   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		V2_STORE(access_gradx(i), V2_ADD(V2_LOAD(access_gradx(i)), dx));
		V2_STORE(access_grady(i), V2_ADD(V2_LOAD(access_grady(i)), dy));
		V2_STORE(access_gradz(i), V2_ADD(V2_LOAD(access_gradz(i)), dz));
		#endif
	}	
	
	if (FMMV->NTargets&1) {
		int i = FMMV->NTargets>>1;
		_FLOAT_ xif = access_tx(i)[0];
		_FLOAT_ yif = access_ty(i)[0];
		_FLOAT_ zif = access_tz(i)[0];

		access_pot(i)[0] += dxf*xif + dyf*yif + dzf*zif + Qf;
		
                #if ((FMM_KIND==FMM_GRAD)||(FMM_KIND==FMM_DIPOLE_GRAD) \
                   ||(FMM_KIND==FMM_ST_GRAD)||(FMM_KIND==FMM_ST_DIPOLE_GRAD))
		access_gradx(i)[0] += dxf;
		access_grady(i)[0] += dyf;
		access_gradz(i)[0] += dzf;
		#endif
	}
}
