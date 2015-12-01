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


#if (defined(C_CHARGE))
	#define FMM_KIND FMM_C_CHARGE	
	#define extrinsic_correction extrinsic_correction_c_charge
#elif (defined(C_DIPOLE))
	#define FMM_KIND FMM_C_DIPOLE	
	#define extrinsic_correction extrinsic_correction_c_dipole
#elif (defined(C_CHARGE_DIPOLE))
	#define FMM_KIND FMM_C_CHARGE_DIPOLE	
	#define extrinsic_correction extrinsic_correction_c_charge_dipole
#elif (defined(C_CHARGE_GRAD))
	#define FMM_KIND FMM_C_CHARGE_GRAD	
	#define extrinsic_correction extrinsic_correction_c_charge_grad
#elif (defined(C_DIPOLE_GRAD))
	#define FMM_KIND FMM_C_DIPOLE_GRAD	
	#define extrinsic_correction extrinsic_correction_c_dipole_grad
#elif (defined(C_CHARGE_DIPOLE_GRAD))
	#define FMM_KIND FMM_C_CHARGE_DIPOLE_GRAD	
	#define extrinsic_correction extrinsic_correction_c_charge_dipole_grad
#elif (defined(ST_C_CHARGE))
	#define FMM_KIND FMM_ST_C_CHARGE	
	#define extrinsic_correction extrinsic_correction_ST_c_charge
#elif (defined(ST_C_DIPOLE))
	#define FMM_KIND FMM_ST_C_DIPOLE	
	#define extrinsic_correction extrinsic_correction_ST_c_dipole
#elif (defined(ST_C_CHARGE_DIPOLE))
	#define FMM_KIND FMM_ST_C_CHARGE_DIPOLE	
	#define extrinsic_correction extrinsic_correction_ST_c_charge_dipole
#elif (defined(ST_C_CHARGE_GRAD))
	#define FMM_KIND FMM_ST_C_CHARGE_GRAD	
	#define extrinsic_correction extrinsic_correction_ST_c_charge_grad
#elif (defined(ST_C_DIPOLE_GRAD))
	#define FMM_KIND FMM_ST_C_DIPOLE_GRAD	
	#define extrinsic_correction extrinsic_correction_ST_c_dipole_grad
#elif (defined(ST_C_CHARGE_DIPOLE_GRAD))
	#define FMM_KIND FMM_ST_C_CHARGE_DIPOLE_GRAD	
	#define extrinsic_correction extrinsic_correction_ST_c_charge_dipole_grad
#endif


#include"fmmv_access.h"

/* THIS IS VER DUBIOUS!!!  HAS TO BE CHECKED !!! */

void extrinsic_correction(FmmvHandle *FMMV)
{
	DEFINE_IDA_LOCAL_ALIASES(FMMV)

	_FLOAT_ dx_re, dx_im, dy_re, dy_im;
	_FLOAT_ Q_re, Q_im, r;
	_FLOAT_ xi, yi;
        #if(defined(access_qx))
        _FLOAT_ qi_re, qi_im;
	#endif
        #if(defined(access_mx))
	_FLOAT_ mxi, myi;
	#endif
	int i;
        _FLOAT_ PI = 3.14159265358979323846;
	
	dx_re = 0.0; 
	dx_im = 0.0; 
	dy_re = 0.0;
	dy_im = 0.0;
	Q_re = 0.0;
	Q_im = 0.0;
	
	for (i=0; i<FMMV->NParticles; i++) {
                xi = access_x(i); 
		yi = access_y(i);
                #if(defined(access_qx))
		qi_re = access_qx(i);
		qi_re = access_qy(i);

		dx_re += qi_re*xi; 
		dx_im += qi_im*xi; 
		dy_re += qi_re*yi;
		dy_im += qi_im*yi;

                r = xi*xi + yi*yi;
		Q_re +=  r*qi_re;
		Q_im +=  r*qi_im;
                #endif
                #if(defined(access_mx))
		mxi = access_mx(i); 
		myi = access_my(i);

		dx_re +=  mxi; 
		dy_re -=  myi; 
		
		Q_re +=  2*(  mxi*xi -  myi*yi); 
		#endif
	}

	Q_re *= -PI;
	Q_im *= -PI;
        dx_re *= 2*PI; 
        dx_im *= 2*PI; 
        dy_re *= 2*PI;
        dy_im *= 2*PI;

	for (i=0; i<FMMV->NTargets; i++) {
		xi = access_tx(i); 
		yi = access_ty(i);
		
		access_potx(i) +=  dx_re*xi + dy_re*yi + Q_re;   /* real part */
		access_poty(i) +=  dy_im*xi + dx_im*yi + Q_im;   /* imaginary part */

                #if(defined(access_gradx))
		access_gradx(i) += dx_re; 
		access_grady(i) += dy_re;
		#endif
	}
}	
		
	
 
