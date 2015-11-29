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

#include "FFT_M2X_X2L_tables_simd2.h"


#include "simd.h"
#if (FMM_PRECISION==0)
   #include"simd2s.h"
#else
   #include"simd2d.h"
#endif

void F_M2X_simd2(FmmvHandle *FMMV, _FLOAT_ *x, _FLOAT_ *y)
{
	int *M = FMMV->M;
	int s_eps = FMMV->s_eps;
	_FLOAT_ *xp = x;
	_FLOAT_ *yp = y;
	int i,j, m; 
	V2_TYPE s;
	SIMD_ALIGN _FLOAT_ h[128*2*2];
	SIMD_ALIGN _FLOAT_ hh[128*2*2];

	for (i=0; i<s_eps; i++) {
		m = M[i];
		
		VEC_COPY((m+2)<<1, xp, h);
                xp += (m+2)<<1;
                for (j=1; j<FMMV->FFT_rep[i]; j++) {
                       VEC_ADD((m+2)<<1, h, xp, h);
                       xp += (m+2)<<1;
                }
		
		FFT_M2X_prepare_simd2[m>>2](h, hh, &s);
		FFT_M2X_simd2[m>>2](yp, hh);
		FFT_M2X_finish_simd2[m>>2](yp, s);
		
		yp += m<<1;
	}
 
}

void F_X2L_simd2(FmmvHandle *FMMV, _FLOAT_ *x, _FLOAT_ *y)
{
	int *M = FMMV->M;
	int s_eps = FMMV->s_eps;
	_FLOAT_ *xp = x;
	_FLOAT_ *yp = y;
	int i, m;
	V2_TYPE s;
	SIMD_ALIGN _FLOAT_ h[128*2*2];

	for (i=0; i<s_eps; i++) {
		m = M[i];
		
		FFT_X2L_prepare_simd2[m>>2](xp, &s);
		FFT_X2L_simd2[m>>2](h, xp);
		FFT_X2L_finish_simd2[m>>2](h, yp, s);
		
		yp += (m+2)<<1;
		xp += m<<1;
	}
}	
		
