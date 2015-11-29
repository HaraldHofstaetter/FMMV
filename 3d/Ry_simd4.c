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

#include "_fmmv.h"
#include "simd.h"
#if (FMM_PRECISION==0)
   #include"simd4s.h"
#else
   #include"simd4d.h"
#endif


void Ry_simd4(int p, _FLOAT_ **blocks, _FLOAT_ *x, _FLOAT_ *y)
{
	int k,j, j1;

	gemv_simd4(1, 1, blocks[0], 1, x, y);
	k = 4*1;
	for (j=1; j<=p; j++) {
		j1 = j+1;
		gemv_simd4(j1, j1, blocks[2*j], j1, x+k, y+k);
		k += 4*j1;
		gemv_simd4(j, j, blocks[2*j+1], j, x+k, y+k);
		k += 4*j;
	}
}	


void Ry_pi_simd4(int p, /*_FLOAT_ *x,*/ _FLOAT_ *y)
{
	int len = 4*(p+1)*(p+1);
	int dk, k0, k1, i;
	
	k0 = 4;
	k1 = 12;
	dk = 8;
	while (k1<=len) {
		for (i=k0; i<k1; i+=4) {
			V4_SETITEM(y, i, V4_NEG(V4_GETITEM(y, i)));
		}
		dk += 8;	
		k0 = k1 + dk;
		dk += 8;
		k1 = k0 + dk;
	}	
	for (i=k0; i<len; i+=4) {
		V4_SETITEM(y, i, V4_NEG(V4_GETITEM(y, i)));
	}
}	


	
