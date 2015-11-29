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

void Rz_simd4(int p, V4_BASETYPE *x, V4_BASETYPE *y, int back)
{
	V4_INIT
	V4_TYPE ppmm;
	V4_TYPE mmpp;
	V4_TYPE mpmp;
	V4_TYPE pmpm;
	V4_TYPE mppm;
	V4_TYPE pmmp;
	V4_TYPE mmmm;
	V4_TYPE s2 = V4_SET1(0.70710678118654752440);

	V4_BASETYPE *yp = y;
	V4_BASETYPE *xp = x;
	int n,k;

	
	mmmm = V4_SET(-0.0, -0.0, -0.0, -0.0);
	ppmm = V4_SET(+0.0, +0.0, -0.0, -0.0);
	mmpp = V4_SET(-0.0, -0.0, +0.0, +0.0);
	if (!back) /* forward: p -> +, m -> - */ {
		mpmp = V4_SET(-0.0, +0.0, -0.0, +0.0);
		pmpm = V4_SET(+0.0, -0.0, +0.0, -0.0);
		mppm = V4_SET(-0.0, +0.0, +0.0, -0.0);
		pmmp = V4_SET(+0.0, -0.0, -0.0, +0.0);
	}
	else /* back: exchange 1<->2 and 3<->4 */ {
		mpmp = V4_SET(+0.0, -0.0, +0.0, -0.0);
		pmpm = V4_SET(-0.0, +0.0, -0.0, +0.0);
		mppm = V4_SET(+0.0, -0.0, -0.0, +0.0);
		pmmp = V4_SET(-0.0, +0.0, +0.0, -0.0);
	}
	for (n=0; n<=p; n++) {
		for (k=8; k<=n; k+=8) {
			V4_SETITEM(yp, 0, V4_GETITEM(xp, 0));
			V4_SETITEM(yp, 4, V4_GETITEM(xp, 4));
			V4_SETITEM(yp, 8, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETITEM(xp, 8), ppmm), V4_XOR(V4_GETITEM(xp, 12), mpmp))));
			V4_SETITEM(yp, 12, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETITEM(xp, 8), pmpm), V4_XOR(V4_GETITEM(xp, 12), ppmm))));
			V4_SETITEM(yp, 16, V4_XOR(V4_GETITEM(xp, 20), mppm));
			V4_SETITEM(yp, 20, V4_XOR(V4_GETITEM(xp, 16), pmmp));
			V4_SETITEM(yp, 24, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETITEM(xp, 24), mmpp), V4_XOR(V4_GETITEM(xp, 28), mpmp))));
			V4_SETITEM(yp, 28, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETITEM(xp, 24), pmpm), V4_XOR(V4_GETITEM(xp, 28), mmpp))));
			V4_SETITEM(yp, 32, V4_XOR(V4_GETITEM(xp, 32), mmmm));
			V4_SETITEM(yp, 36, V4_XOR(V4_GETITEM(xp, 36), mmmm));
			V4_SETITEM(yp, 40, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETITEM(xp, 40), mmpp), V4_XOR(V4_GETITEM(xp, 44), pmpm))));
			V4_SETITEM(yp, 44, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETITEM(xp, 40), mpmp), V4_XOR(V4_GETITEM(xp, 44), mmpp))));
			V4_SETITEM(yp, 48, V4_XOR(V4_GETITEM(xp, 52), pmmp));
			V4_SETITEM(yp, 52, V4_XOR(V4_GETITEM(xp, 48), mppm));
			V4_SETITEM(yp, 56, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETITEM(xp, 56), ppmm), V4_XOR(V4_GETITEM(xp, 60), pmpm))));
			V4_SETITEM(yp, 60, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETITEM(xp, 56), mpmp), V4_XOR(V4_GETITEM(xp, 60), ppmm))));
			yp += 64;
			xp += 64;
		}	
		switch (n+8-k) {/* Note: no break statements! */
		case 7:	
			V4_SETITEM(yp, 56, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETITEM(xp, 56), ppmm), V4_XOR(V4_GETITEM(xp, 60), pmpm))));
			V4_SETITEM(yp, 60, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETITEM(xp, 56), mpmp), V4_XOR(V4_GETITEM(xp, 60), ppmm))));
		case 6:	
			V4_SETITEM(yp, 48, V4_XOR(V4_GETITEM(xp, 52), pmmp));
			V4_SETITEM(yp, 52, V4_XOR(V4_GETITEM(xp, 48), mppm));
		case 5:	
			V4_SETITEM(yp, 40, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETITEM(xp, 40), mmpp), V4_XOR(V4_GETITEM(xp, 44), pmpm))));
			V4_SETITEM(yp, 44, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETITEM(xp, 40), mpmp), V4_XOR(V4_GETITEM(xp, 44), mmpp))));
		case 4:	
			V4_SETITEM(yp, 32, V4_XOR(V4_GETITEM(xp, 32), mmmm));
			V4_SETITEM(yp, 36, V4_XOR(V4_GETITEM(xp, 36), mmmm));
		case 3:	
			V4_SETITEM(yp, 24, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETITEM(xp, 24), mmpp), V4_XOR(V4_GETITEM(xp, 28), mpmp))));
			V4_SETITEM(yp, 28, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETITEM(xp, 24), pmpm), V4_XOR(V4_GETITEM(xp, 28), mmpp))));
		case 2:	
			V4_SETITEM(yp, 16, V4_XOR(V4_GETITEM(xp, 20), mppm));
			V4_SETITEM(yp, 20, V4_XOR(V4_GETITEM(xp, 16), pmmp));
		case 1:	
			V4_SETITEM(yp, 8, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETITEM(xp, 8), ppmm), V4_XOR(V4_GETITEM(xp, 12), mpmp))));
			V4_SETITEM(yp, 12, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETITEM(xp, 8), pmpm), V4_XOR(V4_GETITEM(xp, 12), ppmm))));
		case 0:	
			V4_SETITEM(yp, 0, V4_GETITEM(xp, 0));
			V4_SETITEM(yp, 4, V4_GETITEM(xp, 4));
			
			yp += 8*(n+9-k);
			xp += 8*(n+9-k);
		}
	}
}	

void Rz1_simd4(int p, V4_BASETYPE *x, V4_BASETYPE *y, int back)
{
	V4_INIT
	V4_TYPE ppmm;
	V4_TYPE mmpp;
	V4_TYPE mpmp;
	V4_TYPE pmpm;
	V4_TYPE mppm;
	V4_TYPE pmmp;
	V4_TYPE mmmm;
	V4_TYPE s2 = V4_SET1(0.70710678118654752440);

	V4_BASETYPE *yp = y;
	V4_BASETYPE *xp = x;
	int n,k;

	
	mmmm = V4_SET(-0.0, -0.0, -0.0, -0.0);
	ppmm = V4_SET(+0.0, +0.0, -0.0, -0.0);
	mmpp = V4_SET(-0.0, -0.0, +0.0, +0.0);
	if (!back) /* forward: p -> +, m -> - */ {
		mpmp = V4_SET(-0.0, +0.0, -0.0, +0.0);
		pmpm = V4_SET(+0.0, -0.0, +0.0, -0.0);
		mppm = V4_SET(-0.0, +0.0, +0.0, -0.0);
		pmmp = V4_SET(+0.0, -0.0, -0.0, +0.0);
	}
	else /* back: exchange 1<->2 and 3<->4 */ {
		mpmp = V4_SET(+0.0, -0.0, +0.0, -0.0);
		pmpm = V4_SET(-0.0, +0.0, -0.0, +0.0);
		mppm = V4_SET(+0.0, -0.0, -0.0, +0.0);
		pmmp = V4_SET(-0.0, +0.0, +0.0, -0.0);
	}

	for (n=0; n<=p; n++) {
		for (k=8; k<=n; k+=8) {
			V4_SETITEM(yp, 0, V4_GETBASEITEM(xp, 0));
			V4_SETITEM(yp, 4, V4_GETBASEITEM(xp, 1));
			V4_SETITEM(yp, 8, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETBASEITEM(xp, 2), ppmm), V4_XOR(V4_GETBASEITEM(xp, 3), mpmp))));
			V4_SETITEM(yp, 12, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETBASEITEM(xp, 2), pmpm), V4_XOR(V4_GETBASEITEM(xp, 3), ppmm))));
			V4_SETITEM(yp, 16, V4_XOR(V4_GETBASEITEM(xp, 5), mppm));
			V4_SETITEM(yp, 20, V4_XOR(V4_GETBASEITEM(xp, 4), pmmp));
			V4_SETITEM(yp, 24, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETBASEITEM(xp, 6), mmpp), V4_XOR(V4_GETBASEITEM(xp, 7), mpmp))));
			V4_SETITEM(yp, 28, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETBASEITEM(xp, 6), pmpm), V4_XOR(V4_GETBASEITEM(xp, 7), mmpp))));
			V4_SETITEM(yp, 32, V4_XOR(V4_GETBASEITEM(xp, 8), mmmm));
			V4_SETITEM(yp, 36, V4_XOR(V4_GETBASEITEM(xp, 9), mmmm));
			V4_SETITEM(yp, 40, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETBASEITEM(xp, 10), mmpp), V4_XOR(V4_GETBASEITEM(xp, 11), pmpm))));
			V4_SETITEM(yp, 44, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETBASEITEM(xp, 10), mpmp), V4_XOR(V4_GETBASEITEM(xp, 11), mmpp))));
			V4_SETITEM(yp, 48, V4_XOR(V4_GETBASEITEM(xp, 13), pmmp));
			V4_SETITEM(yp, 52, V4_XOR(V4_GETBASEITEM(xp, 12), mppm));
			V4_SETITEM(yp, 56, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETBASEITEM(xp, 14), ppmm), V4_XOR(V4_GETBASEITEM(xp, 15), pmpm))));
			V4_SETITEM(yp, 60, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETBASEITEM(xp, 14), mpmp), V4_XOR(V4_GETBASEITEM(xp, 15), ppmm))));
			yp += 64;
			xp += 16;
		}	
		switch (n+8-k) {/* Note: no break statements! */
		case 7:	
			V4_SETITEM(yp, 56, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETBASEITEM(xp, 14), ppmm), V4_XOR(V4_GETBASEITEM(xp, 15), pmpm))));
			V4_SETITEM(yp, 60, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETBASEITEM(xp, 14), mpmp), V4_XOR(V4_GETBASEITEM(xp, 15), ppmm))));
		case 6:	
			V4_SETITEM(yp, 48, V4_XOR(V4_GETBASEITEM(xp, 13), pmmp));
			V4_SETITEM(yp, 52, V4_XOR(V4_GETBASEITEM(xp, 12), mppm));
		case 5:	
			V4_SETITEM(yp, 40, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETBASEITEM(xp, 10), mmpp), V4_XOR(V4_GETBASEITEM(xp, 11), pmpm))));
			V4_SETITEM(yp, 44, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETBASEITEM(xp, 10), mpmp), V4_XOR(V4_GETBASEITEM(xp, 11), mmpp))));
		case 4:	
			V4_SETITEM(yp, 32, V4_XOR(V4_GETBASEITEM(xp, 8), mmmm));
			V4_SETITEM(yp, 36, V4_XOR(V4_GETBASEITEM(xp, 9), mmmm));
		case 3:	
			V4_SETITEM(yp, 24, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETBASEITEM(xp, 6), mmpp), V4_XOR(V4_GETBASEITEM(xp, 7), mpmp))));
			V4_SETITEM(yp, 28, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETBASEITEM(xp, 6), pmpm), V4_XOR(V4_GETBASEITEM(xp, 7), mmpp))));
		case 2:	
			V4_SETITEM(yp, 16, V4_XOR(V4_GETBASEITEM(xp, 5), mppm));
			V4_SETITEM(yp, 20, V4_XOR(V4_GETBASEITEM(xp, 4), pmmp));
		case 1:	
			V4_SETITEM(yp, 8, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETBASEITEM(xp, 2), ppmm), V4_XOR(V4_GETBASEITEM(xp, 3), mpmp))));
			V4_SETITEM(yp, 12, V4_MUL(s2, V4_ADD(V4_XOR(V4_GETBASEITEM(xp, 2), pmpm), V4_XOR(V4_GETBASEITEM(xp, 3), ppmm))));
		case 0:	
			V4_SETITEM(yp, 0, V4_GETBASEITEM(xp, 0));
			V4_SETITEM(yp, 4, V4_GETBASEITEM(xp, 1));
			
			yp += 8*(n+9-k);
			xp += 2*(n+9-k);
		}
	}
}	

void Rz_pi2_rrii_simd4(int p, V4_BASETYPE *x, V4_BASETYPE *y)
{
	V4_INIT
	V4_BASETYPE *yr = y;
	V4_BASETYPE *xr = x;
	V4_BASETYPE *yi;
	V4_BASETYPE *xi;
	int n,k;
	
	for (n=0; n<=p; n++) {
		V4_SETITEM(yr, 0, V4_GETITEM(xr, 0));
		yr += 4;
		xr += 4;
		yi = yr + 4*n;
		xi = xr + 4*n;
		for (k=5; k<=n; k+=4) {
			V4_SETITEM(yr, 0, V4_NEG(V4_GETITEM(xi, 0)));
			V4_SETITEM(yi, 0, V4_GETITEM(xr, 0));
			V4_SETITEM(yr, 4, V4_NEG(V4_GETITEM(xr, 4)));
			V4_SETITEM(yi, 4, V4_NEG(V4_GETITEM(xi, 4)));
			V4_SETITEM(yr, 8, V4_GETITEM(xi, 8));
			V4_SETITEM(yi, 8, V4_NEG(V4_GETITEM(xr, 8)));
			V4_SETITEM(yr, 12, V4_GETITEM(xr, 12));
			V4_SETITEM(yi, 12, V4_GETITEM(xi, 12));
			yr += 16;
			xr += 16;
			yi += 16;
			xi += 16;
		}
		switch (n+4-k) {/* Note: no break statements! */
		case 3:	
			V4_SETITEM(yr, 12, V4_GETITEM(xr, 12));
			V4_SETITEM(yi, 12, V4_GETITEM(xi, 12));
		case 2:	
			V4_SETITEM(yr, 8, V4_GETITEM(xi, 8));
			V4_SETITEM(yi, 8, V4_NEG(V4_GETITEM(xr, 8)));
		case 1:	
			V4_SETITEM(yr, 4, V4_NEG(V4_GETITEM(xr, 4)));
			V4_SETITEM(yi, 4, V4_NEG(V4_GETITEM(xi, 4)));
		case 0:
			V4_SETITEM(yr, 0, V4_NEG(V4_GETITEM(xi, 0)));
			V4_SETITEM(yi, 0, V4_GETITEM(xr, 0));
			
			yi += 4*(n+5-k);
			xi += 4*(n+5-k);
		}	
		yr = yi;
		xr = xi;
	}	
}	


void Rz_minus_pi2_rrii_simd4(int p, V4_BASETYPE *x, V4_BASETYPE *y)
{
	V4_INIT
	V4_BASETYPE *yr = y;
	V4_BASETYPE *xr = x;
	V4_BASETYPE *yi;
	V4_BASETYPE *xi;
	int n,k;
	
	for (n=0; n<=p; n++) {
		V4_SETITEM(yr, 0, V4_GETITEM(xr, 0));
		yr += 4;
		xr += 4;
		yi = yr + 4*n;
		xi = xr + 4*n;
		for (k=5; k<=n; k+=4) {
			V4_SETITEM(yr, 0, V4_GETITEM(xi, 0));
			V4_SETITEM(yi, 0, V4_NEG(V4_GETITEM(xr, 0)));
			V4_SETITEM(yr, 4, V4_NEG(V4_GETITEM(xr, 4)));
			V4_SETITEM(yi, 4, V4_NEG(V4_GETITEM(xi, 4)));
			V4_SETITEM(yr, 8, V4_NEG(V4_GETITEM(xi, 8)));
			V4_SETITEM(yi, 8, V4_GETITEM(xr, 8));
			V4_SETITEM(yr, 12, V4_GETITEM(xr, 12));
			V4_SETITEM(yi, 12, V4_GETITEM(xi, 12));
			yr += 16;
			xr += 16;
			yi += 16;
			xi += 16;
		}
		switch (n+4-k) {/* Note: no break statements! */
		case 3:	
			V4_SETITEM(yr, 12, V4_GETITEM(xr, 12));
			V4_SETITEM(yi, 12, V4_GETITEM(xi, 12));
		case 2:	
			V4_SETITEM(yr, 8, V4_NEG(V4_GETITEM(xi, 8)));
			V4_SETITEM(yi, 8, V4_GETITEM(xr, 8));
		case 1:	
			V4_SETITEM(yr, 4, V4_NEG(V4_GETITEM(xr, 4)));
			V4_SETITEM(yi, 4, V4_NEG(V4_GETITEM(xi, 4)));
		case 0:
			V4_SETITEM(yr, 0, V4_GETITEM(xi, 0));
			V4_SETITEM(yi, 0, V4_NEG(V4_GETITEM(xr, 0)));

			yi += 4*(n+5-k);
			xi += 4*(n+5-k);
		}	
		yr = yi;
		xr = xi;
	}	
}	








