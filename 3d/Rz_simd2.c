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
   #include"simd2s.h"
#else
   #include"simd2d.h"
#endif

static V2_BASETYPE PPMM[4] = {+0.0, +0.0, -0.0, -0.0};
static V2_BASETYPE MMPP[4] = {-0.0, -0.0, +0.0, +0.0};
static V2_BASETYPE MPMP[4] = {-0.0, +0.0, -0.0, +0.0};
static V2_BASETYPE PMPM[4] = {+0.0, -0.0, +0.0, -0.0};
static V2_BASETYPE MPPM[4] = {-0.0, +0.0, +0.0, -0.0};
static V2_BASETYPE PMMP[4] = {+0.0, -0.0, -0.0, +0.0};

void Rz_simd2(int p, V2_BASETYPE *x, V2_BASETYPE *y, int k1, int k2, int back)
{
	V2_INIT
	V2_TYPE ppmm;
	V2_TYPE mmpp;
	V2_TYPE mpmp;
	V2_TYPE pmpm;
	V2_TYPE mppm;
	V2_TYPE pmmp;
	V2_TYPE mmmm;
	V2_TYPE s2 = V2_SET1(0.70710678118654752440);

	V2_BASETYPE *yp = y;
	V2_BASETYPE *xp = x;
	int n,k;

	mmmm = V2_SET(-0.0, -0.0);
	ppmm = V2_SET(PPMM[k1], PPMM[k2]);
	mmpp = V2_SET(MMPP[k1], MMPP[k2]);
	if (!back) {
		mpmp = V2_SET(MPMP[k1], MPMP[k2]);
		pmpm = V2_SET(PMPM[k1], PMPM[k2]);
		mppm = V2_SET(MPPM[k1], MPPM[k2]);
		pmmp = V2_SET(PMMP[k1], PMMP[k2]);
	}
	else /* back: exchange 1<->2 and 3<->4 */ {
		pmpm = V2_SET(MPMP[k1], MPMP[k2]);
		mpmp = V2_SET(PMPM[k1], PMPM[k2]);
		pmmp = V2_SET(MPPM[k1], MPPM[k2]);
		mppm = V2_SET(PMMP[k1], PMMP[k2]);
	}

	for (n=0; n<=p; n++) {
		for (k=8; k<=n; k+=8) {
			V2_SETITEM(yp, 0, V2_GETITEM(xp, 0));
			V2_SETITEM(yp, 2, V2_GETITEM(xp, 2));
			V2_SETITEM(yp, 4, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETITEM(xp, 4), ppmm), V2_XOR(V2_GETITEM(xp, 6), mpmp))));
			V2_SETITEM(yp, 6, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETITEM(xp, 4), pmpm), V2_XOR(V2_GETITEM(xp, 6), ppmm))));
			V2_SETITEM(yp, 8, V2_XOR(V2_GETITEM(xp, 10), mppm));
			V2_SETITEM(yp, 10, V2_XOR(V2_GETITEM(xp, 8), pmmp));
			V2_SETITEM(yp, 12, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETITEM(xp, 12), mmpp), V2_XOR(V2_GETITEM(xp, 14), mpmp))));
			V2_SETITEM(yp, 14, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETITEM(xp, 12), pmpm), V2_XOR(V2_GETITEM(xp, 14), mmpp))));
			V2_SETITEM(yp, 16, V2_XOR(V2_GETITEM(xp, 16), mmmm));
			V2_SETITEM(yp, 18, V2_XOR(V2_GETITEM(xp, 18), mmmm));
			V2_SETITEM(yp, 20, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETITEM(xp, 20), mmpp), V2_XOR(V2_GETITEM(xp, 22), pmpm))));
			V2_SETITEM(yp, 22, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETITEM(xp, 20), mpmp), V2_XOR(V2_GETITEM(xp, 22), mmpp))));
			V2_SETITEM(yp, 24, V2_XOR(V2_GETITEM(xp, 26), pmmp));
			V2_SETITEM(yp, 26, V2_XOR(V2_GETITEM(xp, 24), mppm));
			V2_SETITEM(yp, 28, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETITEM(xp, 28), ppmm), V2_XOR(V2_GETITEM(xp, 30), pmpm))));
			V2_SETITEM(yp, 30, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETITEM(xp, 28), mpmp), V2_XOR(V2_GETITEM(xp, 30), ppmm))));
			yp += 32;
			xp += 32;
		}	
		switch (n+8-k) {/* Note: no break statements! */
		case 7:	
			V2_SETITEM(yp, 28, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETITEM(xp, 28), ppmm), V2_XOR(V2_GETITEM(xp, 30), pmpm))));
			V2_SETITEM(yp, 30, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETITEM(xp, 28), mpmp), V2_XOR(V2_GETITEM(xp, 30), ppmm))));
		case 6:	
			V2_SETITEM(yp, 24, V2_XOR(V2_GETITEM(xp, 26), pmmp));
			V2_SETITEM(yp, 26, V2_XOR(V2_GETITEM(xp, 24), mppm));
		case 5:	
			V2_SETITEM(yp, 20, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETITEM(xp, 20), mmpp), V2_XOR(V2_GETITEM(xp, 22), pmpm))));
			V2_SETITEM(yp, 22, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETITEM(xp, 20), mpmp), V2_XOR(V2_GETITEM(xp, 22), mmpp))));
		case 4:	
			V2_SETITEM(yp, 16, V2_XOR(V2_GETITEM(xp, 16), mmmm));
			V2_SETITEM(yp, 18, V2_XOR(V2_GETITEM(xp, 18), mmmm));
		case 3:	
			V2_SETITEM(yp, 12, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETITEM(xp, 12), mmpp), V2_XOR(V2_GETITEM(xp, 14), mpmp))));
			V2_SETITEM(yp, 14, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETITEM(xp, 12), pmpm), V2_XOR(V2_GETITEM(xp, 14), mmpp))));
		case 2:	
			V2_SETITEM(yp, 8, V2_XOR(V2_GETITEM(xp, 10), mppm));
			V2_SETITEM(yp, 10, V2_XOR(V2_GETITEM(xp, 8), pmmp));
		case 1:	
			V2_SETITEM(yp, 4, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETITEM(xp, 4), ppmm), V2_XOR(V2_GETITEM(xp, 6), mpmp))));
			V2_SETITEM(yp, 6, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETITEM(xp, 4), pmpm), V2_XOR(V2_GETITEM(xp, 6), ppmm))));
		case 0:	
			V2_SETITEM(yp, 0, V2_GETITEM(xp, 0));
			V2_SETITEM(yp, 2, V2_GETITEM(xp, 2));
			
			yp += 4*(n+9-k);
			xp += 4*(n+9-k);
		}
	}
}	

	
void Rz1_simd2(int p, V2_BASETYPE *x, V2_BASETYPE *y, int k1, int k2, int back)
{
	V2_INIT
	V2_TYPE ppmm;
	V2_TYPE mmpp;
	V2_TYPE mpmp;
	V2_TYPE pmpm;
	V2_TYPE mppm;
	V2_TYPE pmmp;
	V2_TYPE mmmm;
	V2_TYPE s2 = V2_SET1(0.70710678118654752440);
	
	V2_BASETYPE *yp = y;
	V2_BASETYPE *xp = x;
	int n,k;
	
	mmmm = V2_SET(-0.0, -0.0);
	ppmm = V2_SET(PPMM[k1], PPMM[k2]);
	mmpp = V2_SET(MMPP[k1], MMPP[k2]);
	if (!back) {
		mpmp = V2_SET(MPMP[k1], MPMP[k2]);
		pmpm = V2_SET(PMPM[k1], PMPM[k2]);
		mppm = V2_SET(MPPM[k1], MPPM[k2]);
		pmmp = V2_SET(PMMP[k1], PMMP[k2]);
	}
	else /* back: exchange 1<->2 and 3<->4 */ {
		pmpm = V2_SET(MPMP[k1], MPMP[k2]);
		mpmp = V2_SET(PMPM[k1], PMPM[k2]);
		pmmp = V2_SET(MPPM[k1], MPPM[k2]);
		mppm = V2_SET(PMMP[k1], PMMP[k2]);
	}

	for (n=0; n<=p; n++) {
		for (k=8; k<=n; k+=8) {
			V2_SETITEM(yp, 0, V2_GETBASEITEM(xp, 0));
			V2_SETITEM(yp, 2, V2_GETBASEITEM(xp, 1));
			V2_SETITEM(yp, 4, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETBASEITEM(xp, 2), ppmm), V2_XOR(V2_GETBASEITEM(xp, 3), mpmp))));
			V2_SETITEM(yp, 6, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETBASEITEM(xp, 2), pmpm), V2_XOR(V2_GETBASEITEM(xp, 3), ppmm))));
			V2_SETITEM(yp, 8, V2_XOR(V2_GETBASEITEM(xp, 5), mppm));
			V2_SETITEM(yp, 10, V2_XOR(V2_GETBASEITEM(xp, 4), pmmp));
			V2_SETITEM(yp, 12, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETBASEITEM(xp, 6), mmpp), V2_XOR(V2_GETBASEITEM(xp, 7), mpmp))));
			V2_SETITEM(yp, 14, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETBASEITEM(xp, 6), pmpm), V2_XOR(V2_GETBASEITEM(xp, 7), mmpp))));
			V2_SETITEM(yp, 16, V2_XOR(V2_GETBASEITEM(xp, 8), mmmm));
			V2_SETITEM(yp, 18, V2_XOR(V2_GETBASEITEM(xp, 9), mmmm));
			V2_SETITEM(yp, 20, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETBASEITEM(xp, 10), mmpp), V2_XOR(V2_GETBASEITEM(xp, 11), pmpm))));
			V2_SETITEM(yp, 22, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETBASEITEM(xp, 10), mpmp), V2_XOR(V2_GETBASEITEM(xp, 11), mmpp))));
			V2_SETITEM(yp, 24, V2_XOR(V2_GETBASEITEM(xp, 13), pmmp));
			V2_SETITEM(yp, 26, V2_XOR(V2_GETBASEITEM(xp, 12), mppm));
			V2_SETITEM(yp, 28, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETBASEITEM(xp, 14), ppmm), V2_XOR(V2_GETBASEITEM(xp, 15), pmpm))));
			V2_SETITEM(yp, 30, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETBASEITEM(xp, 14), mpmp), V2_XOR(V2_GETBASEITEM(xp, 15), ppmm))));
			yp += 32;
			xp += 16;
		}	
		switch (n+8-k) {/* Note: no break statements! */
		case 7:	
			V2_SETITEM(yp, 28, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETBASEITEM(xp, 14), ppmm), V2_XOR(V2_GETBASEITEM(xp, 15), pmpm))));
			V2_SETITEM(yp, 30, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETBASEITEM(xp, 14), mpmp), V2_XOR(V2_GETBASEITEM(xp, 15), ppmm))));
		case 6:	
			V2_SETITEM(yp, 24, V2_XOR(V2_GETBASEITEM(xp, 13), pmmp));
			V2_SETITEM(yp, 26, V2_XOR(V2_GETBASEITEM(xp, 12), mppm));
		case 5:	
			V2_SETITEM(yp, 20, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETBASEITEM(xp, 10), mmpp), V2_XOR(V2_GETBASEITEM(xp, 11), pmpm))));
			V2_SETITEM(yp, 22, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETBASEITEM(xp, 10), mpmp), V2_XOR(V2_GETBASEITEM(xp, 11), mmpp))));
		case 4:	
			V2_SETITEM(yp, 16, V2_XOR(V2_GETBASEITEM(xp, 8), mmmm));
			V2_SETITEM(yp, 18, V2_XOR(V2_GETBASEITEM(xp, 9), mmmm));
		case 3:	
			V2_SETITEM(yp, 12, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETBASEITEM(xp, 6), mmpp), V2_XOR(V2_GETBASEITEM(xp, 7), mpmp))));
			V2_SETITEM(yp, 14, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETBASEITEM(xp, 6), pmpm), V2_XOR(V2_GETBASEITEM(xp, 7), mmpp))));
		case 2:	
			V2_SETITEM(yp, 8, V2_XOR(V2_GETBASEITEM(xp, 5), mppm));
			V2_SETITEM(yp, 10, V2_XOR(V2_GETBASEITEM(xp, 4), pmmp));
		case 1:	
			V2_SETITEM(yp, 4, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETBASEITEM(xp, 2), ppmm), V2_XOR(V2_GETBASEITEM(xp, 3), mpmp))));
			V2_SETITEM(yp, 6, V2_MUL(s2, V2_ADD(V2_XOR(V2_GETBASEITEM(xp, 2), pmpm), V2_XOR(V2_GETBASEITEM(xp, 3), ppmm))));
		case 0:	
			V2_SETITEM(yp, 0, V2_GETBASEITEM(xp, 0));
			V2_SETITEM(yp, 2, V2_GETBASEITEM(xp, 1));
			
			yp += 4*(n+9-k);
			xp += 2*(n+9-k);
		}
	}
}	


	
void Rz_pi2_rrii_simd2(int p, V2_BASETYPE *x, V2_BASETYPE *y)
{
	V2_INIT
	V2_BASETYPE *yr = y;
	V2_BASETYPE *xr = x;
	V2_BASETYPE *yi;
	V2_BASETYPE *xi;
	int n,k;
	
	for (n=0; n<=p; n++) {
		V2_SETITEM(yr, 0, V2_GETITEM(xr, 0));
		yr += 2;
		xr += 2;
		yi = yr + 2*n;
		xi = xr + 2*n;
		for (k=5; k<=n; k+=4) {
			V2_SETITEM(yr, 0, V2_NEG(V2_GETITEM(xi, 0)));
			V2_SETITEM(yi, 0, V2_GETITEM(xr, 0));
			V2_SETITEM(yr, 2, V2_NEG(V2_GETITEM(xr, 2)));
			V2_SETITEM(yi, 2, V2_NEG(V2_GETITEM(xi, 2)));
			V2_SETITEM(yr, 4, V2_GETITEM(xi, 4));
			V2_SETITEM(yi, 4, V2_NEG(V2_GETITEM(xr, 4)));
			V2_SETITEM(yr, 6, V2_GETITEM(xr, 6));
			V2_SETITEM(yi, 6, V2_GETITEM(xi, 6));
			yr += 8;
			xr += 8;
			yi += 8;
			xi += 8;
		}
		switch (n+4-k) {/* Note: no break statements! */
		case 3:	
			V2_SETITEM(yr, 6, V2_GETITEM(xr, 6));
			V2_SETITEM(yi, 6, V2_GETITEM(xi, 6));
		case 2:	
			V2_SETITEM(yr, 4, V2_GETITEM(xi, 4));
			V2_SETITEM(yi, 4, V2_NEG(V2_GETITEM(xr, 4)));
		case 1:	
			V2_SETITEM(yr, 2, V2_NEG(V2_GETITEM(xr, 2)));
			V2_SETITEM(yi, 2, V2_NEG(V2_GETITEM(xi, 2)));
		case 0:
			V2_SETITEM(yr, 0, V2_NEG(V2_GETITEM(xi, 0)));
			V2_SETITEM(yi, 0, V2_GETITEM(xr, 0));
			
			yi += 2*(n+5-k);
			xi += 2*(n+5-k);
		}	
		yr = yi;
		xr = xi;
	}	
}	


void Rz_minus_pi2_rrii_simd2(int p, V2_BASETYPE *x, V2_BASETYPE *y){
	V2_INIT
	V2_BASETYPE *yr = y;
	V2_BASETYPE *xr = x;
	V2_BASETYPE *yi;
	V2_BASETYPE *xi;
	
	int n,k;
	
	for (n=0; n<=p; n++) {
		V2_SETITEM(yr, 0, V2_GETITEM(xr, 0));
		yr += 2;
		xr += 2;
		yi = yr + 2*n;
		xi = xr + 2*n;
		for (k=5; k<=n; k+=4) {
			V2_SETITEM(yr, 0, V2_GETITEM(xi, 0));
			V2_SETITEM(yi, 0, V2_NEG(V2_GETITEM(xr, 0)));
			V2_SETITEM(yr, 2, V2_NEG(V2_GETITEM(xr, 2)));
			V2_SETITEM(yi, 2, V2_NEG(V2_GETITEM(xi, 2)));
			V2_SETITEM(yr, 4, V2_NEG(V2_GETITEM(xi, 4)));
			V2_SETITEM(yi, 4, V2_GETITEM(xr, 4));
			V2_SETITEM(yr, 6, V2_GETITEM(xr, 6));
			V2_SETITEM(yi, 6, V2_GETITEM(xi, 6));
			yr += 8;
			xr += 8;
			yi += 8;
			xi += 8;
		}
		switch (n+4-k) {/* Note: no break statements! */
		case 3:	
			V2_SETITEM(yr, 6, V2_GETITEM(xr, 6));
			V2_SETITEM(yi, 6, V2_GETITEM(xi, 6));
		case 2:	
			V2_SETITEM(yr, 4, V2_NEG(V2_GETITEM(xi, 4)));
			V2_SETITEM(yi, 4, V2_GETITEM(xr, 4));
		case 1:	
			V2_SETITEM(yr, 2, V2_NEG(V2_GETITEM(xr, 2)));
			V2_SETITEM(yi, 2, V2_NEG(V2_GETITEM(xi, 2)));
		case 0:
			V2_SETITEM(yr, 0, V2_GETITEM(xi, 0));
			V2_SETITEM(yi, 0, V2_NEG(V2_GETITEM(xr, 0)));

			yi += 2*(n+5-k);
			xi += 2*(n+5-k);
		}	
		yr = yi;
		xr = xi;
	}	
}	

			
