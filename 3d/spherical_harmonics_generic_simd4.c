/*
 * FMMV - the Fastest Multipole Method of Vienna
 * Copyright (c) 2006-2010 Harald Hofstaetter
 * Institute Mathematics
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

#include "_fmmv.h"
#include<math.h>
#include "simd.h"
#if (FMM_PRECISION==0)
   #include"simd4s.h"
#else
   #include"simd4d.h"
#endif


#if (FMM_PRECISION==0)
   #define R R_single
   #define B B_single
#endif

extern _FLOAT_ R[];
extern _FLOAT_ B[];

void compute_spherical_harmonics_generic_simd4(int p, V4_BASETYPE *Y, V4_TYPE sin_phi, V4_TYPE cos_phi, V4_TYPE cos_theta)
{
	V4_INIT

	V4_TYPE	pmm;
	V4_TYPE	pm1;
	V4_TYPE	pm2;
	V4_TYPE	pml;
	V4_TYPE	c;
	V4_TYPE	s;
	V4_TYPE	h;
	V4_TYPE	alpha;
	V4_TYPE	beta;
	V4_TYPE sqrt_1_minus_cos_theta_2 = V4_SQRT(V4_FNMS(cos_theta, cos_theta, V4_SET1(1.0)));

	int m, l, k, kk;

	//pmm = 1.0;
	pmm = V4_SET1(1.0);
	
	/* m==0: ***************/
	//Y[0] = B[0] * pmm;
	V4_SETITEM(Y, 0, V4_SET1(1.0000000000000000000E0));
	//Y[1] = 0.0;
	V4_SETITEM(Y, 4, V4_ZERO);
	pm2 = pmm;
	//pml = pmm*cos_theta;
	pml = V4_MUL(pmm, cos_theta);
	//Y[2] = pml;
	V4_SETITEM(Y, 8, pml);
	//Y[3] = 0.0 ;
	V4_SETITEM(Y, 12, V4_ZERO);
	k = 1;
	for (l=2; l<=p; l++) {
		pm1 = pml;
		//pml = R[l]*((2*l-1)*cos_theta*pm1 - (l-1)*pm2); 
		pml = V4_MUL(V4_GETBASEITEM(R, l), 
		      V4_SUB(V4_MUL(V4_SET1(2*l-1),V4_MUL(cos_theta, pm1)), V4_MUL(V4_SET1(l-1), pm2)));
		pm2 = pm1;
		k += l;
		//Y[2*k] = pml;
		V4_SETITEM(Y, 8*k, pml);
		//Y[2*k+1] = 0.0;
		V4_SETITEM(Y, 8*k+4, V4_ZERO);
	}

	/* m==1: ***************/
	//pmm *= -sqrt_1_minus_cos_theta_2;
	pmm = V4_MUL(pmm, V4_NEG(sqrt_1_minus_cos_theta_2));
	s = sin_phi;
	c = cos_phi;
	//alpha = 1-c;
	alpha = V4_SUB(V4_SET1(1), c);
	beta = s;
	//h = B[2] * pmm;
	h = V4_MUL(V4_GETBASEITEM(B, 2), pmm);
	//Y[4] = c*h;
	//Y[5] = s*h;
	V4_SETITEM(Y, 16, V4_MUL(c, h));
	V4_SETITEM(Y, 20, V4_MUL(s, h));
	pm2 = pmm;
	//pml = 3*pmm*cos_theta;
	pml = V4_MUL(V4_MUL(V4_SET1(3), pmm), cos_theta);
	//h = B[4] * pml;
	h = V4_MUL(V4_GETBASEITEM(B, 4), pml);
	//Y[8] = c*h;
	//Y[9] = s*h;
	V4_SETITEM(Y, 32, V4_MUL(c, h));
	V4_SETITEM(Y, 36, V4_MUL(s, h));
	k=4;
	for (l=3; l<=p; l++) {
		pm1 = pml;
		//pml = R[l-1]*((2*l-1)*cos_theta*pm1 - l*pm2); 
		pml = V4_MUL(V4_GETBASEITEM(R, l-1), 
		      V4_SUB(V4_MUL(V4_SET1(2*l-1),V4_MUL(cos_theta, pm1)), V4_MUL(V4_SET1(l), pm2)));
		pm2 = pm1;
		k += l;
		//h = B[k] * pml;
		h = V4_MUL(V4_GETBASEITEM(B, k), pml);
		//Y[2*k] = c*h;
		//Y[2*k+1] = s*h;
		V4_SETITEM(Y, 8*k, V4_MUL(c, h));
		V4_SETITEM(Y, 8*k+4, V4_MUL(s, h));
	}

	/* 2 <=m <= p-1: ***************/
	kk = 1;
	for (m=2; m<p; m++) {
		//pmm *= (1-2*m)*sqrt_1_minus_cos_theta_2;
		pmm = V4_MUL(pmm, V4_MUL(V4_SET1(1-2*m), sqrt_1_minus_cos_theta_2));
		//h = (alpha*c + beta*s);
		h = V4_FMA(alpha, c, V4_MUL(beta, s));
		//s = s - alpha*s + beta*c;
		s = V4_FMA(beta, c, V4_FNMS(alpha, s, s));
		//c -= h;
		c = V4_SUB(c, h);
		kk += m;
		k = kk + m;
		//h = B[k] * pmm;
		h = V4_MUL(V4_GETBASEITEM(B, k), pmm);
		//Y[2*k] = c*h;
		//Y[2*k+1] = s*h;
		V4_SETITEM(Y, 8*k, V4_MUL(c, h));
		V4_SETITEM(Y, 8*k+4, V4_MUL(s, h));
		pm2 = pmm;
		//pml = (2*m+1)*pmm*cos_theta;
		pml = V4_MUL(V4_MUL(V4_SET1(2*m+1), pmm), cos_theta);
		k += m+1;
		//h = B[k] * pml;
		h = V4_MUL(V4_GETBASEITEM(B, k), pml);
		//Y[2*k] = c*h;
		//Y[2*k+1] = s*h;
		V4_SETITEM(Y, 8*k, V4_MUL(c, h));
		V4_SETITEM(Y, 8*k+4, V4_MUL(s, h));
		for (l=m+2; l<=p; l++) {
			pm1 = pml;
			//pml = R[l-m] * ((2*l-1)*cos_theta*pm1 - (l+m-1)*pm2); 
			pml = V4_MUL(V4_GETBASEITEM(R, l-m), 
			      V4_SUB(V4_MUL(V4_SET1(2*l-1),V4_MUL(cos_theta, pm1)), V4_MUL(V4_SET1(l+m-1), pm2)));
			pm2 = pm1;
			k += l;
			//h = B[k] * pml;
			h = V4_MUL(V4_GETBASEITEM(B, k), pml);
			//Y[2*k] = c*h;
			//Y[2*k+1] = s*h;
			V4_SETITEM(Y, 8*k, V4_MUL(c, h));
			V4_SETITEM(Y, 8*k+4, V4_MUL(s, h));
		}
	}

	/* m==p: ***************/
        //pmm *= (1-2*p)*sqrt_1_minus_cos_theta_2;
	pmm = V4_MUL(pmm, V4_MUL(V4_SET1(1-2*p), sqrt_1_minus_cos_theta_2));
        //h = (alpha*c + beta*s);
	h = V4_FMA(alpha, c, V4_MUL(beta, s));
        //s = s - alpha*s + beta*c;
	s = V4_FMA(beta, c, V4_FNMS(alpha, s, s));
        //c -= h;
	c = V4_SUB(c, h);
	kk += p;
	k = kk + p;
        //h = B[k] * pmm;
	h = V4_MUL(V4_GETBASEITEM(B, k), pmm);
        //Y[2*k] = c*h;
        //Y[2*k+1] = s*h;
	V4_SETITEM(Y, 8*k, V4_MUL(c, h));
	V4_SETITEM(Y, 8*k+4, V4_MUL(s, h));
}


