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
#include<math.h>

#if (FMM_PRECISION==0)
   #define R R_single
   #define B B_single
#endif

extern _FLOAT_ R[];
extern _FLOAT_ B[];

void compute_spherical_harmonics_generic(int p, _FLOAT_ *Y, _FLOAT_ sin_phi, _FLOAT_ cos_phi, _FLOAT_ cos_theta)
{

	_FLOAT_	pmm;
	_FLOAT_	pm1;
	_FLOAT_	pm2;
	_FLOAT_	pml;
	_FLOAT_	c;
	_FLOAT_	s;
	_FLOAT_	h;
	_FLOAT_	alpha;
	_FLOAT_	beta;
	_FLOAT_	sqrt_1_minus_cos_theta_2 = sqrt(1.0 - cos_theta*cos_theta);

	int m, l, k, kk;

	pmm = 1.0;
	
	/* m==0: ***************/
	Y[0] = B[0] * pmm;
	Y[1] = 0.0;
	pm2 = pmm;
	pml = pmm*cos_theta;
	Y[2] = pml;
	Y[3] = 0.0 ;
	k = 1;
	for (l=2; l<=p; l++) {
		pm1 = pml;
		pml = R[l]*((2*l-1)*cos_theta*pm1 - (l-1)*pm2);
		pm2 = pm1;
		k += l;
		Y[2*k] = pml;
		Y[2*k+1] = 0.0;
	}

	/* m==1: ***************/
	pmm *= -sqrt_1_minus_cos_theta_2;
	s = sin_phi;
	c = cos_phi;
	alpha = 1-c;
	beta = s;
	h = B[2] * pmm;
	Y[4] = c*h;
	Y[5] = s*h;
	pm2 = pmm;
	pml = 3*pmm*cos_theta;
	h = B[4] * pml;
	Y[8] = c*h;
	Y[9] = s*h;
	k=4;
	for (l=3; l<=p; l++) {
		pm1 = pml;
		pml = R[l-1]*((2*l-1)*cos_theta*pm1 - l*pm2);
		pm2 = pm1;
		k += l;
		h = B[k] * pml;
		Y[2*k] = c*h;
		Y[2*k+1] = s*h;
	}

	/* 2 <=m <= p-1: ***************/
	kk = 1;
	for (m=2; m<p; m++) {
		pmm *= (1-2*m)*sqrt_1_minus_cos_theta_2;
		h = (alpha*c + beta*s);
		s = s - alpha*s + beta*c;
		c -= h;
		kk += m;
		k = kk + m;
		h = B[k] * pmm;
		Y[2*k] = c*h;
		Y[2*k+1] = s*h;
		pm2 = pmm;
		pml = (2*m+1)*pmm*cos_theta;
		k += m+1;
		h = B[k] * pml;
		Y[2*k] = c*h;
		Y[2*k+1] = s*h;
		for (l=m+2; l<=p; l++) {
			pm1 = pml;
			pml = R[l-m] * ((2*l-1)*cos_theta*pm1 - (l+m-1)*pm2);
			pm2 = pm1;
			k += l;
			h = B[k] * pml;
			Y[2*k] = c*h;
			Y[2*k+1] = s*h;
		}
	}

	/* m==p: ***************/
        pmm *= (1-2*p)*sqrt_1_minus_cos_theta_2;
        h = (alpha*c + beta*s);
        s = s - alpha*s + beta*c;
        c -= h;
	kk += p;
	k = kk + p;
        h = B[k] * pmm;
        Y[2*k] = c*h;
        Y[2*k+1] = s*h;
}


