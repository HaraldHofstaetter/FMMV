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

/* Potential and gradients computed by the direct method.
 * Accumulated sums are computed using extended (double-double)
 * precision.
 */

#include"_fmmv.h"
#include<math.h>

static __inline__ 
void acc(double *hi, double *lo, double b)
{
	double s1 = *hi + b;
	double h = s1 - *hi;
	double s2 = ((*hi - (s1-h)) + (b-h)) + *lo;
	h = s1 + s2;
	*lo = s2 - (h-s1);
	*hi = h;
}


static void exactPotential( unsigned int N, 
	_FLOAT_ particles[][2], _FLOAT_ charges[][2],
	_FLOAT_ targets0[][2],
	unsigned int index, _FLOAT_ pot[])
{
	int j;
	double x, y;
	double hx, hy;
	double qx, qy;
	double potx_hi, potx_lo;
	double poty_hi, poty_lo;
	_FLOAT_ (*targets)[2] = (!targets0 ? particles : targets0);
	
	potx_hi = 0; potx_lo = 0;
	poty_hi = 0; poty_lo = 0;
        for (j=0; j<N; j++) {
   	        if ((!targets0) && (j==index)) continue;
	        x = ((double)targets[index][0]) - ((double)particles[j][0]);
	        y = ((double)targets[index][1]) - ((double)particles[j][1]);
                hx = 0.5 * log(x*x + y*y);
                hy = atan2(y, x);
                qx = (double) charges[j][0];
                qy = (double) charges[j][1];
	        acc(&potx_hi, &potx_lo, hx*qx - hy*qy); 
	        acc(&poty_hi, &poty_lo, hx*qy + hy*qx); 
        }
	pot[0] = (_FLOAT_) potx_hi;
	pot[1] = (_FLOAT_) poty_hi;
}
	

static void exactDipolePotential( unsigned int N, 
	_FLOAT_ particles[][2], _FLOAT_ charges[][2],
	_FLOAT_ dipoleMoments[][2],
	_FLOAT_ targets0[][2],
	unsigned int index, _FLOAT_ pot[])
{
	int j;
	double x, y;
	double hx, hy;
	double qx, qy;
	double wx, wy;
	double potx_hi, potx_lo;
	double poty_hi, poty_lo;
	double r2, one_over_r2;

	_FLOAT_ (*targets)[2] = (!targets0 ? particles : targets0);
	
	potx_hi = 0; potx_lo = 0;
	poty_hi = 0; poty_lo = 0;
            
    	for (j=0; j<N; j++) {
	        if ((!targets0) && (j==index)) continue;
	        x = ((double)targets[index][0]) - ((double)particles[j][0]);
	        y = ((double)targets[index][1]) - ((double)particles[j][1]);
    
	        r2 = x*x + y*y;
                one_over_r2 = 1.0/r2;
	        wx = ( x*((double)dipoleMoments[j][0])
   	              +y*((double)dipoleMoments[j][1]))*one_over_r2;
	        wx = ( x*((double)dipoleMoments[j][1])
   	              -y*((double)dipoleMoments[j][0]))*one_over_r2;
                acc(&potx_hi, &potx_lo, wx);
                acc(&poty_hi, &poty_lo, wy);
                if(charges) {
                    hx = 0.5 * log(r2);
                    hy = atan2(y, x);
                    qx = (double) charges[j][0];
                    qy = (double) charges[j][1];
                    acc(&potx_hi, &potx_lo, hx*qx - hy*qy);
                    acc(&poty_hi, &poty_lo, hx*qy + hy*qx);
                }
        }
	pot[0] = (_FLOAT_) potx_hi;
	pot[1] = (_FLOAT_) poty_hi;
}	


static void exactGradient( unsigned int N, 
	_FLOAT_ particles[][2], _FLOAT_ charges[][2],
	_FLOAT_ targets0[][2],
	unsigned int index, _FLOAT_ pot[], _FLOAT_ grad[])
{
	int j;
	double x, y;
	double hx, hy;
	double qx, qy;
	double ux, uy;
	double potx_hi, potx_lo;
	double poty_hi, poty_lo;
	double gradx_hi, gradx_lo;
	double grady_hi, grady_lo;
        double r2, one_over_r2;
	_FLOAT_ (*targets)[2] = (!targets0 ? particles : targets0);
	
	potx_hi = 0; potx_lo = 0;
	poty_hi = 0; poty_lo = 0;
	gradx_hi = 0.0; gradx_lo = 0.0;
	grady_hi = 0.0; grady_lo = 0.0;


   	for (j=0; j<N; j++) {
	        if ((!targets0) && (j==index)) continue;
	        x = ((double)targets[index][0]) - ((double)particles[j][0]);
	        y = ((double)targets[index][1]) - ((double)particles[j][1]);
	        r2 = x*x + y*y;
                hx = 0.5 * log(r2);
                hy = atan2(y, x);
                qx = (double) charges[j][0];
                qy = (double) charges[j][1];
	        acc(&potx_hi, &potx_lo, hx*qx - hy*qy); 
	        acc(&poty_hi, &poty_lo, hx*qy + hy*qx); 

                one_over_r2 = 1.0/r2;
	        ux = x*one_over_r2;
                uy = y*one_over_r2;
	        acc(&gradx_hi, &gradx_lo, ux*qx + uy*qy);
	        acc(&grady_hi, &grady_lo, ux*qy - uy*qx);
        }
	pot[0] = (_FLOAT_) potx_hi;
	pot[1] = (_FLOAT_) poty_hi;
	grad[0] = (_FLOAT_) gradx_hi;
	grad[1] = (_FLOAT_) grady_hi;
}
	

static void exactDipoleGradient( unsigned int N, 
	_FLOAT_ particles[][2], _FLOAT_ charges[][2],
	_FLOAT_ dipoleMoments[][2],
	_FLOAT_ targets0[][2],
	unsigned int index, _FLOAT_ pot[], _FLOAT_ grad[])
{
	int j;
	double x, y;
	double potx_hi, potx_lo;
	double poty_hi, poty_lo;
	double hx, hy;
	double qx, qy;
	double wx, wy;
	double ux, uy;
	double gradx_hi, gradx_lo;
	double grady_hi, grady_lo;
        double r2, one_over_r2;
	_FLOAT_ (*targets)[2] = (!targets0 ? particles : targets0);
	
	potx_hi = 0; potx_lo = 0;
	poty_hi = 0; poty_lo = 0;
	gradx_hi = 0.0; gradx_lo = 0.0;
	grady_hi = 0.0; grady_lo = 0.0;

	for (j=0; j<N; j++) {
	        if ((!targets0) && (j==index)) continue;
	        x = ((double)targets[index][0]) - ((double)particles[j][0]);
	        y = ((double)targets[index][1]) - ((double)particles[j][1]);
        
	        r2 = x*x + y*y;
                one_over_r2 = 1.0/r2;
	        ux = x*one_over_r2;
                uy = y*one_over_r2;
	        wx = ( ux*((double)dipoleMoments[j][0])
   	              +uy*((double)dipoleMoments[j][1]));
	        wx = ( ux*((double)dipoleMoments[j][1])
   	              -uy*((double)dipoleMoments[j][0]));
                acc(&potx_hi, &potx_lo, wx);
                acc(&poty_hi, &poty_lo, wy);
                acc(&gradx_hi, &gradx_lo, -2.0*(x*wx + y*wy)*one_over_r2);
                acc(&grady_hi, &grady_lo, -2.0*(x*wy - y*wx)*one_over_r2);
                if(charges) {
                    hx = 0.5 * log(r2);
                    hy = atan2(y, x);
                    qx = (double) charges[j][0];
                    qy = (double) charges[j][1];
                    acc(&potx_hi, &potx_lo, hx*qx - hy*qy);
                    acc(&poty_hi, &poty_lo, hx*qy + hy*qx);
	            acc(&gradx_hi, &gradx_lo, ux*qx + uy*qy);
        	    acc(&grady_hi, &grady_lo, ux*qy - uy*qx);
                }
        }
	pot[0] = (_FLOAT_) potx_hi;
	pot[1] = (_FLOAT_) poty_hi;
	grad[0] = (_FLOAT_) gradx_hi;
	grad[1] = (_FLOAT_) grady_hi;
}	





#if (FMM_PRECISION == 0)
#	define errbuf _fmmvf_errbuf_
#else
#	define errbuf _fmmv_errbuf_
#endif

extern char errbuf[1024];

void direct_method_complex_xp
(
        unsigned int NParticles,
        _FLOAT_ particles[][2], 
	_FLOAT_ charges[][2],
	_FLOAT_ dipoleMoments[][2],
        unsigned int NTargets,
	_FLOAT_ targets[][2],
	_FLOAT_ pot[][2],
	_FLOAT_ grad[][2],
        double beta, 
        double *time,
        char **errorMessage
)
{
	int i;
        if (pot&&grad) {
                if (dipoleMoments) {
			for (i=0; i<NTargets; i++) {
        			exactDipoleGradient(NParticles, particles, charges, dipoleMoments, targets, i, pot[i], grad[i]); 
			}
                }
                else if(charges) {
			for (i=0; i<NTargets; i++) {
        			 exactGradient(NParticles, particles, charges, targets, i, pot[i], grad[i]); 
			}
                }
                else {
                        *errorMessage = "charges and dipoleMoments must not both be NULL.";
                        return;
                }
        }
        else if (pot) {
                if (dipoleMoments) {
			for (i=0; i<NTargets; i++) {
        			exactDipolePotential(NParticles, particles, charges, dipoleMoments, targets, i, pot[i]); 
			}
                }
                else if(charges) {
			for (i=0; i<NTargets; i++) {
        			exactPotential(NParticles, particles, charges, targets, i, pot[i]); 
			}
                }
                else {
                        *errorMessage = "charges and dipoleMoments must not both be NULL.";
                        return;
                }
        }
        else {
                *errorMessage = "pot must not be NULL.";
                return;
        }
	*errorMessage = 0;
}		

