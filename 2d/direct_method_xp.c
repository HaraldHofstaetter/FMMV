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

double k0(double x);
double k1(double x);

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


static _FLOAT_ exactPotential(double beta,
        unsigned int N, 
	_FLOAT_ particles[][2], _FLOAT_ charges[],
	_FLOAT_ targets0[][2],
	unsigned int index)
{
	int j;
	double x, y;
	double pot_hi, pot_lo;
	_FLOAT_ (*targets)[2] = (!targets0 ? particles : targets0);
	
	pot_hi = 0.0;
	pot_lo = 0.0;
        if (beta==0) {
	    for (j=0; j<N; j++) {
   	        if ((!targets0) && (j==index)) continue;
	        x = ((double)targets[index][0]) - ((double)particles[j][0]);
	        y = ((double)targets[index][1]) - ((double)particles[j][1]);
	        acc(&pot_hi, &pot_lo, (0.5*(double)charges[j])*log(x*x + y*y)); 
	    }
        }
        else {
	    for (j=0; j<N; j++) {
   	        if ((!targets0) && (j==index)) continue;
	        x = ((double)targets[index][0]) - ((double)particles[j][0]);
	        y = ((double)targets[index][1]) - ((double)particles[j][1]);
	        acc(&pot_hi, &pot_lo, ((double)charges[j])*k0(beta*sqrt(x*x + y*y))); 
	    }
        }
	return (_FLOAT_) pot_hi;
}
	

static _FLOAT_ exactDipolePotential(double beta,
        unsigned int N, 
	_FLOAT_ particles[][2], _FLOAT_ charges[],
	_FLOAT_ dipoleMoments[][2],
	_FLOAT_ targets0[][2],
	unsigned int index)
{
	int j;
	double x, y;
	double pot_hi, pot_lo;
	_FLOAT_ (*targets)[2] = (!targets0 ? particles : targets0);
	
	pot_hi = 0.0;
	pot_lo = 0.0;
            
        if (beta==0) {
	    double r2, one_over_r2, m_times_r;
    	    for (j=0; j<N; j++) {
	        if ((!targets0) && (j==index)) continue;
	        x = ((double)targets[index][0]) - ((double)particles[j][0]);
	        y = ((double)targets[index][1]) - ((double)particles[j][1]);
    
	        r2 = x*x + y*y;
                one_over_r2 = 1.0/r2;
	        m_times_r =  x*((double)dipoleMoments[j][0])
   	                   + y*((double)dipoleMoments[j][1]);
                if(charges) {
	            acc(&pot_hi, &pot_lo, 0.5*((double)charges[j])*log(r2) 
                                     + m_times_r*one_over_r2);
                }
                else {
	            acc(&pot_hi, &pot_lo,  m_times_r*one_over_r2);
                }
            }
        }
        else {
	    double r, one_over_r, m_times_r;
	    for (j=0; j<N; j++) {
   	        if ((!targets0) && (j==index)) continue;
	        x = ((double)targets[index][0]) - ((double)particles[j][0]);
	        y = ((double)targets[index][1]) - ((double)particles[j][1]);
                
	        r = sqrt(x*x + y*y);
                one_over_r = 1.0/r;
	        m_times_r =  x*((double)dipoleMoments[j][0])
   	                   + y*((double)dipoleMoments[j][1]);
                if(charges) {
	            acc(&pot_hi, &pot_lo, ((double)charges[j])*k0(beta*r)
                                          -beta*m_times_r*k1(beta*r)*one_over_r); 
                }
                else {
	            acc(&pot_hi, &pot_lo, -beta*m_times_r*k1(beta*r)*one_over_r); 
                }

	    }
        }
	return pot_hi;
}	


static _FLOAT_ exactGradient(double beta,
        unsigned int N, 
	_FLOAT_ particles[][2], _FLOAT_ charges[],
	_FLOAT_ targets0[][2],
	unsigned int index, _FLOAT_ grad[])
{
	int j;
	double x, y;
	double pot_hi, pot_lo;
	double gradx_hi, gradx_lo;
	double grady_hi, grady_lo;
	_FLOAT_ (*targets)[2] = (!targets0 ? particles : targets0);
	
	pot_hi = 0; pot_lo = 0;
	gradx_hi = 0.0; gradx_lo = 0.0;
	grady_hi = 0.0; grady_lo = 0.0;


        if (beta==0) {
    	    for (j=0; j<N; j++) {
	        double r2, one_over_r2, q_over_r2;
	        if ((!targets0) && (j==index)) continue;
	        x = ((double)targets[index][0]) - ((double)particles[j][0]);
	        y = ((double)targets[index][1]) - ((double)particles[j][1]);
    
	        r2 = x*x + y*y;
                one_over_r2 = 1.0/r2;
	        acc(&pot_hi, &pot_lo, (0.5*(double)charges[j])*log(r2)); 
	        q_over_r2 = ((double)charges[j])*one_over_r2;
	        acc(&gradx_hi, &gradx_lo, x*q_over_r2);
	        acc(&grady_hi, &grady_lo, y*q_over_r2);
            }
        }
        else {
	    double r, q_beta_k1_over_r;
	    for (j=0; j<N; j++) {
   	        if ((!targets0) && (j==index)) continue;
	        x = ((double)targets[index][0]) - ((double)particles[j][0]);
	        y = ((double)targets[index][1]) - ((double)particles[j][1]);

	        r = sqrt(x*x + y*y);

	        acc(&pot_hi, &pot_lo, ((double)charges[j])*k0(beta*r)); 
	        q_beta_k1_over_r = ((double)charges[j])*beta*k1(beta*r)/r;
                acc(&gradx_hi, &gradx_lo, -x*q_beta_k1_over_r);
                acc(&grady_hi, &grady_lo, -y*q_beta_k1_over_r);
	    }
        }
	grad[0] = (_FLOAT_) gradx_hi;
	grad[1] = (_FLOAT_) grady_hi;
	return (_FLOAT_) pot_hi;
}
	

static _FLOAT_ exactDipoleGradient(double beta,
        unsigned int N, 
	_FLOAT_ particles[][2], _FLOAT_ charges[],
	_FLOAT_ dipoleMoments[][2],
	_FLOAT_ targets0[][2],
	unsigned int index, _FLOAT_ grad[])
{
	int j;
	double x, y;
	double pot_hi, pot_lo;
	double gradx_hi, gradx_lo;
	double grady_hi, grady_lo;
	_FLOAT_ (*targets)[2] = (!targets0 ? particles : targets0);
	
	pot_hi = 0; pot_lo = 0;
	gradx_hi = 0.0; gradx_lo = 0.0;
	grady_hi = 0.0; grady_lo = 0.0;

        if (beta==0) {
	    double r2, one_over_r2, q_over_r2, m_times_r, m_times_r_over_r2;
	    for (j=0; j<N; j++) {
	        if ((!targets0) && (j==index)) continue;
	        x = ((double)targets[index][0]) - ((double)particles[j][0]);
	        y = ((double)targets[index][1]) - ((double)particles[j][1]);
        
	        r2 = x*x + y*y;
                one_over_r2 = 1.0/r2;
	        m_times_r =  x*((double)dipoleMoments[j][0])
   	                   + y*((double)dipoleMoments[j][1]);
                m_times_r_over_r2 = m_times_r*one_over_r2;
                if(charges) {
	            acc(&pot_hi, &pot_lo, 0.5*((double)charges[j])*log(r2) 
                                     + m_times_r_over_r2);
        
	            q_over_r2 = ((double)charges[j])*one_over_r2;
	            acc(&gradx_hi, &gradx_lo, x*q_over_r2 + (((double)dipoleMoments[j][0]) 
                    -2.0*x*m_times_r_over_r2)*one_over_r2);
	            acc(&grady_hi, &grady_lo, y*q_over_r2 + (((double)dipoleMoments[j][1]) 
                    -2.0*y*m_times_r_over_r2)*one_over_r2);
                }
                else {
	            acc(&pot_hi, &pot_lo, m_times_r_over_r2);
        
	            acc(&gradx_hi, &gradx_lo, (((double)dipoleMoments[j][0]) 
                    -2.0*x*m_times_r_over_r2)*one_over_r2);
	            acc(&grady_hi, &grady_lo, (((double)dipoleMoments[j][1]) 
                    -2.0*y*m_times_r_over_r2)*one_over_r2);
                }
            }
        }
        else {
            double r, one_over_r, one_over_r2, mx, my, m_times_r;
            double beta2_m_r_k0_over_r2, beta_k1_over_r3;
            double beta_k1_over_r, q_beta_k1_over_r, xmy_m_ymx;
	    for (j=0; j<N; j++) {
   	        if ((!targets0) && (j==index)) continue;
	        x = ((double)targets[index][0]) - ((double)particles[j][0]);
	        y = ((double)targets[index][1]) - ((double)particles[j][1]);

	        r = sqrt(x*x + y*y);
                one_over_r = 1.0/r;
                one_over_r2 = one_over_r*one_over_r;
	        mx = (double)dipoleMoments[j][0];
   	        my = (double)dipoleMoments[j][1];
                m_times_r = x*mx + y*my;
                xmy_m_ymx =  my*x - mx*y;

                beta_k1_over_r = beta*k1(beta*r)*one_over_r;
                beta_k1_over_r3 = beta_k1_over_r*one_over_r2;
                beta2_m_r_k0_over_r2 = beta*beta*m_times_r*k0(beta*r)*one_over_r2;

                if(charges) {
	            acc(&pot_hi, &pot_lo, ((double)charges[j])*k0(beta*r)
                                          -m_times_r*beta_k1_over_r); 
	            q_beta_k1_over_r = ((double)charges[j])*beta_k1_over_r;
	            acc(&gradx_hi, &gradx_lo, x*(beta2_m_r_k0_over_r2-q_beta_k1_over_r)
                               + (x*m_times_r+y*xmy_m_ymx)*beta_k1_over_r3); 
	            acc(&grady_hi, &grady_lo, y*(beta2_m_r_k0_over_r2-q_beta_k1_over_r)
                               + (y*m_times_r-x*xmy_m_ymx)*beta_k1_over_r3); 
                }
                else {
	            acc(&pot_hi, &pot_lo, -m_times_r*beta_k1_over_r); 
	            acc(&gradx_hi, &gradx_lo, x*beta2_m_r_k0_over_r2
                               + (x*m_times_r+y*xmy_m_ymx)*beta_k1_over_r3); 
	            acc(&grady_hi, &grady_lo, y*beta2_m_r_k0_over_r2
                               + (y*m_times_r-x*xmy_m_ymx)*beta_k1_over_r3); 
                }

	    }
        }
	grad[0] = (_FLOAT_) gradx_hi;
	grad[1] = (_FLOAT_) grady_hi;
	return (_FLOAT_) pot_hi;
}	





#if (FMM_PRECISION == 0)
#	define errbuf _fmmvf_errbuf_
#else
#	define errbuf _fmmv_errbuf_
#endif

extern char errbuf[1024];

void direct_method_xp
(
        unsigned int NParticles,
        _FLOAT_ particles[][2], 
        _FLOAT_ charges[],
        _FLOAT_ dipoleMoments[][2],
        unsigned int NTargets,
	_FLOAT_ targets[][2],
        _FLOAT_ pot[], 
	_FLOAT_ grad[][2], /* result of length NTarget */
        double beta, 
        double *time,
        char **errorMessage
)
{
	int i;
        if (grad) {
                if (dipoleMoments) {
			for (i=0; i<NTargets; i++) {
        			pot[i] = exactDipoleGradient(beta, NParticles, particles, charges, dipoleMoments, targets, i, grad[i]); 
			}
                }
                else {
			for (i=0; i<NTargets; i++) {
        			pot[i] = exactGradient(beta, NParticles, particles, charges, targets, i, grad[i]); 
			}
                }
        }
        else {
                if (dipoleMoments) {
			for (i=0; i<NTargets; i++) {
        			pot[i] = exactDipolePotential(beta, NParticles, particles, charges, dipoleMoments, targets, i); 
			}
                }
                else {
			for (i=0; i<NTargets; i++) {
        			pot[i] = exactPotential(beta, NParticles, particles, charges, targets, i); 
			}
                }
        }
	*errorMessage = 0;
}		

