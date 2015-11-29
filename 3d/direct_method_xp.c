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
	_FLOAT_ particles[][3], _FLOAT_ charges[],
	_FLOAT_ targets0[][3],
	unsigned int index)
{
	int j;
	double x, y, z;
	double pot_hi, pot_lo;
	_FLOAT_ (*targets)[3] = (!targets0 ? particles : targets0);
	
	pot_hi = 0.0;
	pot_lo = 0.0;
        if (beta==0) {
	    for (j=0; j<N; j++) {
		if ((!targets0) && (j==index)) continue;
		x = ((double)targets[index][0]) - ((double)particles[j][0]);
		y = ((double)targets[index][1]) - ((double)particles[j][1]);
		z = ((double)targets[index][2]) - ((double)particles[j][2]);
		acc(&pot_hi, &pot_lo, ((double)charges[j])/sqrt(x*x + y*y + z*z)); 
	    }
        }
        else {
            double r;

	    for (j=0; j<N; j++) {
		if ((!targets0) && (j==index)) continue;
		x = ((double)targets[index][0]) - ((double)particles[j][0]);
		y = ((double)targets[index][1]) - ((double)particles[j][1]);
		z = ((double)targets[index][2]) - ((double)particles[j][2]);
                r = sqrt(x*x + y*y +z*z);
		acc(&pot_hi, &pot_lo, ((double)charges[j]) * exp(-beta*r)/r); 
	    }
        }
	return (_FLOAT_) pot_hi;
}
	

static _FLOAT_ exactDipolePotential(double beta, 
        unsigned int N, 
	_FLOAT_ particles[][3], _FLOAT_ charges[],
	_FLOAT_ dipoleMoments[][3],
	_FLOAT_ targets0[][3],
	unsigned int index)
{
	int j;
	double x, y, z;
	double pot_hi, pot_lo;
	_FLOAT_ (*targets)[3] = (!targets0 ? particles : targets0);
	
	pot_hi = 0.0;
	pot_lo = 0.0;
        if (beta==0) {
	    double one_over_r, m_times_rj, one_o_r_3;
            
    	    for (j=0; j<N; j++) {
		if ((!targets0) && (j==index)) continue;
		x = ((double)targets[index][0]) - ((double)particles[j][0]);
		y = ((double)targets[index][1]) - ((double)particles[j][1]);
		z = ((double)targets[index][2]) - ((double)particles[j][2]);

		one_over_r = 1.0/sqrt(x*x + y*y + z*z);
		one_o_r_3 = one_over_r*one_over_r*one_over_r;
		m_times_rj =  x*((double)dipoleMoments[j][0])
   		            + y*((double)dipoleMoments[j][1])
			    + z*((double)dipoleMoments[j][2]);
		acc(&pot_hi, &pot_lo, ((double)charges[j])*one_over_r + m_times_rj*one_o_r_3);
	    }
        }
        else {
            double r2, one_over_r, beta_r, exp_minus_beta_r_over_r, m_times_rj, hh0;

    	    for (j=0; j<N; j++) {
		if ((!targets0) && (j==index)) continue;
		x = ((double)targets[index][0]) - ((double)particles[j][0]);
		y = ((double)targets[index][1]) - ((double)particles[j][1]);
		z = ((double)targets[index][2]) - ((double)particles[j][2]);

                r2 = x*x + y*y + z*z;
		one_over_r = 1.0/sqrt(r2); 
                beta_r = beta*r2*one_over_r;
                exp_minus_beta_r_over_r = exp(-beta_r)*one_over_r;
		m_times_rj =  x*((double)dipoleMoments[j][0])
   		            + y*((double)dipoleMoments[j][1])
			    + z*((double)dipoleMoments[j][2]);
                hh0 = (1.0+beta_r)*one_over_r*one_over_r*exp_minus_beta_r_over_r;
		acc(&pot_hi, &pot_lo, ((double)charges[j])*exp_minus_beta_r_over_r + m_times_rj*hh0);
	    }
        }
	return pot_hi;
}	


static _FLOAT_ exactGradient(double beta, 
        unsigned int N, 
	_FLOAT_ particles[][3], _FLOAT_ charges[],
	_FLOAT_ targets0[][3],
	unsigned int index, _FLOAT_ grad[])
{
	int j;
	double x, y, z;
	double pot_hi, pot_lo;
	double gradx_hi, gradx_lo;
	double grady_hi, grady_lo;
	double gradz_hi, gradz_lo;
	_FLOAT_ (*targets)[3] = (!targets0 ? particles : targets0);
	
	pot_hi = 0; pot_lo = 0;
	gradx_hi = 0.0; gradx_lo = 0.0;
	grady_hi = 0.0; grady_lo = 0.0;
	gradz_hi = 0.0; gradz_lo = 0.0;

        if (beta==0) {
	    double one_over_r, q_one_over_r_3;

    	    for (j=0; j<N; j++) {
		if ((!targets0) && (j==index)) continue;
		x = ((double)targets[index][0]) - ((double)particles[j][0]);
		y = ((double)targets[index][1]) - ((double)particles[j][1]);
		z = ((double)targets[index][2]) - ((double)particles[j][2]);

		one_over_r = 1.0/sqrt(x*x + y*y + z*z);
		acc(&pot_hi, &pot_lo, ((double)charges[j])*one_over_r);
		q_one_over_r_3 = -((double)charges[j])*one_over_r*one_over_r*one_over_r;
		acc(&gradx_hi, &gradx_lo, x*q_one_over_r_3);
		acc(&grady_hi, &grady_lo, y*q_one_over_r_3);
		acc(&gradz_hi, &gradz_lo, z*q_one_over_r_3);
	    }
        }
        else {
            double r2, one_over_r, beta_r, exp_minus_beta_r_over_r, hh1;

    	    for (j=0; j<N; j++) {
		if ((!targets0) && (j==index)) continue;
		x = ((double)targets[index][0]) - ((double)particles[j][0]);
		y = ((double)targets[index][1]) - ((double)particles[j][1]);
		z = ((double)targets[index][2]) - ((double)particles[j][2]);

                r2 = x*x + y*y + z*z;
		one_over_r = 1.0/sqrt(r2); 
                beta_r = beta*r2*one_over_r;
                exp_minus_beta_r_over_r = exp(-beta_r)*one_over_r;
		acc(&pot_hi, &pot_lo, ((double)charges[j])*exp_minus_beta_r_over_r);
                hh1 = -((double)charges[j]) *
                      (1.0+beta_r)*one_over_r*one_over_r*exp_minus_beta_r_over_r;
		acc(&gradx_hi, &gradx_lo, x*hh1);
		acc(&grady_hi, &grady_lo, y*hh1);
		acc(&gradz_hi, &gradz_lo, z*hh1);
            }
        }
	grad[0] = (_FLOAT_) gradx_hi;
	grad[1] = (_FLOAT_) grady_hi;
	grad[2] = (_FLOAT_) gradz_hi;
	return (_FLOAT_) pot_hi;
}
	

static _FLOAT_ exactDipoleGradient(double beta, 
        unsigned int N, 
	_FLOAT_ particles[][3], _FLOAT_ charges[],
	_FLOAT_ dipoleMoments[][3],
	_FLOAT_ targets0[][3],
	unsigned int index, _FLOAT_ grad[])
{
	int j;
	double x, y, z;
	double pot_hi, pot_lo;
	double gradx_hi, gradx_lo;
	double grady_hi, grady_lo;
	double gradz_hi, gradz_lo;
	_FLOAT_ (*targets)[3] = (!targets0 ? particles : targets0);
	
	pot_hi = 0; pot_lo = 0;
	gradx_hi = 0.0; gradx_lo = 0.0;
	grady_hi = 0.0; grady_lo = 0.0;
	gradz_hi = 0.0; gradz_lo = 0.0;

        if (beta==0) {
	    double one_over_r, m_times_rj, one_o_r_3, qj_o_r_3;
	    for (j=0; j<N; j++) {
		if ((!targets0) && (j==index)) continue;
		x = ((double)targets[index][0]) - ((double)particles[j][0]);
		y = ((double)targets[index][1]) - ((double)particles[j][1]);
		z = ((double)targets[index][2]) - ((double)particles[j][2]);

		one_over_r = 1.0/sqrt(x*x + y*y + z*z);
		one_o_r_3 = one_over_r*one_over_r*one_over_r;
		m_times_rj =  x*((double)dipoleMoments[j][0])
  		            + y*((double)dipoleMoments[j][1])
			    + z*((double)dipoleMoments[j][2]);
	        qj_o_r_3 = one_o_r_3*(((double)charges[j]) + 3.0*m_times_rj*one_over_r*one_over_r);

		acc(&pot_hi, &pot_lo, ((double)charges[j])*one_over_r + m_times_rj*one_o_r_3);
		acc(&gradx_hi, &gradx_lo, ((double)dipoleMoments[j][0]) * one_o_r_3 - x*qj_o_r_3);
		acc(&grady_hi, &grady_lo, ((double)dipoleMoments[j][1]) * one_o_r_3 - y*qj_o_r_3);
		acc(&gradz_hi, &gradz_lo, ((double)dipoleMoments[j][2]) * one_o_r_3 - z*qj_o_r_3);
	    }
        }
        else {
            double beta2, r2, one_over_r, beta_r, exp_minus_beta_r_over_r, m_times_rj, hh0, hh1, hh2;

            beta2 = beta*beta;
	    for (j=0; j<N; j++) {
		if ((!targets0) && (j==index)) continue;
		x = ((double)targets[index][0]) - ((double)particles[j][0]);
		y = ((double)targets[index][1]) - ((double)particles[j][1]);
		z = ((double)targets[index][2]) - ((double)particles[j][2]);

                r2 = x*x + y*y + z*z;
		one_over_r = 1.0/sqrt(r2); 
                beta_r = beta*r2*one_over_r;
                exp_minus_beta_r_over_r = EXP2(-beta_r)*one_over_r;
		m_times_rj =  x*((double)dipoleMoments[j][0])
  		            + y*((double)dipoleMoments[j][1])
			    + z*((double)dipoleMoments[j][2]);
                hh0 = (1.0+beta_r)*one_over_r*one_over_r*exp_minus_beta_r_over_r;
		acc(&pot_hi, &pot_lo, ((double)charges[j])*exp_minus_beta_r_over_r + m_times_rj*hh0);
                hh1 = one_over_r*one_over_r*m_times_rj;
                hh2 = -(((double)charges[j]) + 3.0*hh1)*hh0 + beta2*hh1*exp_minus_beta_r_over_r;
  	        acc(&gradx_hi, &gradx_lo, x*hh2 + ((double)dipoleMoments[j][0])*hh0);
	        acc(&grady_hi, &grady_lo, y*hh2 + ((double)dipoleMoments[j][1])*hh0);
	        acc(&gradz_hi, &gradz_lo, z*hh2 + ((double)dipoleMoments[j][2])*hh0);
	    }
        }
	grad[0] = (_FLOAT_) gradx_hi;
	grad[1] = (_FLOAT_) grady_hi;
	grad[2] = (_FLOAT_) gradz_hi;
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
        _FLOAT_ particles[][3], 
        _FLOAT_ charges[],
        _FLOAT_ dipoleMoments[][3],
        unsigned int NTargets,
	_FLOAT_ targets[][3],
        _FLOAT_ pot[], 
	_FLOAT_ grad[][3], /* result of length NTarget */
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

