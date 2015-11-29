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

#include"fmmv3d.h"
#include<math.h>
#include<stdlib.h>
#include <time.h>
#include <sys/time.h>


#define NO_PARTICLES_B 100
#define MAX_NO_PARTICLES_PER_BOX 1000
#define NParticles  2*MAX_NO_PARTICLES_PER_BOX /* 2 boxes */
enum { STANDARD=0, DIPOLE=1, GRAD=2, DIPOLE_GRAD=3 };



static double get_elapsed_time(void)
{
        struct timeval buf;
        gettimeofday(&buf,0);
        return (double)buf.tv_sec
         + ((double)buf.tv_usec/1000000);
}


static void measure_time2(FmmvHandle *FMMV, void (*fun)(FmmvHandle*, Box*, Box*), Box *boxA, Box *boxB, int niter, int n, double *mean, double *min, double *max, double *stdder)
{
	double t;
	double tmin = 1e30; 
	double tmax = 0.0;
	double tmean = 0.0;
	double t2mean = 0.0;
	int i,j;

	for (j=0; j<n; j++) {
		if (niter==1) {
			t = get_elapsed_time();
			fun(FMMV, boxA, boxB);
			t = (get_elapsed_time()-t)/((double) niter);
		}
		else {
			t = get_elapsed_time();
			for (i=0; i<niter; i++) {
				fun(FMMV, boxA, boxB);
			}
			t = (get_elapsed_time()-t)/((double) niter);
		}	
		tmean += t;
		t2mean += t*t;
		if (t<tmin) tmin = t;
		if (t>tmax) tmax = t;
	}

	tmean /= (double) n;
	if (mean) *mean = tmean;
	if (min) *min = tmin;
	if (max) *max = tmax;
	if (stdder) *stdder = sqrt((t2mean - tmean*tmean)/(double) n);
}	
	
static void init_internal_data(int dipole_grad, int n)
{	
	switch (dipole_grad) {
	case STANDARD:
		init_internal_data_standard(n);
		break;
	case DIPOLE:
		init_internal_data_dipole(n);
		break;
	case GRAD:
		init_internal_data_grad(n);
		break;
	case DIPOLE_GRAD:	
		init_internal_data_dipole_grad(n);
		break;
	}	
}	

		
void gen_L_eval_M_grad_p16(Box*, Box*);

int determineDirectEvalThreshold(FmmvHandle *FMMV, int grad, int dipole, int approx_inv_sqrt)
{
	int dipole_grad, err, i;
	Box *boxA, *boxB;
	_FLOAT_ particles[NParticles][3];
	_FLOAT_ charges[NParticles];
	_FLOAT_ dipoleMoments[NParticles][3];
	int perm[NParticles];
	double tmean, tmin, tmax, tt;
	int nparts;

	for (i=0; i<NParticles; i++) {
		perm[i] = i;
	}	


	if (grad) {
		if (dipole) {
			dipole_grad = DIPOLE_GRAD;
		}
		else {
			dipole_grad = GRAD;
		}
	}
	else {
		if (dipole) {
			dipole_grad = DIPOLE;
		}
		else {
			dipole_grad = STANDARD;
		}
	}	

	if (dipole) {
		for (i=0; i<NParticles; i++) {
			charges[i] = 0.0;
			dipoleMoments[i][0] = 2.0*((rand() / (_FLOAT_) RAND_MAX)) - 1.0 ;
			dipoleMoments[i][1] = 2.0*((rand() / (_FLOAT_) RAND_MAX)) - 1.0 ;
			dipoleMoments[i][2] = 2.0*((rand() / (_FLOAT_) RAND_MAX)) - 1.0 ;
		}
	}
	else {
		for (i=0; i<NParticles; i++) {
			charges[i] = 2.0*((rand() / (_FLOAT_) RAND_MAX)) - 1.0 ;
		}	
	}	
	for (i=0; i<MAX_NO_PARTICLES_PER_BOX; i++) {
		particles[i][0] = 0.005 + 0.24*((rand() / (_FLOAT_) RAND_MAX));
		particles[i][1] = 0.005 + 0.24*((rand() / (_FLOAT_) RAND_MAX));
		particles[i][2] = 0.005 + 0.24*((rand() / (_FLOAT_) RAND_MAX));
	}	
	for (i=MAX_NO_PARTICLES_PER_BOX; i<MAX_NO_PARTICLES_PER_BOX; 2*i++) {
		particles[i][0] = 0.755 + 0.24*((rand() / (_FLOAT_) RAND_MAX));
		particles[i][1] = 0.755 + 0.24*((rand() / (_FLOAT_) RAND_MAX));
		particles[i][2] = 0.755 + 0.24*((rand() / (_FLOAT_) RAND_MAX));
	}	
	

	switch (dipole_grad) {
	case STANDARD:
		err = allocate_internal_data_standard(NParticles);
		copy_internal_data_standard(NParticles, perm, particles, charges);
		break;
	case DIPOLE:
		err = allocate_internal_data_dipole(NParticles);
		copy_internal_data_dipole(NParticles, perm, particles, charges, dipoleMoments);
		break;
	case GRAD:
		err = allocate_internal_data_grad(NParticles);
		copy_internal_data_grad(NParticles, perm, particles, charges);
		break;
	case DIPOLE_GRAD:	
		err = allocate_internal_data_dipole_grad(NParticles);
		copy_internal_data_dipole_grad(NParticles, perm, particles, charges, dipoleMoments);
		break;
	}	

	init_internal_data(GRAD, NParticles);

	FMMV->firstSourceBoxOfLevel[2] = 0;
	boxA = createBox(FMMV, 0, 0.875, 0.875, 0.875, MAX_NO_PARTICLES_PER_BOX+1, 2*MAX_NO_PARTICLES_PER_BOX , 2, 0, SOURCE|TARGET);
	boxB = createBox(FMMV, 0, 0.125, 0.125, 0.125, 1, NO_PARTICLES_B , 2, 0, SOURCE|TARGET);

	gen_M_grad_p6(boxA, boxB);
	measure_time2(FMMV, gen_L_eval_M_grad_p16, boxA, boxB, 10, 30, &tmean, &tmin, &tmax, 0);
	free(boxB->M); 
	printf("%4i %.4e %.4e %.4e\n", nparts, tmean, tmin, tmax);
	tt = tmin;
	printf("--------------------------------------------\n");
	for (nparts=10; nparts<=150; nparts+=5) {
		boxA->noOfParticles = nparts;
		measure_time2(FMMV, eval_direct_grad, boxA, boxB, 10, 20, &tmean, &tmin, &tmax, 0);
		printf("%4i %.4e %.4e %.4e\n", nparts, tmean, tmin, tmax);
		//if (tmin>1.0*tt) break;
	}	
			
	


	free(boxA);
	free(boxB);

	switch (dipole_grad) {
	case STANDARD:
		free_internal_data_standard(NParticles);
		break;
	case DIPOLE:
		free_internal_data_dipole(NParticles);
		break;
	case GRAD:
		free_internal_data_grad(NParticles);
		break;
	case DIPOLE_GRAD:	
		free_internal_data_dipole_grad(NParticles);
		break;
	}
}
