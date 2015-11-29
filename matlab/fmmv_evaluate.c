/*
 * FMMV - the Fastest Multipole Method of Vienna
 * Copyright (c) 2006-2010 Harald Hofstaetter
 * Institute of Mathematics 
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

#include "mex.h"
#if (DIM==2)
#include "fmmv2d.h"
#else
#include "fmmv3d.h"
#endif

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{		
	char *err=0;
	struct FmmvOptions options;
	struct FmmvStatistics statistics;
	struct FmmvStatistics *statistics_p;	

	
	int NParticles;
	int NTargets;
	int computeGradient;
	int typeSources;
	int getStatistics;

	double *charges;
	double (*dipoleMoments)[DIM];
	double *pot;
	double (*grad)[DIM];

	void * fmmvHandle;
	
	
	if (nrhs<2) {
		mexErrMsgTxt("Not enough input arguments.");
	}

	if (nrhs>3) {
		mexErrMsgTxt("Too many input arguments.");
	}

	if (nlhs<1) {
		mexErrMsgTxt("Not enough output arguments.");
	}	
	
	if (nlhs>2) {
		mexErrMsgTxt("Too many output arguments.");
	}	
	
	if ((mxGetClassID(prhs[0])!=mxUINT32_CLASS)
	     ||mxIsComplex(prhs[0])
	     ||(mxGetM(prhs[0])!=1) 
	     ||(mxGetN(prhs[0])!=7)
	     ||(*((unsigned int*) mxGetPr(prhs[0]))!=2718281828U)) {
	     mexErrMsgTxt("First argument must be a fmmv handle.");
	}	
	fmmvHandle = (void *) *((unsigned int*) mxGetPr(prhs[0])+1);
	NParticles = *((unsigned int*) mxGetPr(prhs[0])+2);
	NTargets = *((unsigned int*) mxGetPr(prhs[0])+3);
	typeSources = *((unsigned int*) mxGetPr(prhs[0])+4);
	computeGradient = *((unsigned int*) mxGetPr(prhs[0])+5);
	getStatistics = *((unsigned int*) mxGetPr(prhs[0])+6);

	if (typeSources==1) {
		charges = 0; /* only dipole charges */
	}	
	else if	(!(  mxIsDouble(prhs[1])
		   &&!mxIsComplex(prhs[1])
		   &&(  ((mxGetN(prhs[1])==1)&&(mxGetM(prhs[1])==NParticles))
		      ||((mxGetM(prhs[1])==1)&&(mxGetN(prhs[1])==NParticles)))
		  )) {
	     mexErrMsgTxt("Second argument 'charges' must be a noncomplex vector with length = nparticles.");
	}
	else {
		charges = mxGetPr(prhs[1]);
	}

	if (typeSources==0) {
		dipoleMoments = 0; /* only monopole charges */
	}	
	else if ((nrhs<=2) 
	     ||(!mxIsDouble(prhs[2]))
	     ||mxIsComplex(prhs[2])
	     ||(mxGetM(prhs[2])!=DIM)
	     ||(mxGetN(prhs[2])!=NParticles)) {
#if (DIM==2)
	     mexErrMsgTxt("Third argument 'dipolemoments' must be a noncomplex matrix with 2 rows and nparticles columns.");
#else
	     mexErrMsgTxt("Third argument 'dipolemoments' must be a noncomplex matrix with 3 rows and nparticles columns.");
#endif
	}
	else {
		dipoleMoments = (double(*)[DIM]) mxGetPr(prhs[2]);
	}	

	if (computeGradient&&(nlhs<2)) {
	     mexErrMsgTxt("No output argument 'grad' specified.");
	}

	if (!computeGradient&&(nlhs>1)) {
	     mexWarnMsgTxt("Output argument 'grad' not assigned.");
	}
	

	if (getStatistics) {
		statistics_p = &statistics;
	}
	else {
		statistics_p = 0;
	}	

	plhs[0] = mxCreateDoubleMatrix(1, NTargets, mxREAL);
	pot = mxGetPr(plhs[0]);

	if (computeGradient) {
		plhs[1] = mxCreateDoubleMatrix(DIM, NTargets, mxREAL);
		grad =  (double(*)[DIM]) mxGetPr(plhs[1]);
	}	
	else {
		grad = 0;
	}

#if (DIM==2)
        fmmv2d_evaluate
#else
        fmmv3d_evaluate
#endif
                (
		fmmvHandle, 
		charges,  
		dipoleMoments,
		pot,
		grad,
		statistics_p,
                err
		);

	if (err) {
		mexPrintf("%s\n", err);
		mexErrMsgTxt("FMMV error,");
	}

	if(getStatistics) {
		genStatistics(statistics_p);
	}	
}	

       
