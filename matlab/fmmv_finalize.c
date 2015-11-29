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

	void * fmmvHandle;
	
	
	if (nrhs<1) {
		mexErrMsgTxt("Not enough input arguments.");
	}

	if (nrhs>1) {
		mexErrMsgTxt("Too many input arguments.");
	}

	if (nlhs>1) {
		mexErrMsgTxt("No output arguments required.");
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

	if (getStatistics) {
		statistics_p = &statistics;
	}
	else {
		statistics_p = 0;
	}	


#if (DIM==2)
        fmmv2d_finalize
#else
        fmmv3d_finalize
#endif
                (
		fmmvHandle, 
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

       
