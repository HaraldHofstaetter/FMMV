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

	void* fmmvHandle;

	int typeSources = 0;
	int computeGradient = 0;
	int getStatistics = 0;
	
	int NParticles;
	int NTargets;

	double (*particles)[DIM];
	double (*targets)[DIM];
	
	
	if (nrhs<1) {
		mexErrMsgTxt("Not enough input arguments.");
	}

	if (nrhs>3) {
		mexErrMsgTxt("Too many output arguments.");
	}

	if (nlhs<1) {
		mexErrMsgTxt("Not enough output arguments.");
	}	
	
	if (nlhs>1) {
		mexErrMsgTxt("Too many output arguments.");
	}	
	
	if (!mxIsDouble(prhs[0])
	     ||mxIsComplex(prhs[0])
	     ||(mxGetM(prhs[0])!=DIM)) {
#if (DIM==2)
	     mexErrMsgTxt("First argument 'particles' must be a noncomplex matrix with 2 rows.");
#else
	     mexErrMsgTxt("First argument 'particles' must be a noncomplex matrix with 3 rows.");
#endif
	}	
	NParticles = mxGetN(prhs[0]);
	particles = (double(*)[DIM]) mxGetPr(prhs[0]);	

	
	if ((nrhs<=1)||(mxIsDouble(prhs[1])&&(mxGetN(prhs[1])==0)&&(mxGetM(prhs[1])==0))){
		NTargets = NParticles;
		targets = 0; /* no independent target locations */
	}	
	else if (!mxIsDouble(prhs[1])
	     ||mxIsComplex(prhs[1])
	     ||(mxGetM(prhs[1])!=DIM)) {
#if (DIM==2)
	     mexErrMsgTxt("Second argument 'targets' must be a noncomplex matrix with 2 rows.");
#else
	     mexErrMsgTxt("Second argument 'targets' must be a noncomplex matrix with 3 rows.");
#endif
	}
	else {
		NTargets = mxGetN(prhs[1]);
		targets =  (double(*)[DIM]) mxGetPr(prhs[1]);
	}

	options = fmmvGetDefaultOptions();
	if ((nrhs<=2)||(mxIsDouble(prhs[2])&&(mxGetN(prhs[2])==0)&&(mxGetM(prhs[2])==0))){
		/* default options */
	}
	else if (!mxIsStruct(prhs[2])||mxGetNumberOfElements(prhs[2])!=1) {
		mexErrMsgTxt("Fifth argument 'options' must be an options structure.");
	}
	else {
		int ifield;
		char *fname;
		double value;
		mxArray *tmp;
		int nfields = mxGetNumberOfFields(prhs[2]);

		for (ifield = 0; ifield < nfields; ifield++) {
			fname = mxGetFieldNameByNumber(prhs[2], ifield);
			tmp = mxGetFieldByNumber(prhs[2], 0, ifield);
			if ((tmp==0)|| 
			 (mxIsDouble(tmp)&&(mxGetN(tmp)==0)&&(mxGetM(tmp)==0))){
				//mexPrintf("%s: EMPTY\n", fname);
				value = -1;
			}
			else if (mxIsChar(tmp)) {
				char svalue[2];
				mxGetString(tmp, svalue, 4);
				//mexPrintf("%s: %s\n", fname, svalue);
				if (!strcmp(svalue, "on")) {
					value = 1;
				}	
				else if (!strcmp(svalue, "off")) {
					value = 0;
				}	
				else {
			     		mexPrintf("FMMV option: '%s'\n", fname);	
				        mexErrMsgTxt("Option must be a non-complex scalar.");
				}	
			}
			else if (!mxIsDouble(tmp)
	     			||mxIsComplex(tmp)
		  	        ||(mxGetM(tmp)!=1)
 	                        ||(mxGetN(tmp)!=1)) {
			     mexPrintf("FMMV option: '%s'\n", fname);	
			     mexErrMsgTxt("Option must be a non-complex scalar.");
			}
			else {
			     value = mxGetScalar(tmp);	

			}
			//mexPrintf("OPTION: %s VALUE: %g\n", fname, value);	

			     if (!strcmp(fname, "pM"))
			     	options.pM = value;
			else if (!strcmp(fname, "pL"))
			     	options.pL = value;
			else if (!strcmp(fname, "s"))
			     	options.s = value;
			else if (!strcmp(fname, "ws"))
			     	options.ws = value;
			else if (!strcmp(fname, "reducedScheme"))
			     	options.ws = value;
			else if (!strcmp(fname, "scale"))
			     	options.scale = value;
			else if (!strcmp(fname, "beta"))
			     	options.beta = value;
			else if (!strcmp(fname, "splitThreshold"))
			     	options.splitThreshold = value;
			else if (!strcmp(fname, "splitTargetThreshold" ))
			     	options.splitTargetThreshold = value;
			else if (!strcmp(fname, "levels"))
			     	options.levels = value;
			else if (!strcmp(fname, "directEvalThreshold"))
			     	options.directEvalThreshold = value;
			else if (!strcmp(fname, "periodicBoundaryConditions"))
			     	options.periodicBoundaryConditions = value;
			else if (!strcmp(fname, "extrinsicCorrection"))
			     	options.extrinsicCorrection = value;
			else if (!strcmp(fname, "useHilbertOrder"))
			     	options.useHilbertOrder = value;
			else if (!strcmp(fname, "directEvalAccuracy"))
			     	options.directEvalAccuracy= value;
			else if (!strcmp(fname, "useFarfieldNearfieldThreads"))
			     	options.useFarfieldNearfieldThreads = value;
			else if (!strcmp(fname, "getStatistics"))
				getStatistics = 1;
			else  {
			     mexPrintf("FMMV option: '%s'\n", fname);	
			     mexWarnMsgTxt("Unknown option.");
			}
		}	
	}


	if (getStatistics) {
		statistics_p = &statistics;
	}
	else {
		statistics_p = 0;
	}	

#if (DIM==2)
        fmmv2d_initialize
#else
        fmmv3d_initialize
#endif
                (
		&fmmvHandle,
		NParticles, 
		particles,
		NTargets,		
		targets,
		typeSources,
		computeGradient,
	        &options,
		statistics_p,
                err
		);

	if (err) {
		mexPrintf("%s\n", err);
		mexErrMsgTxt("FMMV error,");
	}

	{
	int dims[2]= { 1, 7};
	plhs[0] = mxCreateNumericArray(2, dims, mxUINT32_CLASS, mxREAL);
	}
	*((unsigned int*)mxGetData(plhs[0])) = 2718281828U; /* magic number */
	*((unsigned int*)mxGetData(plhs[0])+1) = (unsigned int) fmmvHandle;
	*((unsigned int*)mxGetData(plhs[0])+2) = (unsigned int) NParticles;
	*((unsigned int*)mxGetData(plhs[0])+3) = (unsigned int) NTargets;
	*((unsigned int*)mxGetData(plhs[0])+4) = (unsigned int) typeSources;
	*((unsigned int*)mxGetData(plhs[0])+5) = (unsigned int) computeGradient;
	*((unsigned int*)mxGetData(plhs[0])+6) = (unsigned int) getStatistics;

	if(getStatistics) {
		genStatistics(statistics_p);
	}	
}	

       
