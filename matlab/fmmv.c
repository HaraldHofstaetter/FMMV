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

#include <stdint.h>  /* uint16_t etc. */
#include <uchar.h>  /* uchar16_t etc. */
#include "mex.h"
#if (DIM==2)
#include "fmmv2d.h"
#else
#include "fmmv3d.h"
#endif
//#include<fenv.h> 


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{		
	char *err=0;
	struct FmmvOptions options;
	struct FmmvStatistics statistics;
	struct FmmvStatistics *statistics_p;	

	int computeGradient = (nlhs>1);
	int getStatistics = 0;
	
	int NParticles;
	int NTargets;

	double (*particles)[DIM];
	double *charges;
	double (*dipoleMoments)[DIM];
	double (*targets)[DIM];
	double *pot;
	double (*grad)[DIM];

//  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );


	
	
	if (nrhs<2) {
		mexErrMsgTxt("Not enough input arguments.");
	}

	if (nrhs>5) {
		mexErrMsgTxt("Too many output arguments.");
	}

	if (nlhs<1) {
		mexErrMsgTxt("Not enough output arguments.");
	}	
	
	if (nlhs>2) {
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

	if (mxIsDouble(prhs[1])&&(mxGetN(prhs[1])==0)&&(mxGetM(prhs[1])==0)){
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

	if ((nrhs<=2)||(mxIsDouble(prhs[2])&&(mxGetN(prhs[2])==0)&&(mxGetM(prhs[2])==0))){
		dipoleMoments = 0; /* only monopole charges */
	}	
	else if (!mxIsDouble(prhs[2])
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
	
	if ((nrhs<=3)||(mxIsDouble(prhs[3])&&(mxGetN(prhs[3])==0)&&(mxGetM(prhs[3])==0))){
		NTargets = NParticles;
		targets = 0; /* no independent target locations */
	}	
	else if (!mxIsDouble(prhs[3])
	     ||mxIsComplex(prhs[3])
	     ||(mxGetM(prhs[3])!=DIM)) {
#if (DIM==2)
	     mexErrMsgTxt("Fourth argument 'targets' must be a noncomplex matrix with 2 rows.");
#else
	     mexErrMsgTxt("Fourth argument 'targets' must be a noncomplex matrix with 3 rows.");
#endif
	}
	else {
		NTargets = mxGetN(prhs[3]);
		targets =  (double(*)[DIM]) mxGetPr(prhs[3]);
	}

	options = fmmvGetDefaultOptions();
	if ((nrhs<=4)||(mxIsDouble(prhs[4])&&(mxGetN(prhs[4])==0)&&(mxGetM(prhs[4])==0))){
		/* default options */
	}
	else if (!mxIsStruct(prhs[4])||mxGetNumberOfElements(prhs[4])!=1) {
		mexErrMsgTxt("Fifth argument 'options' must be an options structure.");
	}
	else {
		int ifield;
		char *fname;
		double value;
		mxArray *tmp;
		int nfields = mxGetNumberOfFields(prhs[4]);

		for (ifield = 0; ifield < nfields; ifield++) {
			fname = mxGetFieldNameByNumber(prhs[4], ifield);
			tmp = mxGetFieldByNumber(prhs[4], 0, ifield);
			if ((tmp==0)|| 
			 (mxIsDouble(tmp)&&(mxGetN(tmp)==0)&&(mxGetM(tmp)==0))){
				//mexPrintf("%s: EMPTY\n", fname);
				value = -1;
			}
			else if (mxIsChar(tmp)) {
				char svalue[4];
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
			     	options.directEvalAccuracy = value;
			else if (!strcmp(fname, "useFarfieldNearfieldThreads"))
			     	options.useFarfieldNearfieldThreads = value;
			else if (!strcmp(fname, "getStatistics"))
				getStatistics = 1;
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
			     	options.directEvalAccuracy = value;
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
	fmmv2d
#else
	fmmv3d
#endif
                (
                NParticles, 
		particles,
		charges,  
		dipoleMoments,
		NTargets,		
		targets,
		pot,
		grad,
	        &options,
		statistics_p,
                &err
		);

	if (err) {
		mexPrintf("%s\n", err);
		mexErrMsgTxt("FMMV error,");
	}

	if(getStatistics) {
		genStatistics(statistics_p);
	}	
}	

       
