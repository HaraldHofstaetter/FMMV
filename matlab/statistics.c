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

#include <stdint.h>  /* uint16_t etc. */
#include <uchar.h>  /* uchar16_t */
#include "mex.h"
#if (DIM==2)
#include "fmmv2d.h"
#else
#include "fmmv3d.h"
#endif

void genStatistics(struct FmmvStatistics *statistics)
{
	static char *stat_names[] = {
		"totalTime",
		"buildTreeTime",
		"genMTime",
		"M2MTime",
		"M2LTime",
		"L2LTime",
		"evalLTime",
		"list1Time",
		"list3Time",
		"list4Time",
		"list34Time",
		"farfieldTime",
		"nearfieldTime",
		"initializeTime",
		"evaluateTime",
		"finalizeTime",

		"noOfParticles",
		"noOfTargets",
		"maxNoOfStoredXin",
		"noOfSourceLevels",
		"noOfTargetLevels",
		"noOfSourceBoxes", 
		"noOfTargetBoxes",
		"noOfSourceLeafBoxes",
		"noOfTargetLeafBoxes",
		"averageNoOfParticlesPerLeafBox",
		"averageNoOfTargetsPerLeafBox",
		"maxAllocatedMemory",
		"noOfDirectInteractions",

		"noOfParticlesInLevel",
		"noOfTargetsInLevel",
		"noOfSourceBoxesInLevel",
		"noOfTargetBoxesInLevel",
		"noOfSourceLeafBoxesInLevel",
		"noOfTargetLeafBoxesInLevel",
		"averageNoOfParticlesPerLeafBoxInLevel",
		"averageNoOfTargetsPerLeafBoxInLevel",
	};
	
	int i;
	mxArray *stat = mxCreateStructMatrix(1, 1, 36, stat_names);
	mxArray *tmp;
	double *p;

	for (i=0; i<_STAT_LAST_; i++) {
		if (statistics->time[i]>=0.0) {
			mxSetFieldByNumber(stat, 0, i, mxCreateDoubleScalar(statistics->etime[i]));
		}	
	}	
	if (statistics->noOfParticles >=0)
		mxSetField(stat, 0, "noOfParticles", mxCreateDoubleScalar(statistics->noOfParticles));
	if (statistics->noOfTargets >=0)
		mxSetField(stat, 0, "noOfTargets", mxCreateDoubleScalar(statistics->noOfTargets));
	if (statistics->maxNoOfStoredXin >=0)
		mxSetField(stat, 0, "maxNoOfStoredXin", mxCreateDoubleScalar(statistics->maxNoOfStoredXin));
	if (statistics->noOfSourceLevels >=0)
		mxSetField(stat, 0, "noOfSourceLevels", mxCreateDoubleScalar(statistics->noOfSourceLevels));
	if (statistics->noOfTargetLevels >=0)
		mxSetField(stat, 0, "noOfTargetLevels", mxCreateDoubleScalar(statistics->noOfTargetLevels));
	if (statistics->noOfSourceBoxes >=0)
		mxSetField(stat, 0, "noOfSourceBoxes", mxCreateDoubleScalar(statistics->noOfSourceBoxes));
	if (statistics->noOfTargetBoxes >=0)
		mxSetField(stat, 0, "noOfTargetBoxes", mxCreateDoubleScalar(statistics->noOfTargetBoxes));
	if (statistics->noOfSourceLeafBoxes >=0)
		mxSetField(stat, 0, "noOfSourceLeafBoxes", mxCreateDoubleScalar(statistics->noOfSourceLeafBoxes));
	if (statistics->noOfTargetLeafBoxes >=0)
		mxSetField(stat, 0, "noOfTargetLeafBoxes", mxCreateDoubleScalar(statistics->noOfTargetLeafBoxes));
	if (statistics->averageNoOfParticlesPerLeafBox >=0)
		mxSetField(stat, 0, "averageNoOfParticlesPerLeafBox", mxCreateDoubleScalar(statistics->averageNoOfParticlesPerLeafBox));
	if (statistics->averageNoOfTargetsPerLeafBox >=0)
		mxSetField(stat, 0, "averageNoOfTargetsPerLeafBox", mxCreateDoubleScalar(statistics->averageNoOfTargetsPerLeafBox));
	if (statistics->maxAllocatedMemory >=0)
		mxSetField(stat, 0, "maxAllocatedMemory", mxCreateDoubleScalar(statistics->maxAllocatedMemory));
	if (statistics->noOfDirectInteractions >=0)
		mxSetField(stat, 0, "noOfDirectInteractions", mxCreateDoubleScalar(statistics->noOfDirectInteractions));

	if (statistics->noOfParticles >=0) {
		tmp = mxCreateDoubleMatrix(1, statistics->noOfSourceLevels+1, mxREAL);
		p = mxGetPr(tmp);
		for (i=0; i<=statistics->noOfSourceLevels; i++) 
			p[i] = statistics->noOfParticlesInLevel[i];
		mxSetField(stat, 0, "noOfParticlesInLevel", tmp);	
	}	
	if (statistics->noOfTargets >=0) {
		tmp = mxCreateDoubleMatrix(1, statistics->noOfTargetLevels+1, mxREAL);
		p = mxGetPr(tmp);
		for (i=0; i<=statistics->noOfTargetLevels; i++) 
			p[i] = statistics->noOfTargetsInLevel[i];
		mxSetField(stat, 0, "noOfTargetsInLevel", tmp);	
	}	
	if (statistics->noOfSourceBoxes >=0) {
		tmp = mxCreateDoubleMatrix(1, statistics->noOfSourceLevels+1, mxREAL);
		p = mxGetPr(tmp);
		for (i=0; i<=statistics->noOfSourceLevels; i++) 
			p[i] = statistics->noOfSourceBoxesInLevel[i];
		mxSetField(stat, 0, "noOfSourceBoxesInLevel", tmp);	
	}	
	if (statistics->noOfTargetBoxes >=0) {
		tmp = mxCreateDoubleMatrix(1, statistics->noOfTargetLevels+1, mxREAL);
		p = mxGetPr(tmp);
		for (i=0; i<=statistics->noOfTargetLevels; i++) 
			p[i] = statistics->noOfTargetBoxesInLevel[i];
		mxSetField(stat, 0, "noOfTargetBoxesInLevel", tmp);
	}	
	if (statistics->noOfSourceLeafBoxes >=0) {
		tmp = mxCreateDoubleMatrix(1, statistics->noOfSourceLevels+1, mxREAL);
		p = mxGetPr(tmp);
		for (i=0; i<=statistics->noOfSourceLevels; i++) 
			p[i] = statistics->noOfSourceLeafBoxesInLevel[i];
		mxSetField(stat, 0, "noOfSourceLeafBoxesInLevel", tmp);	
	}	
	if (statistics->noOfTargetLeafBoxes >=0) {
		tmp = mxCreateDoubleMatrix(1, statistics->noOfTargetLevels+1, mxREAL);
		p = mxGetPr(tmp);
		for (i=0; i<=statistics->noOfTargetLevels; i++) 
			p[i] = statistics->noOfTargetLeafBoxesInLevel[i];
		mxSetField(stat, 0, "noOfTargetLeafBoxesInLevel", tmp);
	}	
	if (statistics->noOfSourceLeafBoxes >=0) {
		tmp = mxCreateDoubleMatrix(1, statistics->noOfSourceLevels+1, mxREAL);
		p = mxGetPr(tmp);
		for (i=0; i<=statistics->noOfSourceLevels; i++) 
			p[i] = statistics->averageNoOfParticlesPerLeafBoxInLevel[i];
		mxSetField(stat, 0, "averageNoOfParticlesPerLeafBoxInLevel", tmp);	
	}	
	if (statistics->noOfTargetLeafBoxes >=0) {
		tmp = mxCreateDoubleMatrix(1, statistics->noOfTargetLevels+1, mxREAL);
		p = mxGetPr(tmp);
		for (i=0; i<=statistics->noOfTargetLevels; i++) 
			p[i] = statistics->averageNoOfParticlesPerLeafBoxInLevel[i];
		mxSetField(stat, 0, "averageNoOfParticlesPerLeafBoxInLevel", tmp);
	}	
	
	


	mexPutVariable("base", "fmmvstatistics", stat);
	//mxSetName(stat, "fmmvstatistics");
	//mexPutArray(stat, "base");
	
}

