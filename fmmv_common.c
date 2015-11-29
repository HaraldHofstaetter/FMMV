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

#include"_fmmv.h"
#include<math.h>
#include<stdlib.h>
#include<assert.h>
#include<stdio.h>

#ifdef USE_PTHREADS	
  #include<pthread.h>
#endif

char errbuf[1024] = "generic error";

struct FmmvOptions 
	fmmvGetDefaultOptions
(void)
{
	struct FmmvOptions options;
	options.pM = 6;
	options.pL = 6;
	options.s = 8;
	options.ws = 1;

	options.splitThreshold = -1; 
	options.splitTargetThreshold = -1; 
	options.levels = -1;
	options.directEvalThreshold = -1; 
	
	options.periodicBoundaryConditions = 0;
	options.extrinsicCorrection = 0;
	options.scale = 1.0;
	options.beta = 0.0;
	options.useHilbertOrder = 0;
#if (FMM_PRECISION==0)
	options.directEvalAccuracy = 1;
#else
	options.directEvalAccuracy = 2;
#endif
	options.useFarfieldNearfieldThreads = 0;
	options.reducedScheme = 0;
	options.PAPIeventSet = -1;
/* TODO: handle addtitional options
	options.x = 0;
	options.w = 0;
	options.M = 0;
*/	
	return options;
}

/* 	
//TODO: no global variables
_FLOAT_ user_defined_x[FMM_S_EPS_MAX];
_FLOAT_ user_defined_w[FMM_S_EPS_MAX];
int user_defined_M[FMM_S_EPS_MAX];
*/
	
/***** INITIALIZE ********************************************************/	
char* fmmv_initialize(FmmvHandle **fh, struct FmmvOptions *options, struct FmmvStatistics *statistics)	
{
	int err;
	struct FmmvOptions _options;
	FmmvHandle *FMMV = *fh;
	
	if (!options) {
		_options = fmmvGetDefaultOptions();
		options = &_options;
	}	
	FMMV->statistics.PAPIeventSet = options->PAPIeventSet;
	stat_init(FMMV);
	stat_start(FMMV, STAT_TOTAL);
	stat_start(FMMV, STAT_INITIALIZE);
	if (!FMMV) goto _err;

	/* TODO: check options */
	FMMV->pM = options->pM;
	FMMV->pL = options->pL;
	FMMV->s_eps = options->s;
	FMMV->ws = options->ws;
	FMMV->splitThreshold = options->splitThreshold;
	FMMV->splitTargetThreshold = options->splitTargetThreshold;
	FMMV->maxLevel = options->levels;
	FMMV->directEvalThreshold = options->directEvalThreshold;
	FMMV->periodicBoundaryConditions = options->periodicBoundaryConditions;
	FMMV->extrinsicCorrection = options->extrinsicCorrection;
	FMMV->scale = (_FLOAT_) options->scale;
	FMMV->useHilbertOrder = options->useHilbertOrder;
	FMMV->directEvalAccuracy = options->directEvalAccuracy;
	FMMV->useFarfieldNearfieldThreads = options->useFarfieldNearfieldThreads;
	FMMV->reducedScheme = options->reducedScheme;

	FMMV->lambda = 0;
/* TODO: handle addtitional options */
/*
	if (options->x) { // TODO: no global variables!!!
		for (i=0; i<FMMV->s_eps; i++) {
			user_defined_x[i] = (_FLOAT_) options->x[i];
			user_defined_w[i] = (_FLOAT_) options->w[i];
			user_defined_M[i] = (_FLOAT_) options->M[i];
		}
		FMMV->lambda = user_defined_x;
		FMMV->w = user_defined_w;
		FMMV->M = user_defined_M;
	}	
*/

	/* TODO: search for best values: */
	if (FMMV->splitThreshold==-1) FMMV->splitThreshold = 200; /* TODO: find precision-, dipole_grad-, etc.- dependent better value */
	if (FMMV->splitTargetThreshold==-1) FMMV->splitTargetThreshold = 200; /* TODO: find precision-, dipole_grad-, etc.- dependent better value */
	if (FMMV->splitThreshold==0) FMMV->splitTargetThreshold = 0; /* non adaptive FMM! */

	#if (FMM_PRECISION==0)
	if (FMMV->maxLevel==-1) FMMV->maxLevel = 22;
	#elif (FMM_PRECISION==1)
	if (FMMV->maxLevel==-1) FMMV->maxLevel = 51;
	#endif
	if (FMMV->directEvalThreshold==-1) FMMV->directEvalThreshold = 200; /* TODO: find precision-, dipole_grad-, etc.- dependent better value */
	
	FMMV->allocatedMemory = 0;
	FMMV->maxAllocatedMemory = 0;

	FMMV->noOfDirectInteractions = 0;

	init_all(FMMV); 

	stat_start(FMMV, STAT_BUILD_TREE); 

	FMMV->perm = (int *) FMMV_MALLOC(FMMV, FMMV->NParticles*sizeof(int));
	if (FMMV->perm==0) goto _err;

	if (FMMV->targets) {
		FMMV->permTargets = (int *) FMMV_MALLOC(FMMV, FMMV->NTargets*sizeof(int));
		if (FMMV->permTargets==0) goto _err;
	
		err = buildTree_ST(FMMV);
		if (err) goto _err;
	}
	else {
		FMMV->permTargets = FMMV->perm;
		err = buildTree(FMMV);
		if (err) goto _err;
	}	
	stat_stop(FMMV, STAT_BUILD_TREE); 

	if (FMMV->targets) {
        	FMMV->maxTargetLevel = FMMV->maxLevel;
        	for (; FMMV->firstSourceBoxOfLevel[FMMV->maxLevel] == 0; FMMV->maxLevel--)
			;
        	for (; FMMV->firstTargetBoxOfLevel[FMMV->maxTargetLevel] == 0; FMMV->maxTargetLevel--)
                	;

		genTreeStatistics_ST(FMMV);
	}	
	else {
        	for (; FMMV->firstSourceBoxOfLevel[FMMV->maxLevel] == 0; FMMV->maxLevel--)
			;

		genTreeStatistics(FMMV);
		FMMV->statistics.noOfTargets = -1;
		FMMV->statistics.noOfTargetLevels = -1;
		FMMV->statistics.noOfTargetBoxes = -1;
		FMMV->statistics.noOfTargetLeafBoxes = -1;
		FMMV->statistics.averageNoOfTargetsPerLeafBox = -1;
	}	

	ida_allocate(FMMV);
	copy_particles(FMMV);

	stat_stop(FMMV, STAT_INITIALIZE);
	if (statistics) {
		*statistics = FMMV->statistics;
		statistics->maxAllocatedMemory = FMMV->maxAllocatedMemory;
		statistics->noOfDirectInteractions = FMMV->noOfDirectInteractions;
		statistics->PAPIeventSet = FMMV->statistics.PAPIeventSet;
	}

	*fh = (void*) FMMV;
	return 0;
_err:
	return errbuf;	
}	

/***** EVALUATE ********************************************************/
char* fmmv_evaluate(FmmvHandle *FMMV, struct FmmvStatistics *statistics)
{
	GenericFmmThreadArg STANDARD_THREAD, FARFIELD_THREAD, NEARFIELD_THREAD;
#ifdef USE_PTHREADS
	pthread_t nearfield_thread;
#endif	
	int err;
	
	stat_start(FMMV, STAT_EVALUATE);

	FMMV->noOfStoredXin = 0;
	FMMV->maxNoOfStoredXin = 0;

	STANDARD_THREAD.fh = FMMV;
	FARFIELD_THREAD.fh = FMMV;
	NEARFIELD_THREAD.fh = FMMV;
	STANDARD_THREAD.thread = 0;
	FARFIELD_THREAD.thread = 1;
	NEARFIELD_THREAD.thread = 2;

	copy_charges(FMMV);
	zero_pot(FMMV);
	
	if (FMMV->useFarfieldNearfieldThreads) {
#ifdef USE_PTHREADS	
	   if (FMMV->useFarfieldNearfieldThreads<0) { /* test mode */
#endif	   
		if (FMMV->periodicBoundaryConditions) {
			if (FMMV->splitThreshold == 0) {
			    if (FMMV->ws == 1) {	
					non_adaptive_fmm_periodic(&FARFIELD_THREAD);
					non_adaptive_fmm_periodic(&NEARFIELD_THREAD);
			    }
			    else { /* FMMV->ws == 2 */
					non_adaptive_fmm_periodic_ws2(&FARFIELD_THREAD);
					non_adaptive_fmm_periodic_ws2(&NEARFIELD_THREAD);
			    }

			}
			else {
				assert(0);
				/* adaptive_fmm_periodic - not yet implemented... */
			}
		} 
		else {
			if (FMMV->splitThreshold == 0) {
			    if (FMMV->ws == 1) {	
					non_adaptive_fmm(&FARFIELD_THREAD);
					non_adaptive_fmm(&NEARFIELD_THREAD);
			    }
			    else { /* FMMV->ws == 2 */
					non_adaptive_fmm_ws2(&FARFIELD_THREAD);
					non_adaptive_fmm_ws2(&NEARFIELD_THREAD);
			    }
			}
			else {
			    if (FMMV->ws == 1) {	
				if (FMMV->targets) {
					adaptive_fmm_ST(&FARFIELD_THREAD);
					adaptive_fmm_ST(&NEARFIELD_THREAD);
				}
				else {
					adaptive_fmm(&FARFIELD_THREAD);
					adaptive_fmm(&NEARFIELD_THREAD);
				}
			    }
			    else { /* FMMV->ws == 2 */
				assert(0); 
				/*adaptive_fmm_ws2 not yet implemented */
			    }
			}
		}
#ifdef USE_PTHREADS	
	   }
	   else {
		if (FMMV->periodicBoundaryConditions) {
			if (FMMV->splitThreshold == 0) {
			    if (FMMV->ws == 1) {	
					/* TODO error handling ... */
					err =  pthread_create(&nearfield_thread, 0, 
							(void* (*)(void *))
							non_adaptive_fmm_periodic, (void*) &NEARFIELD_THREAD);
					non_adaptive_fmm_periodic(&FARFIELD_THREAD);
			    }
			    else { /* FMMV->ws == 2 */
					/* TODO error handling ... */
					err =  pthread_create(&nearfield_thread, 0, 
							(void* (*)(void *))
							non_adaptive_fmm_periodic_ws2, (void*) &NEARFIELD_THREAD);
					non_adaptive_fmm_periodic_ws2(&FARFIELD_THREAD);
			    }
			}
			else {
				assert(0);
				/* adaptive_fmm_periodic - not yet implemented... */
			}
		} 
		else {
			if (FMMV->splitThreshold == 0) {
			    if (FMMV->ws == 1) {	
					err =  pthread_create(&nearfield_thread, 0, 
							(void* (*)(void *))
							non_adaptive_fmm, (void*) &NEARFIELD_THREAD);
					non_adaptive_fmm(&FARFIELD_THREAD);
			    }
			    else { /* FMMV->ws == 2 */
					err =  pthread_create(&nearfield_thread, 0, 
							(void* (*)(void *))
							non_adaptive_fmm_ws2, (void*) &NEARFIELD_THREAD);
					non_adaptive_fmm_ws2(&FARFIELD_THREAD);
			    }
			}
			else {
			    if (FMMV->ws == 1) {	
				if (FMMV->targets) {
					err =  pthread_create(&nearfield_thread, 0, 
							(void* (*)(void *))
							adaptive_fmm_ST, (void*) &NEARFIELD_THREAD);
					adaptive_fmm_ST(&FARFIELD_THREAD);
				}
				else {
					err =  pthread_create(&nearfield_thread, 0, 
							(void* (*)(void *))
							adaptive_fmm, (void*) &NEARFIELD_THREAD);
					adaptive_fmm(&FARFIELD_THREAD);
				}
			    }
			    else { /* FMMV->ws == 2 */
				assert(0); 
				/* adaptive_fmm_ws2 not yet implemented */
			    }
			}
		}
		pthread_join(nearfield_thread, 0);
	    }	
#endif	   
	}
	else {
		if (FMMV->periodicBoundaryConditions) {
			if (FMMV->splitThreshold == 0) {
			    if (FMMV->ws == 1) {	
					non_adaptive_fmm_periodic(&STANDARD_THREAD);
			    }
			    else { /* FMMV->ws == 2 */
					non_adaptive_fmm_periodic_ws2(&STANDARD_THREAD);
			    }
			}
			else {
				assert(0);
				/* adaptive_fmm_periodic - not yet implemented... */
			}
		} 
		else {
			if (FMMV->splitThreshold == 0) {
			    if (FMMV->ws == 1) {	
					non_adaptive_fmm(&STANDARD_THREAD);
			    }
			    else { /* FMMV->ws == 2 */
					non_adaptive_fmm_ws2(&STANDARD_THREAD);
			    }
			}
			else {
			    if (FMMV->ws == 1) {	
				if (FMMV->targets) {
					adaptive_fmm_ST(&STANDARD_THREAD);
				}	
				else {	
					adaptive_fmm(&STANDARD_THREAD);
				}
			    }
			    else { /* FMMV->ws == 2 */
				assert(0); 
				/* adaptive_fmm_ws2 not yet implemented */
			    }
			}
		}
	}
	
	if (FMMV->periodicBoundaryConditions && FMMV->extrinsicCorrection) {
		FMMV->extrinsic_correction(FMMV);
	}

	backcopy_pot(FMMV);

	if (FMMV->targets==0) {
		FMMV->noOfDirectInteractions *= 2;
	}

	stat_stop(FMMV, STAT_EVALUATE);
	if (statistics) {
		*statistics = FMMV->statistics;
		statistics->maxNoOfStoredXin = FMMV->maxNoOfStoredXin;
		statistics->maxAllocatedMemory = FMMV->maxAllocatedMemory;
		statistics->noOfDirectInteractions = FMMV->noOfDirectInteractions;
		statistics->PAPIeventSet = FMMV->statistics.PAPIeventSet;
	}
	return 0;
/*  //TODO: error handling
_err:
	return errbuf;	
*/
	
}	

/***** FINALIZE ********************************************************/	
char* fmmv_finalize(FmmvHandle *FMMV, struct FmmvStatistics *statistics)
{
	ida_free(FMMV);
	
	if (FMMV->targets) { 
	  	FMMV_FREE(FMMV, FMMV->perm, FMMV->NParticles*sizeof(int));
  		FMMV_FREE(FMMV, FMMV->permTargets, FMMV->NTargets*sizeof(int));
		freeTree_ST(FMMV, FMMV->firstSourceBoxOfLevel[0]);
	}
	else {
  		FMMV_FREE(FMMV, FMMV->perm, FMMV->NParticles*sizeof(int));
		freeTree(FMMV, FMMV->firstSourceBoxOfLevel[0]);
	}	

	finish_all(FMMV); 
	
	if (FMMV->allocatedMemory!=0) {
	       printf("*** ALLOCATED MEMORY NOT FREED: %u ***\n", FMMV->allocatedMemory);	
	}       
	stat_stop(FMMV, STAT_TOTAL);
	if (statistics) {
		*statistics = FMMV->statistics;
		statistics->maxNoOfStoredXin = FMMV->maxNoOfStoredXin;
		statistics->maxAllocatedMemory = FMMV->maxAllocatedMemory;
		statistics->noOfDirectInteractions = FMMV->noOfDirectInteractions;
		statistics->PAPIeventSet = FMMV->statistics.PAPIeventSet;
	}	
	
	free(FMMV);
		
	return 0;
/*  //TODO: error handling
_err:
	return errbuf;	
*/
}

