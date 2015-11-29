/*
 * FMMV - the Fastest Multipole Method of Vienna
 * Copyright (c) 2006 Harald Hofstaetter
 * Institute for Analysis and Scientific Computing
 * Vienna University of Technology
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

//#define USE_PAPI
#undef USE_PAPI

#include"_fmmv.h"
#include<time.h>
#include<stdio.h>
#ifndef _MSC_VER
  #include<unistd.h>
  #include<sys/time.h>
  #include<sys/times.h>
//#else
//#include "timeval.h"
#endif
#include<assert.h>

#ifdef USE_PAPI
#include<papi.h>
#endif


//static double clockticks;

static double get_cpu_time(void) 
{
	//struct tms buf;
	//times(&buf);
	//return (double)buf.tms_utime / clockticks;
	return (double) clock()/CLOCKS_PER_SEC;
}

static double get_elapsed_time(void) 
{
#ifndef _MSC_VER
	struct timeval buf;
	gettimeofday(&buf,0);
	return (double)buf.tv_sec
        + ((double)buf.tv_usec/1000000); 
#else	 
	return (double) clock()/CLOCKS_PER_SEC;
#endif	
}

void stat_init(FmmvHandle *FMMV)
{
	int i,j;
	//clockticks = (double)sysconf(_SC_CLK_TCK);
	for (i=0; i<STAT_MAX; i++) {
		FMMV->statistics.time[i] = -1.0;
		FMMV->statistics.etime[i] = -1.0;
		for( j=0; j<MAX_NUM_PAPI_EVENTS; j++) {
			FMMV->statistics.PAPIvalues[i][j] = -1;
		}	
	}	
}

void stat_start(FmmvHandle *FMMV,enum StatStep step)
{
	FMMV->statistics.time[step] = get_cpu_time();
	FMMV->statistics.etime[step] = get_elapsed_time();
#ifdef USE_PAPI
	if (FMMV->statistics.PAPIeventSet!=PAPI_NULL) {
		if (step==STAT_TOTAL) {
			int retval = PAPI_start(FMMV->statistics.PAPIeventSet);
			if (retval != PAPI_OK) 
				fprintf(stderr, "PAPI error %d: %s\n",retval,PAPI_strerror(retval));
			assert(retval == PAPI_OK);
		}
		else {
			int retval = PAPI_read(FMMV->statistics.PAPIeventSet, FMMV->statistics.PAPIvalues[step]);
			if (retval != PAPI_OK) 
				fprintf(stderr, "PAPI error %d: %s\n",retval,PAPI_strerror(retval));
			assert(retval == PAPI_OK);
		}	
	}
#endif
}
	
void stat_stop(FmmvHandle *FMMV, enum StatStep step)
{
	FMMV->statistics.time[step] = get_cpu_time() - FMMV->statistics.time[step];
	FMMV->statistics.etime[step] = get_elapsed_time() - FMMV->statistics.etime[step];
#ifdef USE_PAPI
    	if (FMMV->statistics.PAPIeventSet!=PAPI_NULL) {
		int retval, i;
		long_long values[MAX_NUM_PAPI_EVENTS];
		
		if (step==STAT_TOTAL) 
			retval = PAPI_stop(FMMV->statistics.PAPIeventSet, values);	
		else {
			retval = PAPI_read(FMMV->statistics.PAPIeventSet, values);	
		}	
		if (retval != PAPI_OK) 
			fprintf(stderr, "PAPI error %d: %s\n",retval,PAPI_strerror(retval));
		assert(retval == PAPI_OK);
		if (step==STAT_TOTAL) {
			for (i=0; i<PAPI_num_events(FMMV->statistics.PAPIeventSet); i++) {
				FMMV->statistics.PAPIvalues[step][i] = values[i];
			}

		}	
		else {			
			for (i=0; i<PAPI_num_events(FMMV->statistics.PAPIeventSet); i++) {
				FMMV->statistics.PAPIvalues[step][i] = values[i]-FMMV->statistics.PAPIvalues[step][i];
			}
		}	
    	}
#endif    
}

void genTreeStatistics(FmmvHandle *FMMV)
{
	int level;
	Box *box;
	
	FMMV->statistics.pM = FMMV->pM;
	FMMV->statistics.pL = FMMV->pL;
	FMMV->statistics.s_eps = FMMV->s_eps;
	FMMV->statistics.s_exp = FMMV->s_exp;
	FMMV->statistics.maxNoOfStoredXin = FMMV->maxNoOfStoredXin;
	FMMV->statistics.maxAllocatedMemory = FMMV->maxAllocatedMemory;

	for (; FMMV->firstSourceBoxOfLevel[FMMV->maxLevel] == 0; FMMV->maxLevel--)
		;
	FMMV->statistics.noOfSourceLevels = FMMV->maxLevel;
	FMMV->statistics.noOfTargetLevels = FMMV->maxLevel;
	
	FMMV->statistics.noOfParticles = 0;
	FMMV->statistics.noOfSourceBoxes = 0;
	FMMV->statistics.noOfSourceLeafBoxes = 0;
	for (level=0; level<=FMMV->maxLevel; level++) {
		FMMV->statistics.noOfParticlesInLevel[level] = 0;
		FMMV->statistics.noOfSourceBoxesInLevel[level] = 0;
		FMMV->statistics.noOfSourceLeafBoxesInLevel[level] = 0;
		for (box=FMMV->firstSourceBoxOfLevel[level]; box!=NULL; box=box->nextSourceBox) {
			++FMMV->statistics.noOfSourceBoxes;
			++FMMV->statistics.noOfSourceBoxesInLevel[level];
			if (isChildless(box)) {
				++FMMV->statistics.noOfSourceLeafBoxes;
				++FMMV->statistics.noOfSourceLeafBoxesInLevel[level];
				FMMV->statistics.noOfParticlesInLevel[level] += box->noOfParticles;
				FMMV->statistics.noOfParticles += box->noOfParticles;
			}
		}
		if (FMMV->statistics.noOfSourceLeafBoxesInLevel[level]>0) {
			FMMV->statistics.averageNoOfParticlesPerLeafBoxInLevel[level] =
				((float) FMMV->statistics.noOfParticlesInLevel[level])
				 / ((float) FMMV->statistics.noOfSourceLeafBoxesInLevel[level]);
		}	
		else {
			FMMV->statistics.averageNoOfParticlesPerLeafBoxInLevel[level] = 0.0;
		}
	}
	for (level=FMMV->maxLevel+1; level<52; level++) {
		FMMV->statistics.noOfParticlesInLevel[level] = 0;
		FMMV->statistics.noOfSourceBoxesInLevel[level] = 0;
		FMMV->statistics.noOfSourceLeafBoxesInLevel[level] = 0;
		FMMV->statistics.averageNoOfParticlesPerLeafBoxInLevel[level] = 0.0;
	}	
		
	FMMV->statistics.averageNoOfParticlesPerLeafBox =
		((float) FMMV->statistics.noOfParticles) / ((float) FMMV->statistics.noOfSourceLeafBoxes); 

	FMMV->statistics.noOfTargetLevels = FMMV->statistics.noOfSourceLevels;
	FMMV->statistics.noOfTargets = FMMV->statistics.noOfParticles;
	FMMV->statistics.noOfTargetBoxes = FMMV->statistics.noOfSourceBoxes;
	FMMV->statistics.noOfTargetLeafBoxes = FMMV->statistics.noOfSourceLeafBoxes;
	FMMV->statistics.averageNoOfTargetsPerLeafBox = FMMV->statistics.averageNoOfParticlesPerLeafBox;
	for (level=0; level<52; level++) {
		FMMV->statistics.noOfTargetsInLevel[level] = FMMV->statistics.noOfParticlesInLevel[level];
		FMMV->statistics.noOfTargetBoxesInLevel[level] = FMMV->statistics.noOfSourceBoxesInLevel[level];
		FMMV->statistics.noOfTargetLeafBoxesInLevel[level] = FMMV->statistics.noOfSourceLeafBoxesInLevel[level];
		FMMV->statistics.averageNoOfTargetsPerLeafBoxInLevel[level] = FMMV->statistics.averageNoOfParticlesPerLeafBoxInLevel[level];
	}	
}

void genTreeStatistics_ST(FmmvHandle *FMMV)
{
	int level;
	Box *box;
	
	FMMV->statistics.pM = FMMV->pM;
	FMMV->statistics.pL = FMMV->pL;
	FMMV->statistics.s_eps = FMMV->s_eps;
	FMMV->statistics.s_exp = FMMV->s_exp;
	FMMV->statistics.maxNoOfStoredXin = FMMV->maxNoOfStoredXin;
	FMMV->statistics.maxAllocatedMemory = FMMV->maxAllocatedMemory;

		;
	FMMV->statistics.noOfSourceLevels = FMMV->maxLevel;
	
	FMMV->statistics.noOfParticles = 0;
	FMMV->statistics.noOfSourceBoxes = 0;
	FMMV->statistics.noOfSourceLeafBoxes = 0;
	for (level=0; level<=FMMV->maxLevel; level++) {
		FMMV->statistics.noOfParticlesInLevel[level] = 0;
		FMMV->statistics.noOfSourceBoxesInLevel[level] = 0;
		FMMV->statistics.noOfSourceLeafBoxesInLevel[level] = 0;
		for (box=FMMV->firstSourceBoxOfLevel[level]; box!=0; box=box->nextSourceBox) {
		    if (isSource(box)) {
			++FMMV->statistics.noOfSourceBoxes;
			++FMMV->statistics.noOfSourceBoxesInLevel[level];
			if (isChildlessSource(box)) {
				++FMMV->statistics.noOfSourceLeafBoxes;
				++FMMV->statistics.noOfSourceLeafBoxesInLevel[level];
				FMMV->statistics.noOfParticlesInLevel[level] += box->noOfParticles;
				FMMV->statistics.noOfParticles += box->noOfParticles;
			}
		     }	
		}
		if (FMMV->statistics.noOfSourceLeafBoxesInLevel[level]>0) {
			FMMV->statistics.averageNoOfParticlesPerLeafBoxInLevel[level] =
				((float) FMMV->statistics.noOfParticlesInLevel[level])
				 / ((float) FMMV->statistics.noOfSourceLeafBoxesInLevel[level]);
		}	
		else {
			FMMV->statistics.averageNoOfParticlesPerLeafBoxInLevel[level] = 0.0;
		}
	}
	for (level=FMMV->maxLevel+1; level<52; level++) {
		FMMV->statistics.noOfParticlesInLevel[level] = 0;
		FMMV->statistics.noOfSourceBoxesInLevel[level] = 0;
		FMMV->statistics.noOfSourceLeafBoxesInLevel[level] = 0;
		FMMV->statistics.averageNoOfParticlesPerLeafBoxInLevel[level] = 0.0;
	}	
		
	FMMV->statistics.averageNoOfParticlesPerLeafBox =
		((float) FMMV->statistics.noOfParticles) / ((float) FMMV->statistics.noOfSourceLeafBoxes); 



		;
	FMMV->statistics.noOfTargetLevels = FMMV->maxTargetLevel;
	
	FMMV->statistics.noOfTargets = 0;
	FMMV->statistics.noOfTargetBoxes = 0;
	FMMV->statistics.noOfTargetLeafBoxes = 0;
	for (level=0; level<=FMMV->maxTargetLevel; level++) {
		FMMV->statistics.noOfTargetsInLevel[level] = 0;
		FMMV->statistics.noOfTargetBoxesInLevel[level] = 0;
		FMMV->statistics.noOfTargetLeafBoxesInLevel[level] = 0;
		for (box=FMMV->firstTargetBoxOfLevel[level]; box!=0; box=box->nextTargetBox) {
		    if (isTarget(box)) {
			++FMMV->statistics.noOfTargetBoxes;
			++FMMV->statistics.noOfTargetBoxesInLevel[level];
			if (isChildlessTarget(box)) {
				++FMMV->statistics.noOfTargetLeafBoxes;
				++FMMV->statistics.noOfTargetLeafBoxesInLevel[level];
				FMMV->statistics.noOfTargetsInLevel[level] += box->noOfTargets;
				FMMV->statistics.noOfTargets += box->noOfTargets;
			}
		    }	
		}
		if (FMMV->statistics.noOfTargetLeafBoxesInLevel[level]>0) {
			FMMV->statistics.averageNoOfTargetsPerLeafBoxInLevel[level] =
				((float) FMMV->statistics.noOfTargetsInLevel[level])
				 / ((float) FMMV->statistics.noOfTargetLeafBoxesInLevel[level]);
		}	
		else {
			FMMV->statistics.averageNoOfTargetsPerLeafBoxInLevel[level] = 0.0;
		}
	}
	for (level=FMMV->maxTargetLevel+1; level<52; level++) {
		FMMV->statistics.noOfTargetsInLevel[level] = 0;
		FMMV->statistics.noOfTargetBoxesInLevel[level] = 0;
		FMMV->statistics.noOfTargetLeafBoxesInLevel[level] = 0;
		FMMV->statistics.averageNoOfTargetsPerLeafBoxInLevel[level] = 0.0;
	}	
		
	FMMV->statistics.averageNoOfTargetsPerLeafBox =
		((float) FMMV->statistics.noOfTargets) / ((float) FMMV->statistics.noOfTargetLeafBoxes); 


}

static char *stat_names[_STAT_LAST_] = {
	"total",
	"buildTree",
	"genM",
	"M2M",
	"M2L",
	"L2L",
	"evalL",
	/* "list1", */
	"evalDirect",
	"list3",
	"list4",
	"list34",
	"farfield",
	"nearfield",
	"initialize",
	"evaluate",
	"finalize",
};	
	


void printFmmvStatistics(struct FmmvStatistics *stat)
{
  int level;
  int i;
  double noOfAllInteractions;
	
  printf("p_M = %i, p_L = %i,  s_eps = %i, S_exp = %i\n", stat->pM, stat->pL, stat->s_eps, stat->s_exp);	

  if (stat->noOfTargets>=0) {
	  printf("\n level    sources  boxes  leafs  sour/leaf    targets  boxes  leafs  targ/leaf\n" );
	  printf("-------------------------------------------------------------------------------\n");
	  for (level=0; level<=(stat->noOfSourceLevels>stat->noOfTargetLevels?stat->noOfSourceLevels:stat->noOfTargetLevels); level++) {
		printf("%5i %10i %6i %6i %10.2f %10i %6i %6i %10.2f\n",
			level,
			stat->noOfParticlesInLevel[level],
			stat->noOfSourceBoxesInLevel[level],
			stat->noOfSourceLeafBoxesInLevel[level],
			stat->averageNoOfParticlesPerLeafBoxInLevel[level],
			stat->noOfTargetsInLevel[level],
			stat->noOfTargetBoxesInLevel[level],
			stat->noOfTargetLeafBoxesInLevel[level],
			stat->averageNoOfTargetsPerLeafBoxInLevel[level]);
	  }		
	  printf("-------------------------------------------------------------------------------\n");
	  printf(" total%10i %6i %6i %10.2f %10i %6i %6i %10.2f\n",
		stat->noOfParticles,
		stat->noOfSourceBoxes,
		stat->noOfSourceLeafBoxes,
		stat->averageNoOfParticlesPerLeafBox,
		stat->noOfTargets,
		stat->noOfTargetBoxes,
		stat->noOfTargetLeafBoxes,
		stat->averageNoOfTargetsPerLeafBox);
  }
  else {
	  printf("\n level    sources  boxes  leafs  sour/leaf\n" );
	  printf("--------------------------------------------\n");
	  for (level=0; level<=(stat->noOfSourceLevels>stat->noOfTargetLevels?stat->noOfSourceLevels:stat->noOfTargetLevels); level++) {
		printf("%5i %10i %6i %6i %10.2f\n",
			level,
			stat->noOfParticlesInLevel[level],
			stat->noOfSourceBoxesInLevel[level],
			stat->noOfSourceLeafBoxesInLevel[level],
			stat->averageNoOfParticlesPerLeafBoxInLevel[level]);
	  }		
	  printf("--------------------------------------------\n");
	  printf(" total%10i %6i %6i %10.2f\n",
		stat->noOfParticles,
		stat->noOfSourceBoxes,
		stat->noOfSourceLeafBoxes,
		stat->averageNoOfParticlesPerLeafBox);
  }
  
  printf("\nmaxNoOfStoredXin = %i\n", stat->maxNoOfStoredXin);
  printf("maxAllocatedMemory = %i\n", stat->maxAllocatedMemory);

  if (stat->noOfTargets>=0) {
  	noOfAllInteractions = ((double) stat->noOfParticles) * stat->noOfTargets;
  }
  else {
  	noOfAllInteractions = ((double) stat->noOfParticles) * (stat->noOfParticles-1);
  }
  printf("noOfDirectInteractions = %lli (%5.2f%% of all interactions)\n", stat->noOfDirectInteractions, ((double) stat->noOfDirectInteractions)/noOfAllInteractions*100.0);
  printf("\nTime:\n"); 
  for (i=1; i<_STAT_LAST_; i++) {
	  if (stat->time[i]>=0.0) {
		  printf(" %10s: %10.2fs %10.2fs\n",  stat_names[i], stat->time[i], stat->etime[i]);
	  }
  }
  /* total time last: */
  printf(" %10s: %10.2fs %10.2fs\n",  stat_names[0], stat->time[0], stat->etime[0]);

#ifdef USE_PAPI
    if (stat->PAPIeventSet!=PAPI_NULL) { 
    	int retval, ev;
	int events[MAX_NUM_PAPI_EVENTS];
	int number = MAX_NUM_PAPI_EVENTS;
	PAPI_event_info_t info;

	retval = PAPI_list_events(stat->PAPIeventSet, events, &number);
	if (retval != PAPI_OK) 
		fprintf(stderr, "PAPI error %d: %s\n",retval,PAPI_strerror(retval));
	assert(retval == PAPI_OK);
	for (ev=0; ev<number; ev++) {
		PAPI_get_event_info(events[ev], &info);
		printf("\n%s (%s):\n", info.symbol, info.long_descr);
	  	for (i=1; i<_STAT_LAST_; i++) {
		  	if (stat->time[i]>=0.0) {
				  printf(" %10s:%20lli\n",  stat_names[i], stat->PAPIvalues[i][ev]);
	       		}
		}	
		printf(" %10s:%20lli\n",  stat_names[0], stat->PAPIvalues[0][ev]);
  	}

    }	
#endif
}

