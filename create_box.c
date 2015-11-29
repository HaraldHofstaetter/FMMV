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

#include "_fmmv.h"
#include<stdlib.h>
#include<string.h>  /* memset */
#include<assert.h>
#ifdef USE_PTHREADS
   #include<pthread.h>
#endif

Box* 
createBox(FmmvHandle *FMMV, Box *parent, 
	  _FLOAT_ x, _FLOAT_ y, 
	  #if (FMM_DIM>=3)
	  _FLOAT_ z,
   	  #endif
	  unsigned int first, 
	  unsigned int noOf,
	  unsigned int level, 
	  unsigned int whichChild,
	  enum KindOfBox kindOfBox)
{
    Box * box;
    int err;

   if ((kindOfBox==TARGET) && parent && (parent->child[whichChild])) {
    	    /* box already exists as a SOURCE box */
	    box = parent->child[whichChild];
	    box->firstTarget = first;
	    box->noOfTargets = noOf;
	    /* insert box into list */
	    box->nextTargetBox = FMMV->firstTargetBoxOfLevel[level];
	    FMMV->firstTargetBoxOfLevel[level] = box;
#ifdef USE_PTHREADS
	    err = pthread_mutex_init(&(box->mutex), 0);
	    if (err) return 0; 
#endif    
	    return box;
    }    
    
    box = (Box*) FMMV_MALLOC(FMMV, sizeof(Box));
    if (!box) return 0; 

    box->firstParticle = 0;
    box->noOfParticles = 0;
    box->M = 0;
    box->firstTarget = 0;
    box->noOfTargets = 0;
    box->L = 0;
    
    box->X[XE] = 0;    
    box->X[XN] = 0;
    box->X[XW] = 0;
    box->X[XS] = 0;
    #if (FMM_DIM>=3)
    box->X[XU] = 0;
    box->X[XD] = 0;
    #endif
    box->noOfXin = -1;

    box->X2[XE] = 0;    
    box->X2[XN] = 0;
    box->X2[XW] = 0;
    box->X2[XS] = 0;
    #if (FMM_DIM>=3)
    box->X2[XU] = 0;
    box->X2[XD] = 0;
    #endif
    box->noOfXin2 = -1;
    
    box->x = x;
    box->y = y;
    #if (FMM_DIM>=3)
    box->z = z;
    #endif
    box->parent = parent;
    box->level = level;
    box->whichChild = whichChild;  
    box->nextSourceBox = 0; 
    box->nextTargetBox = 0; 

    box->child[0] = 0;
    box->child[1] = 0;
    box->child[2] = 0;
    box->child[3] = 0;
    #if (FMM_DIM>=3)
    box->child[4] = 0;
    box->child[5] = 0;
    box->child[6] = 0;
    box->child[7] = 0;
    #endif

    box->neighbor[0] = 0;
    box->neighbor[1] = 0;
    box->neighbor[2] = 0;
    box->neighbor[3] = 0;
    box->neighbor[4] = 0;
    box->neighbor[5] = 0;
    box->neighbor[6] = 0;
    box->neighbor[7] = 0;
    #if (FMM_DIM>=3)
    box->neighbor[8] = 0;
    box->neighbor[9] = 0;
    box->neighbor[10] = 0;
    box->neighbor[11] = 0;
    box->neighbor[12] = 0;
    box->neighbor[13] = 0;
    box->neighbor[14] = 0;
    box->neighbor[15] = 0;
    box->neighbor[16] = 0;
    box->neighbor[17] = 0;
    box->neighbor[18] = 0;
    box->neighbor[19] = 0;
    box->neighbor[20] = 0;
    box->neighbor[21] = 0;
    box->neighbor[22] = 0;
    box->neighbor[23] = 0;
    box->neighbor[24] = 0;
    box->neighbor[25] = 0;
    #endif

    if (FMMV->periodicBoundaryConditions) {
    #if (FMM_DIM==2)	    
    if (parent==0) { /* root */
        box->atBoundary = B_S|B_N|B_W|B_E;
    }	
    else switch (whichChild) {
        case SW:
            box->atBoundary = (B_S|B_W) & parent->atBoundary; 	      
        break;	    
	case NW:
            box->atBoundary = (B_N|B_W) & parent->atBoundary; 	      
        break;	    
	case SE:
            box->atBoundary = (B_S|B_E) & parent->atBoundary; 	      
        break;	    
	case NE:
            box->atBoundary = (B_N|B_E) & parent->atBoundary; 	      
        break;	    
    }	
    #elif (FMM_DIM==3)	    
    if (parent==0) { /* root */
        box->atBoundary = B_S|B_N|B_W|B_E|B_D|B_U;
    }	
    else switch (whichChild) {
        case SWD:
            box->atBoundary = (B_S|B_W|B_D) & parent->atBoundary; 	      
        break;	    
	case NWD:
            box->atBoundary = (B_N|B_W|B_D) & parent->atBoundary; 	      
        break;	    
	case SED:
            box->atBoundary = (B_S|B_E|B_D) & parent->atBoundary; 	      
        break;	    
	case NED:
            box->atBoundary = (B_N|B_E|B_D) & parent->atBoundary; 	      
        break;	    
	case SWU:
            box->atBoundary = (B_S|B_W|B_U) & parent->atBoundary; 	      
        break;	    
	case NWU:
            box->atBoundary = (B_N|B_W|B_U) & parent->atBoundary; 	      
        break;	    
	case SEU:
            box->atBoundary = (B_S|B_E|B_U) & parent->atBoundary; 	      
        break;	    
	case NEU:
            box->atBoundary = (B_N|B_E|B_U) & parent->atBoundary; 	      
        break;	    
    }
    #endif
    }

    if (kindOfBox & SOURCE) {
            box->firstParticle = first;
            box->noOfParticles = noOf;
	    /* insert box into list */
	    box->nextSourceBox = FMMV->firstSourceBoxOfLevel[level];
	    FMMV->firstSourceBoxOfLevel[level] = box;
    }	    

    if (kindOfBox & TARGET) {
    	    box->firstTarget = first;
            box->noOfTargets = noOf;
	    box->nextTargetBox = FMMV->firstTargetBoxOfLevel[level];
	    FMMV->firstTargetBoxOfLevel[level] = box;
#ifdef USE_PTHREADS
	    err = pthread_mutex_init(&(box->mutex), 0);
	    if (err) return 0; 
#endif    
    }
    
    return box;
}	


