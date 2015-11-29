/*
 * FMMV - the Fastest Multipole Method of Vienna
 * Copyright (c) 2006 Harald Hofstaetter
 * Institute for Analysis and Scientific Computing
 * Vienna University of Technology
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

#if (FMM_DIM == 2)

static void 
SW2NE(Box *bSW, Box *bNE)
{
    if ((bSW==NULL) || (bNE==NULL)) return;
    
    bSW->neighbor[NE_] = bNE;
    bNE->neighbor[SW_] = bSW;

    SW2NE(bSW->child[NE], bNE->child[SW]);
}	

static void 
NW2SE(Box *bNW, Box *bSE)
{
    if ((bNW==NULL) || (bSE==NULL)) return;
    
    bNW->neighbor[SE_] = bSE;
    bSE->neighbor[NW_] = bNW;

    NW2SE(bNW->child[SE], bSE->child[NW]);
}	


static void
S2N(Box *bS, Box *bN)
{
    if ((bS==NULL) || (bN==NULL)) return;

    bS->neighbor[NM_] = bN;
    bN->neighbor[SM_] = bS;    
    
    SW2NE(bS->child[NW], bN->child[SE]);
    NW2SE(bN->child[SW], bS->child[NE]);    
    
    S2N(bS->child[NW], bN->child[SW]);
    S2N(bS->child[NE], bN->child[SE]);
}

static void
W2E(Box *bW, Box *bE)
{
    if ((bW==NULL) || (bE==NULL)) return;
 
    bW->neighbor[ME_] = bE;
    bE->neighbor[MW_] = bW;
    
    SW2NE(bW->child[SE], bE->child[NW]);
    NW2SE(bW->child[NE], bE->child[SW]);	   
    
    W2E(bW->child[SE], bE->child[SW]);
    W2E(bW->child[NE], bE->child[NW]);
}

void
genNeighborRelations(Box *box)
{
    if (box==NULL) return;    
    
    SW2NE(box->child[SW], box->child[NE]);
    NW2SE(box->child[NW], box->child[SE]);
    
    S2N(box->child[SW], box->child[NW]);
    S2N(box->child[SE], box->child[NE]);
    
    W2E(box->child[SW], box->child[SE]);
    W2E(box->child[NW], box->child[NE]);

    genNeighborRelations(box->child[0]);
    genNeighborRelations(box->child[1]);
    genNeighborRelations(box->child[2]);
    genNeighborRelations(box->child[3]);
}

void
genPeriodicNeighborRelations(Box *root)
{
    int i;
    
    for (i=0; i<FMM_NEIGHBORS_PER_BOX; i++) {
        root->neighbor[i] = root;
    }
    
    SW2NE(root->child[NE], root->child[SW]);
    SW2NE(root->child[NW], root->child[SE]);
    SW2NE(root->child[SE], root->child[NW]);

    NW2SE(root->child[SE], root->child[NW]);
    NW2SE(root->child[NE], root->child[SW]);
    NW2SE(root->child[SW], root->child[NE]);
    
    S2N(root->child[NW], root->child[SW]);
    S2N(root->child[NE], root->child[SE]);
    
    W2E(root->child[SE], root->child[SW]);
    W2E(root->child[NE], root->child[NW]);
}


#elif (FMM_DIM == 3)

static void 
SWD2NEU(Box *bSWD, Box *bNEU)
{
    if ((bSWD==0) || (bNEU==0)) return;
    
    bSWD->neighbor[NEU_] = bNEU;
    bNEU->neighbor[SWD_] = bSWD;

    SWD2NEU(bSWD->child[NEU], bNEU->child[SWD]);
}	

static void 
NWD2SEU(Box *bNWD, Box *bSEU)
{
    if ((bNWD==0) || (bSEU==0)) return;
    
    bNWD->neighbor[SEU_] = bSEU;
    bSEU->neighbor[NWD_] = bNWD;

    NWD2SEU(bNWD->child[SEU], bSEU->child[NWD]);
}	

static void 
SED2NWU(Box *bSED, Box *bNWU)
{
    if ((bSED==0) || (bNWU==0)) return;
    
    bSED->neighbor[NWU_] = bNWU;
    bNWU->neighbor[SED_] = bSED;

    SED2NWU(bSED->child[NWU], bNWU->child[SED]);
}

static void 
NED2SWU(Box *bNED, Box *bSWU)
{
    if ((bNED==0) || (bSWU==0)) return;
    
    bNED->neighbor[SWU_] = bSWU;
    bSWU->neighbor[NED_] = bNED;

    NED2SWU(bNED->child[SWU], bSWU->child[NED]);
}
	
static void
SW2NE(Box *bSW, Box *bNE)
{
    if ((bSW==0) || (bNE==0)) return;

    bSW->neighbor[NEM_] = bNE;
    bNE->neighbor[SWM_] = bSW;

    SWD2NEU(bSW->child[NED], bNE->child[SWU]);
    NED2SWU(bNE->child[SWD], bSW->child[NEU]);    
    
    SW2NE(bSW->child[NED], bNE->child[SWD]);
    SW2NE(bSW->child[NEU], bNE->child[SWU]);        
}	

static void
NW2SE(Box *bNW, Box *bSE)
{
    if ((bNW==0) || (bSE==0)) return;

    bNW->neighbor[SEM_] = bSE;
    bSE->neighbor[NWM_] = bNW;

    NWD2SEU(bNW->child[SED], bSE->child[NWU]);
    SED2NWU(bSE->child[NWD], bNW->child[SEU]);    
    
    NW2SE(bNW->child[SED], bSE->child[NWD]);
    NW2SE(bNW->child[SEU], bSE->child[NWU]);        
}

static void
SD2NU(Box *bSD, Box *bNU)
{
    if ((bSD==0) || (bNU==0)) return;	
    
    bSD->neighbor[NMU_] = bNU;
    bNU->neighbor[SMD_] = bSD;

    SWD2NEU(bSD->child[NWU], bNU->child[SED]);
    SED2NWU(bSD->child[NEU], bNU->child[SWD]);    
    
    SD2NU(bSD->child[NWU], bNU->child[SWD]);    
    SD2NU(bSD->child[NEU], bNU->child[SED]);
}	

static void
ND2SU(Box *bND, Box *bSU)
{
    if ((bND==0) || (bSU==0)) return;	
    
    bND->neighbor[SMU_] = bSU;
    bSU->neighbor[NMD_] = bND;

    NWD2SEU(bND->child[SWU], bSU->child[NED]);
    NED2SWU(bND->child[SEU], bSU->child[NWD]);    
    
    ND2SU(bND->child[SWU], bSU->child[NWD]);    
    ND2SU(bND->child[SEU], bSU->child[NED]);
}	

static void
WD2EU(Box *bWD, Box *bEU)
{
    if ((bWD==0) || (bEU==0)) return;	
    
    bWD->neighbor[MEU_] = bEU;
    bEU->neighbor[MWD_] = bWD;

    SWD2NEU(bWD->child[SEU], bEU->child[NWD]);
    NWD2SEU(bWD->child[NEU], bEU->child[SWD]);    
    
    WD2EU(bWD->child[SEU], bEU->child[SWD]);
    WD2EU(bWD->child[NEU], bEU->child[NWD]);        
}	

static void
ED2WU(Box *bED, Box *bWU)
{
    if ((bED==0) || (bWU==0)) return;	
    
    bED->neighbor[MWU_] = bWU;
    bWU->neighbor[MED_] = bED;

    SED2NWU(bED->child[SWU], bWU->child[NED]);
    NED2SWU(bED->child[NWU], bWU->child[SED]);    
    
    ED2WU(bED->child[SWU], bWU->child[SED]);
    ED2WU(bED->child[NWU], bWU->child[NED]);        
}	


static void
S2N(Box *bS, Box *bN)
{
    if ((bS==0) || (bN==0)) return;

    bS->neighbor[NMM_] = bN;
    bN->neighbor[SMM_] = bS;

    SWD2NEU(bS->child[NWD], bN->child[SEU]);
    NWD2SEU(bN->child[SWD], bS->child[NEU]);    
    SED2NWU(bS->child[NED], bN->child[SWU]);
    NED2SWU(bN->child[SED], bS->child[NWU]);
    
    SW2NE(bS->child[NWD], bN->child[SED]);
    SW2NE(bS->child[NWU], bN->child[SEU]);
    NW2SE(bN->child[SWD], bS->child[NED]);    
    NW2SE(bN->child[SWU], bS->child[NEU]);    
    SD2NU(bS->child[NWD], bN->child[SWU]);    
    SD2NU(bS->child[NED], bN->child[SEU]);    
    ND2SU(bN->child[SWD], bS->child[NWU]);    
    ND2SU(bN->child[SED], bS->child[NEU]);    

    S2N(bS->child[NWD], bN->child[SWD]);
    S2N(bS->child[NED], bN->child[SED]);    
    S2N(bS->child[NWU], bN->child[SWU]);
    S2N(bS->child[NEU], bN->child[SEU]);
}


static void
W2E(Box *bW, Box *bE)
{
    if ((bW==0) || (bE==0)) return;
    
    bW->neighbor[MEM_] = bE;
    bE->neighbor[MWM_] = bW;

    SWD2NEU(bW->child[SED], bE->child[NWU]);
    NWD2SEU(bW->child[NED], bE->child[SWU]);    
    SED2NWU(bE->child[SWD], bW->child[NEU]);
    NED2SWU(bE->child[NWD], bW->child[SEU]);    
    
    SW2NE(bW->child[SED], bE->child[NWD]);
    SW2NE(bW->child[SEU], bE->child[NWU]);
    NW2SE(bW->child[NED], bE->child[SWD]);	   
    NW2SE(bW->child[NEU], bE->child[SWU]);	   
    WD2EU(bW->child[SED], bE->child[SWU]); 
    WD2EU(bW->child[NED], bE->child[NWU]);  
    ED2WU(bE->child[SWD], bW->child[SEU]); 
    ED2WU(bE->child[NWD], bW->child[NEU]);   
    
    W2E(bW->child[SED], bE->child[SWD]);
    W2E(bW->child[NED], bE->child[NWD]);    
    W2E(bW->child[SEU], bE->child[SWU]);
    W2E(bW->child[NEU], bE->child[NWU]);    
}	

static void
D2U(Box *bD, Box *bU)
{
    if ((bD==0) || (bU==0)) return;

    bD->neighbor[MMU_] = bU;    
    bU->neighbor[MMD_] = bD;    

    SWD2NEU(bD->child[SWU], bU->child[NED]);
    NWD2SEU(bD->child[NWU], bU->child[SED]);    
    SED2NWU(bD->child[SEU], bU->child[NWD]);
    NED2SWU(bD->child[NEU], bU->child[SWD]);    
    
    SD2NU(bD->child[SWU], bU->child[NWD]);
    SD2NU(bD->child[SEU], bU->child[NED]);
    ND2SU(bD->child[NWU], bU->child[SWD]);   
    ND2SU(bD->child[NEU], bU->child[SED]);    
    WD2EU(bD->child[SWU], bU->child[SED]);
    WD2EU(bD->child[NWU], bU->child[NED]);
    ED2WU(bD->child[SEU], bU->child[SWD]);   
    ED2WU(bD->child[NEU], bU->child[NWD]);   
    
    D2U(bD->child[SWU], bU->child[SWD]);
    D2U(bD->child[NWU], bU->child[NWD]);
    D2U(bD->child[SEU], bU->child[SED]);
    D2U(bD->child[NEU], bU->child[NED]);
}	

void
genNeighborRelations(Box *box)
{
    if (box==0) return;
    
    SWD2NEU(box->child[SWD], box->child[NEU]);
    NWD2SEU(box->child[NWD], box->child[SEU]);
    SED2NWU(box->child[SED], box->child[NWU]);
    NED2SWU(box->child[NED], box->child[SWU]);

    SW2NE(box->child[SWD], box->child[NED]);
    SW2NE(box->child[SWU], box->child[NEU]);
    NW2SE(box->child[NWD], box->child[SED]);
    NW2SE(box->child[NWU], box->child[SEU]);
    SD2NU(box->child[SWD], box->child[NWU]);
    SD2NU(box->child[SED], box->child[NEU]);
    ND2SU(box->child[NWD], box->child[SWU]);
    ND2SU(box->child[NED], box->child[SEU]);
    WD2EU(box->child[SWD], box->child[SEU]);
    WD2EU(box->child[NWD], box->child[NEU]);    	        
    ED2WU(box->child[SED], box->child[SWU]);    	        
    ED2WU(box->child[NED], box->child[NWU]);    	        
    
    S2N(box->child[SWD], box->child[NWD]);
    S2N(box->child[SED], box->child[NED]);
    S2N(box->child[SWU], box->child[NWU]);
    S2N(box->child[SEU], box->child[NEU]);
    
    W2E(box->child[SWD], box->child[SED]);
    W2E(box->child[NWD], box->child[NED]);
    W2E(box->child[SWU], box->child[SEU]);
    W2E(box->child[NWU], box->child[NEU]);

    D2U(box->child[SWD], box->child[SWU]);
    D2U(box->child[NWD], box->child[NWU]);
    D2U(box->child[SED], box->child[SEU]);
    D2U(box->child[NED], box->child[NEU]);

    genNeighborRelations(box->child[0]);
    genNeighborRelations(box->child[1]);
    genNeighborRelations(box->child[2]);
    genNeighborRelations(box->child[3]);
    genNeighborRelations(box->child[4]);
    genNeighborRelations(box->child[5]);
    genNeighborRelations(box->child[6]);
    genNeighborRelations(box->child[7]);    
}


void
genPeriodicNeighborRelations(Box *root)
{
    int i;
    
    for (i=0; i<FMM_NEIGHBORS_PER_BOX; i++) {
        root->neighbor[i] = root;
    }
    
    SWD2NEU(root->child[NWD], root->child[SEU]);
    SWD2NEU(root->child[SED], root->child[NWU]);
    SWD2NEU(root->child[NED], root->child[SWU]);
    SWD2NEU(root->child[SWU], root->child[NED]);
    SWD2NEU(root->child[NWU], root->child[SED]);
    SWD2NEU(root->child[SEU], root->child[NWD]);    
    SWD2NEU(root->child[NEU], root->child[SWD]);
	    
    NWD2SEU(root->child[SWD], root->child[NEU]);
    NWD2SEU(root->child[SED], root->child[NWU]);
    NWD2SEU(root->child[NED], root->child[SWU]);
    NWD2SEU(root->child[SWU], root->child[NED]);
    NWD2SEU(root->child[NWU], root->child[SED]);
    NWD2SEU(root->child[SEU], root->child[NWD]);
    NWD2SEU(root->child[NEU], root->child[SWD]);
	    
    SED2NWU(root->child[SWD], root->child[NEU]);
    SED2NWU(root->child[NWD], root->child[SEU]);
    SED2NWU(root->child[NED], root->child[SWU]);
    SED2NWU(root->child[SWU], root->child[NED]);
    SED2NWU(root->child[NWU], root->child[SED]);
    SED2NWU(root->child[SEU], root->child[NWD]);
    SED2NWU(root->child[NEU], root->child[SWD]);
	    
    NED2SWU(root->child[SWD], root->child[NEU]);
    NED2SWU(root->child[NWD], root->child[SEU]);
    NED2SWU(root->child[SED], root->child[NWU]);
    NED2SWU(root->child[SWU], root->child[NED]);
    NED2SWU(root->child[NWU], root->child[SED]);
    NED2SWU(root->child[SEU], root->child[NWD]);
    NED2SWU(root->child[NEU], root->child[SWD]);
    
    SW2NE(root->child[NED], root->child[SWD]);
    SW2NE(root->child[NEU], root->child[SWU]);
    SW2NE(root->child[NWD], root->child[SED]);
    SW2NE(root->child[NWU], root->child[SEU]);
    SW2NE(root->child[SED], root->child[NWD]);
    SW2NE(root->child[SEU], root->child[NWU]);

    NW2SE(root->child[SED], root->child[NWD]);
    NW2SE(root->child[SEU], root->child[NWU]);
    NW2SE(root->child[NED], root->child[SWD]);
    NW2SE(root->child[NEU], root->child[SWU]);
    NW2SE(root->child[SWD], root->child[NED]);
    NW2SE(root->child[SWU], root->child[NEU]);

    SD2NU(root->child[NWU], root->child[SWD]);
    SD2NU(root->child[NEU], root->child[SED]);
    SD2NU(root->child[NWD], root->child[SWU]);
    SD2NU(root->child[NED], root->child[SEU]);
    SD2NU(root->child[SWU], root->child[NWD]);
    SD2NU(root->child[SEU], root->child[NED]);

    ND2SU(root->child[SWU], root->child[NWD]);
    ND2SU(root->child[SEU], root->child[NED]);
    ND2SU(root->child[NWU], root->child[SWD]);
    ND2SU(root->child[NEU], root->child[SED]);
    ND2SU(root->child[SWD], root->child[NWU]);
    ND2SU(root->child[SED], root->child[NEU]);

    WD2EU(root->child[SEU], root->child[SWD]);
    WD2EU(root->child[NEU], root->child[NWD]);
    WD2EU(root->child[SWU], root->child[SED]);
    WD2EU(root->child[NWU], root->child[NED]);
    WD2EU(root->child[SED], root->child[SWU]);
    WD2EU(root->child[NED], root->child[NWU]);
	    
    ED2WU(root->child[SWU], root->child[SED]);    	        
    ED2WU(root->child[NWU], root->child[NED]);    	         
    ED2WU(root->child[SWD], root->child[SEU]);    	        
    ED2WU(root->child[NWD], root->child[NEU]);    	         
    ED2WU(root->child[SEU], root->child[SWD]);    	        
    ED2WU(root->child[NEU], root->child[NWD]);    	         
    
    S2N(root->child[NWD], root->child[SWD]);
    S2N(root->child[NED], root->child[SED]);
    S2N(root->child[NWU], root->child[SWU]);
    S2N(root->child[NEU], root->child[SEU]);
    
    W2E(root->child[SED], root->child[SWD]);
    W2E(root->child[NED], root->child[NWD]);
    W2E(root->child[SEU], root->child[SWU]);
    W2E(root->child[NEU], root->child[NWU]);

    D2U(root->child[SWU], root->child[SWD]);
    D2U(root->child[NWU], root->child[NWD]);
    D2U(root->child[SEU], root->child[SED]);
    D2U(root->child[NEU], root->child[NED]);    
}

#endif /* FMM_DIM == 3 */



