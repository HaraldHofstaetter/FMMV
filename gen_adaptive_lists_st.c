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

#define GENLIST1_1(DIRECTION)    \
static void genList1_##DIRECTION(Box *box, Box **list1, int *n) \
{                                \
    if (!isSource(box)) return;       \
    if (!hasSourceChilds(box)) { \
	list1[*n] = box;         \
	(*n)++;                  \
	return;                  \
    }                            \
    genList1_##DIRECTION(box->child[DIRECTION], list1, n); \
}

#if (FMM_DIM==2)
GENLIST1_1(SW)
GENLIST1_1(NW)
GENLIST1_1(SE)
GENLIST1_1(NE)
#elif (FMM_DIM==3)
GENLIST1_1(SWD)
GENLIST1_1(NWD)
GENLIST1_1(SED)
GENLIST1_1(NED)
GENLIST1_1(SWU)
GENLIST1_1(NWU)
GENLIST1_1(SEU)
GENLIST1_1(NEU)
#endif	
	
#define GENLIST1_2(DIRECTION,C1,C2) \
static void genList1_##DIRECTION(Box *box, Box **list1, int *n) \
{                                \
    if (!isSource(box)) return;       \
    if (!hasSourceChilds(box)) { \
	list1[*n] = box;         \
	(*n)++;                  \
	return;                  \
    }                            \
    genList1_##DIRECTION(box->child[C1], list1, n); \
    genList1_##DIRECTION(box->child[C2], list1, n); \
}


#if (FMM_DIM==2)
GENLIST1_2(SM,SW,SE)
GENLIST1_2(NM,NW,NE)
GENLIST1_2(MW,SW,NW)
GENLIST1_2(ME,SE,NE)
#elif (FMM_DIM==3)
GENLIST1_2(SWM,SWD,SWU)
GENLIST1_2(NWM,NWD,NWU)
GENLIST1_2(SEM,SED,SEU)
GENLIST1_2(NEM,NED,NEU)
GENLIST1_2(SMD,SWD,SED)
GENLIST1_2(NMD,NWD,NED)
GENLIST1_2(SMU,SWU,SEU)
GENLIST1_2(NMU,NWU,NEU)
GENLIST1_2(MWD,SWD,NWD)
GENLIST1_2(MED,SED,NED)
GENLIST1_2(MWU,SWU,NWU)
GENLIST1_2(MEU,SEU,NEU)
#endif	

	
#define GENLIST1_3(DIRECTION,C1,C2,C3,C4) \
static void genList1_##DIRECTION(Box *box, Box **list1, int *n) \
{                                \
    if (!isSource(box)) return;       \
    if (!hasSourceChilds(box)) {    \
	list1[*n] = box;         \
	(*n)++;                  \
	return;                  \
    }                            \
    genList1_##DIRECTION(box->child[C1], list1, n); \
    genList1_##DIRECTION(box->child[C2], list1, n); \
    genList1_##DIRECTION(box->child[C3], list1, n); \
    genList1_##DIRECTION(box->child[C4], list1, n); \
}		

#if (FMM_DIM==3)
GENLIST1_3(SMM,SWD,SWU,SED,SEU)
GENLIST1_3(NMM,NWD,NWU,NED,NEU)
GENLIST1_3(MWM,SWD,NWD,SWU,NWU)
GENLIST1_3(MEM,SED,NED,SEU,NEU)
GENLIST1_3(MMD,SWD,NWD,SED,NED)
GENLIST1_3(MMU,SWU,NWU,SEU,NEU)
#endif	

#define PARTLIST1_1(DIRECTION,OPPOSITE_DIRECTION) \
if (isSource(box->neighbor[DIRECTION##_])) { \
    genList1_##OPPOSITE_DIRECTION(box->neighbor[DIRECTION##_], list1, n); \
}	                           \
else {                             \
    c = box;                       \
    nb = 0;                     \
    while (c->whichChild == DIRECTION) { \
        c = c->parent;             \
        if (!c) break;             \
        nb = c->neighbor[DIRECTION##_];     \
        if (isSource(nb)) break;             \
    }                              \
    if (isChildlessSource(nb)) { \
	list1[*n] = nb;            \
	(*n)++;                    \
    }                              \
}		

#define PARTLIST1_2(DIRECTION,OPPOSITE_DIRECTION,C1,C2) \
if (isSource(box->neighbor[DIRECTION##_])) { \
genList1_##OPPOSITE_DIRECTION(box->neighbor[DIRECTION##_], list1, n); \
}	                           \
else {                             \
    c = box;                       \
    nb = 0;                     \
    while ((c->whichChild == C1)||(c->whichChild == C2)) { \
        c = c->parent;             \
        if (!c) break;             \
        nb = c->neighbor[DIRECTION##_];     \
        if (isSource(nb)) break;             \
    }                              \
    if (isChildlessSource(nb)) { \
	list1[*n] = nb;            \
	(*n)++;                    \
    }                              \
}

#define PARTLIST1_3(DIRECTION,OPPOSITE_DIRECTION,C1,C2,C3,C4) \
if (isSource(box->neighbor[DIRECTION##_])) { \
    genList1_##OPPOSITE_DIRECTION(box->neighbor[DIRECTION##_], list1, n); \
}	                           \
else {                             \
    c = box;                       \
    nb = 0;                     \
    while ((c->whichChild == C1)||(c->whichChild == C2) \
         ||(c->whichChild == C3)||(c->whichChild == C4)) { \
        c = c->parent;             \
        if (!c) break;             \
        nb = c->neighbor[DIRECTION##_];     \
        if (isSource(nb)) break;             \
    }                              \
    if (isChildlessSource(nb)) { \
	list1[*n] = nb;            \
	(*n)++;                    \
    }                              \
}


static void AddContainedChildlessSourceBoxes(Box *box, Box **list1, int *n)
{
	if (isChildlessSource(box)) {
		list1[*n] = box; 
		(*n)++;
		return;
	}

	if (isSource(box->child[0])) {
		AddContainedChildlessSourceBoxes(box->child[0], list1, n);
	}		
	if (isSource(box->child[1])) {
		AddContainedChildlessSourceBoxes(box->child[1], list1, n);
	}		
	if (isSource(box->child[2])) {
		AddContainedChildlessSourceBoxes(box->child[2], list1, n);
	}		
	if (isSource(box->child[3])) {
		AddContainedChildlessSourceBoxes(box->child[3], list1, n);
	}	
#if (FMM_DIM==3)	
	if (isSource(box->child[4])) {
		AddContainedChildlessSourceBoxes(box->child[4], list1, n);
	}		
	if (isSource(box->child[5])) {
		AddContainedChildlessSourceBoxes(box->child[5], list1, n);
	}		
	if (isSource(box->child[6])) {
		AddContainedChildlessSourceBoxes(box->child[6], list1, n);
	}		
	if (isSource(box->child[7])) {
		AddContainedChildlessSourceBoxes(box->child[7], list1, n);
	}	
#endif	
}



void genList1_ST(Box *box, Box** list1, int *n)
{
   Box* c;
   Box* nb;

   *n = 0;
   
   if (!isChildlessTarget(box)) return;
#if (FMM_DIM==2)
   PARTLIST1_1(SW,NE)
   PARTLIST1_1(NW,SE)
   PARTLIST1_1(SE,NW)
   PARTLIST1_1(NE,SW)

   PARTLIST1_2(SM,NM,SW,SE)
   PARTLIST1_2(NM,SM,NW,NE)
   PARTLIST1_2(MW,ME,SW,NW)
   PARTLIST1_2(ME,MW,SE,NE)
#elif (FMM_DIM==3)   
   PARTLIST1_1(SWD,NEU)
   PARTLIST1_1(NWD,SEU)
   PARTLIST1_1(SED,NWU)
   PARTLIST1_1(NED,SWU)
   PARTLIST1_1(SWU,NED)
   PARTLIST1_1(NWU,SED)
   PARTLIST1_1(SEU,NWD)
   PARTLIST1_1(NEU,SWD)

   PARTLIST1_2(SWM,NEM,SWD,SWU)
   PARTLIST1_2(NWM,SEM,NWD,NWU)
   PARTLIST1_2(SEM,NWM,SED,SEU)
   PARTLIST1_2(NEM,SWM,NED,NEU)
   PARTLIST1_2(SMD,NMU,SWD,SED)
   PARTLIST1_2(NMD,SMU,NWD,NED)
   PARTLIST1_2(SMU,NMD,SWU,SEU)
   PARTLIST1_2(NMU,SMD,NWU,NEU)
   PARTLIST1_2(MWD,MEU,SWD,NWD)
   PARTLIST1_2(MED,MWU,SED,NED)
   PARTLIST1_2(MWU,MED,SWU,NWU)
   PARTLIST1_2(MEU,MWD,SEU,NEU)   
  
   PARTLIST1_3(SMM,NMM,SWD,SWU,SED,SEU)
   PARTLIST1_3(NMM,SMM,NWD,NWU,NED,NEU)
   PARTLIST1_3(MWM,MEM,SWD,NWD,SWU,NWU)
   PARTLIST1_3(MEM,MWM,SED,NED,SEU,NEU)
   PARTLIST1_3(MMD,MMU,SWD,NWD,SED,NED)
   PARTLIST1_3(MMU,MMD,SWU,NWU,SEU,NEU)   
#endif   
   
   if (isSource(box)) {
	AddContainedChildlessSourceBoxes(box, list1, n);
   } 
   else {
	for (box=box->parent; ((box!=0) && !isSource(box)); box=box->parent)
		;
	if (isChildlessSource(box)) {	
		list1[*n] = box; 
		(*n)++;
	}	
   }	
   
   return;
}   

/*********************************************/
#if (FMM_DIM==2)

#define GENLIST3_1(DIRECTION,C1,C2,C3)    \
static void genList3_##DIRECTION(Box *box, Box **list3, int *n) \
{                                  \
    Box *c;	                   \
    if (!isSource(box)) return;    \
    c = box->child[C1];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
    c = box->child[C2];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
    c = box->child[C3];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
                                   \
    genList3_##DIRECTION(box->child[DIRECTION], list3, n); \
}	

#define GENLIST3_2(DIRECTION,D1,D2,C1,C2)    \
static void genList3_##DIRECTION(Box *box, Box **list3, int *n) \
{                                  \
    Box *c;	                   \
    if (!isSource(box)) return;    \
    c = box->child[C1];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
    c = box->child[C2];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
                                   \
    genList3_##DIRECTION(box->child[D1], list3, n); \
    genList3_##DIRECTION(box->child[D2], list3, n); \
}	

GENLIST3_1(SW,NW,SE,NE)
GENLIST3_1(NW,SW,SE,NE)
GENLIST3_1(SE,SW,NW,NE)
GENLIST3_1(NE,SW,NW,SE)	

GENLIST3_2(SM, SW,SE, NW,NE)
GENLIST3_2(NM, NW,NE, SW,SE)
GENLIST3_2(MW, SW,NW, SE,NE)
GENLIST3_2(ME, SE,NE, SW,NW)


void genList3_ST(Box *box, Box **list3, int *n)
{
   *n = 0;

   if (!isChildlessTarget(box)) return;
   genList3_SW(box->neighbor[NE_], list3, n);
   genList3_NW(box->neighbor[SE_], list3, n);
   genList3_SE(box->neighbor[NW_], list3, n);
   genList3_NE(box->neighbor[SW_], list3, n);   
   
   genList3_SM(box->neighbor[NM_], list3, n);
   genList3_NM(box->neighbor[SM_], list3, n);
   genList3_MW(box->neighbor[ME_], list3, n);
   genList3_ME(box->neighbor[MW_], list3, n);
  
}   	


#elif (FMM_DIM==3)

#define GENLIST3_1(DIRECTION,C1,C2,C3,C4,C5,C6,C7)    \
static void genList3_##DIRECTION(Box *box, Box **list3, int *n) \
{                                  \
    Box *c;	                   \
    if (!isSource(box)) return;    \
    c = box->child[C1];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
    c = box->child[C2];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
    c = box->child[C3];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
    c = box->child[C4];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
    c = box->child[C5];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
    c = box->child[C6];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
    c = box->child[C7];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
                                   \
    genList3_##DIRECTION(box->child[DIRECTION], list3, n); \
}	

#define GENLIST3_2(DIRECTION,D1,D2,C1,C2,C3,C4,C5,C6)    \
static void genList3_##DIRECTION(Box *box, Box **list3, int *n) \
{                                  \
    Box *c;	                   \
    if (!isSource(box)) return;    \
    c = box->child[C1];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
    c = box->child[C2];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
    c = box->child[C3];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
    c = box->child[C4];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
    c = box->child[C5];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
    c = box->child[C6];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
                                   \
    genList3_##DIRECTION(box->child[D1], list3, n); \
    genList3_##DIRECTION(box->child[D2], list3, n); \
}	

#define GENLIST3_3(DIRECTION,D1,D2,D3,D4,C1,C2,C3,C4)    \
static void genList3_##DIRECTION(Box *box, Box **list3, int *n) \
{                                  \
    Box *c;	                   \
    if (!isSource(box)) return;    \
    c = box->child[C1];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
    c = box->child[C2];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
    c = box->child[C3];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
    c = box->child[C4];            \
    if (isSource(c)) {             \
	list3[*n] = c;             \
	(*n)++;                    \
    }	                           \
                                   \
    genList3_##DIRECTION(box->child[D1], list3, n); \
    genList3_##DIRECTION(box->child[D2], list3, n); \
    genList3_##DIRECTION(box->child[D3], list3, n); \
    genList3_##DIRECTION(box->child[D4], list3, n); \
}	

GENLIST3_1(SWD,NWD,SED,NED,SWU,NWU,SEU,NEU)
GENLIST3_1(NWD,SWD,SED,NED,SWU,NWU,SEU,NEU)
GENLIST3_1(SED,SWD,NWD,NED,SWU,NWU,SEU,NEU)
GENLIST3_1(NED,SWD,NWD,SED,SWU,NWU,SEU,NEU)	
GENLIST3_1(SWU,NWU,SEU,NEU,SWD,NWD,SED,NED)
GENLIST3_1(NWU,SWU,SEU,NEU,SWD,NWD,SED,NED)
GENLIST3_1(SEU,SWU,NWU,NEU,SWD,NWD,SED,NED)
GENLIST3_1(NEU,SWU,NWU,SEU,SWD,NWD,SED,NED)	

GENLIST3_2(SWM,SWD,SWU,NWD,SED,NED,NWU,SEU,NEU)
GENLIST3_2(NWM,NWD,NWU,SWD,SED,NED,SWU,SEU,NEU)
GENLIST3_2(SEM,SED,SEU,SWD,NWD,NED,SWU,NWU,NEU)
GENLIST3_2(NEM,NED,NEU,SWD,NWD,SED,SWU,NWU,SEU)
GENLIST3_2(SMD,SWD,SED,NWD,NED,SWU,NWU,SEU,NEU)
GENLIST3_2(NMD,NWD,NED,SWD,SED,SWU,NWU,SEU,NEU)
GENLIST3_2(SMU,SWU,SEU,SWD,NWD,SED,NED,NWU,NEU)
GENLIST3_2(NMU,NWU,NEU,SWD,NWD,SED,NED,SWU,SEU)
GENLIST3_2(MWD,SWD,NWD,SED,NED,SWU,NWU,SEU,NEU)
GENLIST3_2(MED,SED,NED,SWD,NWD,SWU,NWU,SEU,NEU)
GENLIST3_2(MWU,SWU,NWU,SWD,NWD,SED,NED,SEU,NEU)
GENLIST3_2(MEU,SEU,NEU,SWD,NWD,SED,NED,SWU,NWU)

GENLIST3_3(SMM,SWD,SED,SWU,SEU,NWD,NED,NWU,NEU)
GENLIST3_3(NMM,NWD,NED,NWU,NEU,SWD,SED,SWU,SEU)
GENLIST3_3(MWM,SWD,NWD,SWU,NWU,SED,NED,SEU,NEU)
GENLIST3_3(MEM,SED,NED,SEU,NEU,SWD,NWD,SWU,NWU)
GENLIST3_3(MMD,SWD,NWD,SED,NED,SWU,NWU,SEU,NEU)
GENLIST3_3(MMU,NEU,SWU,NWU,SEU,NED,SWD,NWD,SED)


void genList3_ST(Box *box, Box **list3, int *n)
{
   *n = 0;

   if (!isChildlessTarget(box)) return;
   genList3_SWD(box->neighbor[NEU_], list3, n);
   genList3_NWD(box->neighbor[SEU_], list3, n);
   genList3_SED(box->neighbor[NWU_], list3, n);
   genList3_NED(box->neighbor[SWU_], list3, n);   
   genList3_SWU(box->neighbor[NED_], list3, n);
   genList3_NWU(box->neighbor[SED_], list3, n);
   genList3_SEU(box->neighbor[NWD_], list3, n);
   genList3_NEU(box->neighbor[SWD_], list3, n);      
   
   genList3_SWM(box->neighbor[NEM_], list3, n);
   genList3_NWM(box->neighbor[SEM_], list3, n);
   genList3_SEM(box->neighbor[NWM_], list3, n);
   genList3_NEM(box->neighbor[SWM_], list3, n);
   genList3_SMD(box->neighbor[NMU_], list3, n);
   genList3_NMD(box->neighbor[SMU_], list3, n);
   genList3_MWD(box->neighbor[MEU_], list3, n);
   genList3_MED(box->neighbor[MWU_], list3, n);
   genList3_SMU(box->neighbor[NMD_], list3, n);
   genList3_NMU(box->neighbor[SMD_], list3, n);
   genList3_MWU(box->neighbor[MED_], list3, n);
   genList3_MEU(box->neighbor[MWD_], list3, n);
  
   genList3_SMM(box->neighbor[NMM_], list3, n);
   genList3_NMM(box->neighbor[SMM_], list3, n);
   genList3_MWM(box->neighbor[MEM_], list3, n);
   genList3_MEM(box->neighbor[MWM_], list3, n);   
   genList3_MMD(box->neighbor[MMU_], list3, n);
   genList3_MMU(box->neighbor[MMD_], list3, n);
}   

#endif /* FMM_DIM==3 */

/*********************************************/

#define PARTLIST4_1(DIRECTION)     \
   if (box->whichChild != DIRECTION) { \
       c = box;                        \
       nb = 0;                      \
       do {                            \
           c = c->parent;              \
           if (!c) break;              \
           nb = c->neighbor[DIRECTION##_]; \
           if (isSource(nb)) break;              \
       } while (c->whichChild == DIRECTION);  \
       if (isChildlessSource(nb)) {  \
	   list4[*n] = nb;             \
	   (*n)++;                     \
       }                               \
   }      

#define PARTLIST4_2(DIRECTION,C1,C2) \
   if ((box->whichChild != C1)&&(box->whichChild != C2)) { \
       c = box;                        \
       nb = 0;                      \
       do {                            \
           c = c->parent;              \
           if (!c) break;              \
           nb = c->neighbor[DIRECTION##_]; \
           if (isSource(nb)) break;              \
       } while ((c->whichChild == C1)||(c->whichChild == C2));  \
       if (isChildlessSource(nb)) {  \
	   list4[*n] = nb;             \
	   (*n)++;                     \
       }                               \
   }    

#define PARTLIST4_3(DIRECTION,C1,C2,C3,C4) \
   if ((box->whichChild != C1)&&(box->whichChild != C2)    \
    && (box->whichChild != C3)&&(box->whichChild != C4)) { \
       c = box;                        \
       nb = 0;                      \
       do {                            \
           c = c->parent;              \
           if (!c) break;              \
           nb = c->neighbor[DIRECTION##_]; \
           if (isSource(nb)) break;              \
       } while ((c->whichChild == C1)||(c->whichChild == C2)   \
             || (c->whichChild == C3)||(c->whichChild == C4)); \
       if (isChildlessSource(nb)) {  \
	   list4[*n] = nb;             \
	   (*n)++;                     \
       }                               \
   }    


void genList4_ST(Box *box, Box **list4, int *n)
{
   Box* c;
   Box* nb; 

   *n = 0;
   
   if (!isTarget(box)) return;
   if (box->parent==0) return; /* box == root */

#if (FMM_DIM==2)
   PARTLIST4_1(NE)
   PARTLIST4_1(SE)
   PARTLIST4_1(NW)
   PARTLIST4_1(SW)

   PARTLIST4_2(SM,SW,SE)
   PARTLIST4_2(NM,NW,NE)
   PARTLIST4_2(MW,SW,NW)
   PARTLIST4_2(ME,SE,NE)
#elif (FMM_DIM==3)
   PARTLIST4_1(NEU)
   PARTLIST4_1(SEU)
   PARTLIST4_1(NWU)
   PARTLIST4_1(SWU)
   PARTLIST4_1(NED)
   PARTLIST4_1(SED)
   PARTLIST4_1(NWD)
   PARTLIST4_1(SWD)

   PARTLIST4_2(SWM,SWD,SWU)
   PARTLIST4_2(NWM,NWD,NWU)
   PARTLIST4_2(SEM,SED,SEU)
   PARTLIST4_2(NEM,NED,NEU)
   PARTLIST4_2(SMD,SWD,SED)
   PARTLIST4_2(NMD,NWD,NED)
   PARTLIST4_2(SMU,SWU,SEU)
   PARTLIST4_2(NMU,NWU,NEU)
   PARTLIST4_2(MWD,SWD,NWD)
   PARTLIST4_2(MED,SED,NED)
   PARTLIST4_2(MWU,SWU,NWU)
   PARTLIST4_2(MEU,SEU,NEU)   
  
   PARTLIST4_3(SMM,SWD,SWU,SED,SEU)
   PARTLIST4_3(NMM,NWD,NWU,NED,NEU)
   PARTLIST4_3(MWM,SWD,NWD,SWU,NWU)
   PARTLIST4_3(MEM,SED,NED,SEU,NEU)
   PARTLIST4_3(MMD,SWD,NWD,SED,NED)
   PARTLIST4_3(MMU,SWU,NWU,SEU,NEU)      
#endif   
}   
