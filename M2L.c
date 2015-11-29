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

#include "_fmmv.h"
#include "M2L.h"
#include <string.h> /* memset */
#include <stdlib.h> /* malloc, free */

#define FREE_X(FMMV, x) {FMMV_FREE(FMMV, x, FMMV->s_exp*sizeof(_FLOAT_)); x=0;}


static int noOfSourceNeighborsWithSourceChilds(Box *box)
{
        int i;
        int n=0;

        for (i=0; i<26; i++) {
                if (isSource(box->neighbor[i])&& hasSourceChilds(box->neighbor[i])) n++;
        }
        return n;
}

static void gen_M2L_interaction_list(FmmvHandle *FMMV, Box *box, Box **list)
{
        int i,j,k;
        Box *nb;

        memset(list, 0, 8*26*sizeof(Box*));

        k = 0;
        for (i=0; i<26; i++) {
                nb = box->neighbor[i];
                if (isTarget(nb) && hasTargetChilds(nb)) {
                        if (nb->noOfXin == -1) {
                                nb->noOfXin = noOfSourceNeighborsWithSourceChilds(nb);
                                for (j=0; j<8; j++) {
                                        if (isTarget(nb->child[j])) {
                                                list[k] = nb->child[j];
                                                ++(FMMV->noOfStoredXin);
                                                if (FMMV->noOfStoredXin>FMMV->maxNoOfStoredXin) {
                                                        FMMV->maxNoOfStoredXin = FMMV->noOfStoredXin;
                                                }
                                        }
                                        k++;
                                }
                        }
                        else {
                                for (j=0; j<8; j++) {
                                        if (isTarget(nb->child[j])) {
                                                list[k] = nb->child[j];
                                        }
                                        k++;
                                }
                        }
                }
                else {
                        k += 8;
                }
        }
}

static int noOfSourceNeighborsWithSourceChilds2(Box **neighbor2_list, int i0, int i1)
{
        int i;
        int n=0;

        for (i=i0; i<i1; i++) {
                if (isSource(neighbor2_list[i])&& hasSourceChilds(neighbor2_list[i])) n++;
        }
        return n;
}


static void gen_M2L_interaction_list2(FmmvHandle *FMMV, Box *box, Box **list, Box **neighbor2_list)
{
        int i,j,k;
        Box *nb;
        Box *nb_nb2_list[124];
        int non_reduced_neighbors;

	if (FMMV->periodicBoundaryConditions && (box->level==0)) {
		for (i=0; i<124; i++) {
			neighbor2_list[i] = box;
			nb_nb2_list[i] = box;
		}
	}
	else {	
	        gen_neighbor2_list(box, neighbor2_list);
	}	
        non_reduced_neighbors = (FMMV->reducedScheme?26:124);

        memset(list, 0, 8*non_reduced_neighbors*sizeof(Box*));
        k = 0;
        for (i=0; i<non_reduced_neighbors; i++) {
                nb = neighbor2_list[i];
                if (isTarget(nb) && hasTargetChilds(nb)) {
                        if (nb->noOfXin == -1) {
				if (!(FMMV->periodicBoundaryConditions && (box->level==0))) { 
                                	gen_neighbor2_list(nb, nb_nb2_list);
				}	
                                nb->noOfXin = noOfSourceNeighborsWithSourceChilds2(nb_nb2_list, 0, non_reduced_neighbors);
                                for (j=0; j<8; j++) {
                                        if (isTarget(nb->child[j])) {
                                                list[k] = nb->child[j];
                                                ++(FMMV->noOfStoredXin);
                                                if (FMMV->noOfStoredXin>FMMV->maxNoOfStoredXin) {
                                                        FMMV->maxNoOfStoredXin = FMMV->noOfStoredXin;
                                                }
                                        }
                                        k++;
                                }
                        }
                        else {
                                for (j=0; j<8; j++) {
                                        if (isTarget(nb->child[j])) {
                                                list[k] = nb->child[j];
                                        }
                                        k++;
                                }
                        }
                }
                else {
                        k += 8;
                }
        }
}





static void handle_T0(FmmvHandle *FMMV, Box *box, int T0[][9], _FLOAT_ *X, _FLOAT_ *X1, _FLOAT_ *X2)
{
	int i;
	_FLOAT_ *D_X2X = FMMV->D_X2X;
	int s_exp = FMMV->s_exp;
	int s_exp2 = FMMV->s_exp/2;
	
	for (i=0; i<8; i++) {
	        if (isSource(box->child[T0[i][0]])) {
		    if (T0[i][2]) { /* conjugate */
       	        	 VEC_MUL_CJ(s_exp2, D_X2X + s_exp*T0[i][1], X1 + s_exp*T0[i][0], X + s_exp*T0[i][3]);
		    }	 
		    else {
       	        	 VEC_MUL_C(s_exp2, D_X2X + s_exp*T0[i][1], X1 + s_exp*T0[i][0], X + s_exp*T0[i][3]);
		    }	 
                    VEC_ADD(s_exp, X + s_exp*T0[i][3], X + s_exp*T0[i][4], X + s_exp*T0[i][4]);
		    if (T0[i][6]) { /* conjugate */
       	        	 VEC_MUL_CJ(s_exp2, D_X2X + s_exp*T0[i][5], X2 + s_exp*T0[i][0], X + s_exp*T0[i][7]);
		    }	 
		    else {
       	        	 VEC_MUL_C(s_exp2, D_X2X + s_exp*T0[i][5], X2 + s_exp*T0[i][0], X + s_exp*T0[i][7]);
		    }	 
                    VEC_ADD(s_exp, X + s_exp*T0[i][7], X + s_exp*T0[i][8], X + s_exp*T0[i][8]);
		}
        }
        VEC_ADD(s_exp, X + s_exp*XU_D, X + s_exp*XU_all, X + s_exp*XU_all);
        VEC_ADD(s_exp, X + s_exp*XD_D, X + s_exp*XD_all, X + s_exp*XD_all);
}

static void handle_T1(FmmvHandle *FMMV, int T1[][3], int n, _FLOAT_ *X, int index, int dir, Box **IL)
{
	int i;
	_FLOAT_ *D_X2X = FMMV->D_X2X;
	int s_exp = FMMV->s_exp;
	int s_exp2 = FMMV->s_exp/2;
	_FLOAT_ *source = X + s_exp*index;

        SIMD_ALIGN _FLOAT_ z1[FMM_S_EXP_MAX];
        SIMD_ALIGN _FLOAT_ z2[FMM_S_EXP_MAX];

        for (i=0; i<n; i++) {
	       if (IL[T1[i][0]]) {
               	       if (IL[T1[i][1]]) {
				VEC_MUL_CCJ(s_exp2, D_X2X + s_exp*T1[i][2], source, z1, z2);
				IL[T1[i][1]]->X[dir] = VEC_ADD2(FMMV, s_exp, z2, IL[T1[i][1]]->X[dir]);
			}
			else {
				VEC_MUL_C(s_exp2, D_X2X + s_exp*T1[i][2], source, z1);
			}
			IL[T1[i][0]]->X[dir] = VEC_ADD2(FMMV, s_exp, z1, IL[T1[i][0]]->X[dir]);
		}
               	else if (IL[T1[i][1]]) {
			VEC_MUL_CJ(s_exp2, D_X2X + s_exp*T1[i][2], source, z2);
			IL[T1[i][1]]->X[dir] = VEC_ADD2(FMMV, s_exp, z2, IL[T1[i][1]]->X[dir]);
		}
	}	
}	

static void handle_T1_reduced(FmmvHandle *FMMV, int T1[][3], int n, _FLOAT_ *X, int index, int dir, Box **IL)
{
	int i;
	_FLOAT_ *D_X2X = FMMV->D_X2X;
	int s_exp = FMMV->s_exp;
	int s_exp2 = FMMV->s_exp/2;
	_FLOAT_ *source = X + s_exp*index;

        SIMD_ALIGN _FLOAT_ z1[FMM_S_EXP_MAX];
        SIMD_ALIGN _FLOAT_ z2[FMM_S_EXP_MAX];
	
        for (i=0; i<n; i++) {
	       if (IL[T1[i][0]]) {
               	       if (IL[T1[i][1]]) {
				VEC_MUL_CCJ(s_exp2, D_X2X + s_exp*T1[i][2], source, z1, z2);
				IL[T1[i][1]]->X2[dir] = VEC_ADD2(FMMV, s_exp, z2, IL[T1[i][1]]->X2[dir]);
			}
			else {
				VEC_MUL_C(s_exp2, D_X2X + s_exp*T1[i][2], source, z1);
			}
			IL[T1[i][0]]->X2[dir] = VEC_ADD2(FMMV, s_exp, z1, IL[T1[i][0]]->X2[dir]);
		}
               	else if (IL[T1[i][1]]) {
			VEC_MUL_CJ(s_exp2, D_X2X + s_exp*T1[i][2], source, z2);
			IL[T1[i][1]]->X2[dir] = VEC_ADD2(FMMV, s_exp, z2, IL[T1[i][1]]->X2[dir]);
		}
	}	
}	

static void handle_T2(FmmvHandle *FMMV, int T2[][4], int n, _FLOAT_ *X, int dir, Box **IL)
{
	int i;
	_FLOAT_ *D_X2X = FMMV->D_X2X;
	int s_exp = FMMV->s_exp;
	int s_exp2 = FMMV->s_exp/2;

        SIMD_ALIGN _FLOAT_ z1[FMM_S_EXP_MAX];
	
	for (i=0; i<n; i++) {
		if (IL[T2[i][0]]) {
			if (T2[i][2]) { /* conjugate */
       	        		VEC_MUL_CJ(s_exp2, D_X2X + s_exp*T2[i][1], X + s_exp*T2[i][3], z1);
			}
			else {
       	        		VEC_MUL_C(s_exp2, D_X2X + s_exp*T2[i][1], X + s_exp*T2[i][3], z1);
			}	
			IL[T2[i][0]]->X[dir] = VEC_ADD2(FMMV, s_exp, z1, IL[T2[i][0]]->X[dir]);
		}
	}	
}


void M2L(FmmvHandle *FMMV, Box *box)
{
	int s_exp = FMMV->s_exp;
        SIMD_ALIGN _FLOAT_ X1[8*FMM_S_EXP_MAX];
        SIMD_ALIGN _FLOAT_ X2[8*FMM_S_EXP_MAX];
        SIMD_ALIGN _FLOAT_ X[13*FMM_S_EXP_MAX];

        Box *nb;
        int i, j, k;
        Box *IL[208];

        if (!isSource(box) || !hasSourceChilds(box)) return;

       	gen_M2L_interaction_list(FMMV, box, IL);

        /*** U/D lists ***/

        memset(X + s_exp*XU_all, 0, s_exp*sizeof(_FLOAT_));
        memset(X + s_exp*XD_all, 0, s_exp*sizeof(_FLOAT_));
        memset(X + s_exp*XU_D, 0, s_exp*sizeof(_FLOAT_));
        memset(X + s_exp*XD_D, 0, s_exp*sizeof(_FLOAT_));

        M2X(FMMV, 0, box, X1, X2);

	handle_T0(FMMV, box, T0_UD, X, X1, X2);

	handle_T1(FMMV, T1_U, 8, X, XU_D, XU, IL);
	handle_T1(FMMV, T1_U+8, 18, X, XU_all, XU, IL);
	handle_T1(FMMV, T1_D, 8, X, XD_D, XD, IL);
	handle_T1(FMMV, T1_D+8, 18, X, XD_all, XD, IL);

        /*** N/S lists ***/

        memset(X, 0, 12*s_exp*sizeof(_FLOAT_));

        M2X(FMMV, 1, box, X1, X2);

	handle_T0(FMMV, box, T0_NS, X, X1, X2);

	handle_T1(FMMV, T1_N, 4, X, XU_D, XN, IL);
	handle_T1(FMMV, T1_N+4, 12, X, XU_all, XN, IL);
	handle_T1(FMMV, T1_S, 4, X, XD_D, XS, IL);
	handle_T1(FMMV, T1_S+4, 12, X, XD_all, XS, IL);

	handle_T2(FMMV, T2_N, 16, X, XN, IL);
	handle_T2(FMMV, T2_S, 16, X, XS, IL);

        /*** E/W lists ***/

        memset(X, 0, 12*s_exp*sizeof(_FLOAT_));

        M2X(FMMV, 2, box, X1, X2);

	handle_T0(FMMV, box, T0_EW, X, X1, X2);
	
	handle_T1(FMMV, T1_E, 2, X, XU_D, XE, IL);
	handle_T1(FMMV, T1_E+2, 8, X, XU_all, XE, IL);
	handle_T1(FMMV, T1_W, 2, X, XD_D, XW, IL);
	handle_T1(FMMV, T1_W+2, 8, X, XD_all, XW, IL);

	handle_T2(FMMV, T2_E, 20, X, XE, IL);
	handle_T2(FMMV, T2_W, 20, X, XW, IL);


       for (i=0; i<26; i++) {
                nb = box->neighbor[i];

                if (isTarget(nb) && hasTargetChilds(nb)) {
                        --(nb->noOfXin);
                        --(FMMV->noOfStoredXin);
                        if (nb->noOfXin == 0) {
                                X2L(FMMV, 0, nb, 0);
                                X2L(FMMV, 1, nb, 0);
                                X2L(FMMV, 2, nb, 0);
                                for (j=0; j<8; j++) if (nb->child[j]) {
					for (k=0; k<6; k++) {
                                            FREE_X(FMMV, nb->child[j]->X[k]);
                                            nb->child[j]->X[k] = 0;
					}    
                                }
                        }
                }
        }

}

void M2L_ws2(FmmvHandle *FMMV, Box *box)
{
	int s_exp = FMMV->s_exp;
        SIMD_ALIGN _FLOAT_ X1[8*FMM_S_EXP_MAX];
        SIMD_ALIGN _FLOAT_ X2[8*FMM_S_EXP_MAX];
        SIMD_ALIGN _FLOAT_ X[13*FMM_S_EXP_MAX];

        Box *nb;
        int i, j, k;

        Box *IL[992];
        Box *neighbor2_list[124];
		    	
        if (!isSource(box) || !hasSourceChilds(box)) return;

	gen_M2L_interaction_list2(FMMV, box, IL, neighbor2_list);

        /*** U/D lists ***/

        memset(X + s_exp*XU_all, 0, s_exp*sizeof(_FLOAT_));
        memset(X + s_exp*XD_all, 0, s_exp*sizeof(_FLOAT_));
        memset(X + s_exp*XU_D, 0, s_exp*sizeof(_FLOAT_));
        memset(X + s_exp*XD_D, 0, s_exp*sizeof(_FLOAT_));

        M2X(FMMV, 0, box, X1, X2);

	handle_T0(FMMV, box, T0_UD, X, X1, X2);

	handle_T1(FMMV, T1_U_ws2, 18, X, XU_D, XU, IL);
	handle_T1(FMMV, T1_D_ws2, 18, X, XD_D, XD, IL);
	handle_T1(FMMV, T1_U_ws2+18, 100, X, XU_all, XU, IL);
	handle_T1(FMMV, T1_D_ws2+18, 100, X, XD_all, XD, IL);
	
        /*** N/S lists ***/

        memset(X, 0, 12*s_exp*sizeof(_FLOAT_));

        M2X(FMMV, 1, box, X1, X2);

	handle_T0(FMMV, box, T0_NS, X, X1, X2);

	handle_T1(FMMV, T1_N_ws2, 12, X, XU_D, XN, IL);
	handle_T1(FMMV, T1_S_ws2, 12, X, XD_D, XS, IL);
	handle_T1(FMMV, T1_N_ws2+12, 60, X, XU_all, XN, IL);
	handle_T1(FMMV, T1_S_ws2+12, 60, X, XD_all, XS, IL);

	handle_T2(FMMV, T2_N_ws2, 24, X, XN, IL);
	handle_T2(FMMV, T2_S_ws2, 24, X, XS, IL);

        /*** E/W lists ***/

        memset(X, 0, 12*s_exp*sizeof(_FLOAT_));

        M2X(FMMV, 2, box, X1, X2);

	handle_T0(FMMV, box, T0_EW, X, X1, X2);
	
	handle_T1(FMMV, T1_E_ws2, 8, X, XU_D, XE, IL);
	handle_T1(FMMV, T1_W_ws2, 8, X, XD_D, XW, IL);
	handle_T1(FMMV, T1_E_ws2+8, 36, X, XU_all, XE, IL);
	handle_T1(FMMV, T1_W_ws2+8, 36, X, XD_all, XW, IL);

	handle_T2(FMMV, T2_E_ws2, 36, X, XE, IL);
	handle_T2(FMMV, T2_W_ws2, 36, X, XW, IL);


        for (i=0; i< 124; i++) {
		nb = neighbor2_list[i];

                if (isTarget(nb) && hasTargetChilds(nb)) {
                        --(nb->noOfXin);
                        --(FMMV->noOfStoredXin);
                        if (nb->noOfXin == 0) {
                                X2L(FMMV, 0, nb, 0);
                                X2L(FMMV, 1, nb, 0);
                                X2L(FMMV, 2, nb, 0);
                                for (j=0; j<8; j++) if (nb->child[j]) {
					for (k=0; k<6; k++) {
                                            FREE_X(FMMV, nb->child[j]->X[k]);
                                            nb->child[j]->X[k] = 0;
					}    
                                }
                        }
                }
        }

}

void M2L_ws2_reduced(FmmvHandle *FMMV, Box *box0)
{
	int s_exp = FMMV->s_exp;
        SIMD_ALIGN _FLOAT_ X1[8*FMM_S_EXP_MAX];
        SIMD_ALIGN _FLOAT_ X2[8*FMM_S_EXP_MAX];
        SIMD_ALIGN _FLOAT_ X[13*FMM_S_EXP_MAX];
        SIMD_ALIGN _FLOAT_ z1[FMM_S_EXP_MAX];

        Box *nb;
        int i, j, k, kk;
        Box *IL[208]; 
	Box *box;
	Box *ILR[98];
        Box *neighbor2_list[124];
		    	
	if (!((FMMV->periodicBoundaryConditions)&&(box0->level==-1))) {
	        if (!isSource(box0) || !hasSourceChilds(box0)) return;

		for (i=0; i<26; i++) {
			nb = box0->neighbor[i];
			if (isTarget(nb) && hasTargetChilds(nb) && (nb->noOfXin2 == -1)) {
				nb->noOfXin2 = noOfSourceNeighborsWithSourceChilds(nb);
				++(FMMV->noOfStoredXin);
				if (FMMV->noOfStoredXin>FMMV->maxNoOfStoredXin) {
					 FMMV->maxNoOfStoredXin = FMMV->noOfStoredXin;
				}
			}	
		}
	}	
	
	for (k=0; k<8; k++) {
		if ((FMMV->periodicBoundaryConditions)&&(box0->level==-1)&&(k==1)) {	
			break;
		}		
		box = box0->child[k];
	        if (!isSource(box) || !hasSourceChilds(box)) continue;

		gen_M2L_interaction_list2(FMMV, box, IL, neighbor2_list);
			
	 	for (i=0; i<124-26; i++) {
			nb = neighbor2_list[i+26];
			if (isTarget(nb) && hasTargetChilds(nb)) {
				ILR[i] = nb;
			}
			else {
				ILR[i] = 0;
			}
		}
		
        	/*** U/D lists ***/
        	memset(X + s_exp*XU_all, 0, s_exp*sizeof(_FLOAT_));
	        memset(X + s_exp*XD_all, 0, s_exp*sizeof(_FLOAT_));
        	memset(X + s_exp*XU_D, 0, s_exp*sizeof(_FLOAT_));
	        memset(X + s_exp*XD_D, 0, s_exp*sizeof(_FLOAT_));

	        M2X(FMMV, 0, box, X1, X2);

		handle_T0(FMMV, box, T0_UD, X, X1, X2);

		handle_T1(FMMV, T1_U_ws2, 18, X, XU_D, XU, IL);
		handle_T1(FMMV, T1_D_ws2, 18, X, XD_D, XD, IL);

		handle_T1_reduced(FMMV, T1_U_ws2_reduced, 12, X, XU_all, XU, ILR);
		handle_T1_reduced(FMMV, T1_D_ws2_reduced, 12, X, XD_all, XD, ILR);

       		if (ILR[ir_0_0_8]) {
                   	VEC_MUL(s_exp, FMMV->D_X2X + s_exp*dr_0_0_8, X + s_exp*XU_all, z1);
		        ILR[ir_0_0_8]->X2[XU] = VEC_ADD2(FMMV, s_exp, z1, ILR[ir_0_0_8]->X2[XU]);
		}
       		if (ILR[ir_0_0_m8]) {
                    	VEC_MUL(s_exp, FMMV->D_X2X + s_exp*dr_0_0_8, X + s_exp*XD_all, z1);
		        ILR[ir_0_0_m8]->X2[XD] = VEC_ADD2(FMMV, s_exp, z1, ILR[ir_0_0_m8]->X2[XD]);
		}
	
        	/*** N/S lists ***/

        	memset(X, 0, 12*s_exp*sizeof(_FLOAT_));

       		M2X(FMMV, 1, box, X1, X2);

		handle_T0(FMMV, box, T0_NS, X, X1, X2);

		handle_T1(FMMV, T1_N_ws2, 12, X, XU_D, XN, IL);
		handle_T1(FMMV, T1_S_ws2, 12, X, XD_D, XS, IL);

		handle_T1_reduced(FMMV, T1_N_ws2_reduced, 7, X, XU_all, XN, ILR);
		handle_T1_reduced(FMMV, T1_S_ws2_reduced, 7, X, XD_all, XS, ILR);

       		if (ILR[ir_0_8_0]) {
                      	VEC_MUL(s_exp, FMMV->D_X2X + s_exp*dr_0_0_8, X + s_exp*XU_all, z1);
		        ILR[ir_0_8_0]->X2[XN] = VEC_ADD2(FMMV, s_exp, z1, ILR[ir_0_8_0]->X2[XN]);
		}
       		if (ILR[ir_0_m8_0]) {
                	VEC_MUL(s_exp, FMMV->D_X2X + s_exp*dr_0_0_8, X + s_exp*XD_all, z1);
			ILR[ir_0_m8_0]->X2[XS] = VEC_ADD2(FMMV, s_exp, z1, ILR[ir_0_m8_0]->X2[XS]);
		}
		
		handle_T2(FMMV, T2_N_ws2, 24, X, XN, IL);
		handle_T2(FMMV, T2_S_ws2, 24, X, XS, IL);

     		/*** E/W lists ***/

        	memset(X, 0, 12*s_exp*sizeof(_FLOAT_));

	        M2X(FMMV, 2, box, X1, X2);

		handle_T0(FMMV, box, T0_EW, X, X1, X2);
	
		handle_T1(FMMV, T1_E_ws2, 8, X, XU_D, XE, IL);
		handle_T1(FMMV, T1_W_ws2, 8, X, XD_D, XW, IL);

		handle_T1_reduced(FMMV, T1_E_ws2_reduced, 4, X, XU_all, XE, ILR);
		handle_T1_reduced(FMMV, T1_W_ws2_reduced, 4, X, XD_all, XW, ILR);
		
       		if (ILR[ir_8_0_0]) {
                	VEC_MUL(s_exp, FMMV->D_X2X + s_exp*dr_0_0_8, X + s_exp*XU_all, z1);
			ILR[ir_8_0_0]->X2[XE] = VEC_ADD2(FMMV, s_exp, z1, ILR[ir_8_0_0]->X2[XE]);
		}
       		if (ILR[ir_m8_0_0]) {
                	VEC_MUL(s_exp, FMMV->D_X2X + s_exp*dr_0_0_8, X + s_exp*XD_all, z1);
		        ILR[ir_m8_0_0]->X2[XW] = VEC_ADD2(FMMV, s_exp, z1, ILR[ir_m8_0_0]->X2[XW]);
		}
		
		handle_T2(FMMV, T2_E_ws2, 36, X, XE, IL);
		handle_T2(FMMV, T2_W_ws2, 36, X, XW, IL);


		if (!((FMMV->periodicBoundaryConditions)&&(box0->level==-1))) {	
       		    for (i=0; i<26; i++) { //TODO: check 26!
			nb = neighbor2_list[i];

                	if (isTarget(nb) && hasTargetChilds(nb)) {
                        	--(nb->noOfXin);
                        	--(FMMV->noOfStoredXin);
                        	if (nb->noOfXin == 0) {
                                	X2L(FMMV, 0, nb, 0);
                                	X2L(FMMV, 1, nb, 0);
                                	X2L(FMMV, 2, nb, 0);
                                    	for (j=0; j<8; j++) if (nb->child[j]) {
						for (kk=0; kk<6; kk++) {
	                                            FREE_X(FMMV, nb->child[j]->X[kk]);
        	                                    nb->child[j]->X[kk] = 0;
						}    
                                    	}
                        	}
                	}
		    }	
        	}
	}

	if ((FMMV->periodicBoundaryConditions)&&(box0->level==-1)) {	
		box = box0->child[0]; /* box == root */
		X2L(FMMV, 0, box, 0);
		X2L(FMMV, 1, box, 0);
		X2L(FMMV, 2, box, 0);
                for (j=0; j<8; j++) if (box->child[j]) {
			for (k=0; k<6; k++) {
	                       	FREE_X(FMMV, box->child[j]->X[k]);
        	                       box->child[j]->X[k] = 0;
			}    
                }
			
	        for (j=1; j<8; j++) {
		    box0->child[j]= 0;
		}
		
        	X2L(FMMV, 0, box0, 1);
	        X2L(FMMV, 1, box0, 1);
        	X2L(FMMV, 2, box0, 1);

		for (k=0; k<6; k++) {
			FREE_X(FMMV, box->X2[k]);
		   	box->X2[k] = 0;
		}
	}
	else {
		for (i=0; i<26; i++) {
			nb = box0->neighbor[i];

			if (isTarget(nb) && hasTargetChilds(nb)) {
				--(nb->noOfXin2);
				--(FMMV->noOfStoredXin);
				if (nb->noOfXin2 == 0) {
					X2L(FMMV, 0, nb, 1);
					X2L(FMMV, 1, nb, 1);
					X2L(FMMV, 2, nb, 1);
					for (j=0; j<8; j++) if (nb->child[j]) {
						for (k=0; k<6; k++) {
	                                            FREE_X(FMMV, nb->child[j]->X2[k]);
        	                                    nb->child[j]->X2[k] = 0;
						}    
					}
				}
			}
		}
	}
}
