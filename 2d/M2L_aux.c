#include "_fmmv.h"
/* only preliminary, will be part of M2L.c when M2L.h is implemented ... */

static int noOfSourceNeighborsWithSourceChilds(Box *box)
{
        int i;
        int n=0;

        for (i=0; i<8; i++) {
                if (isSource(box->neighbor[i])&& hasSourceChilds(box->neighbor[i])) n++;
        }
        return n;
}

/* static */ void gen_M2L_interaction_list(FmmvHandle *FMMV, Box *box, Box **dummy)
{
        int i,j;
        Box *nb;

        for (i=0; i<8; i++) {
                nb = box->neighbor[i];
                if (isTarget(nb) && hasTargetChilds(nb)) {
                        if (nb->noOfXin == -1) {
                                nb->noOfXin = noOfSourceNeighborsWithSourceChilds(nb);
                                for (j=0; j<4; j++) {
                                        if (isTarget(nb->child[j])) {
                                                ++(FMMV->noOfStoredXin);
                                                if (FMMV->noOfStoredXin>FMMV->maxNoOfStoredXin) {
                                                        FMMV->maxNoOfStoredXin = FMMV->noOfStoredXin;
                                                }
                                        }
                                }
                        }
                }
        }
}


/*
static void gen_M2L_interaction_list(FmmvHandle *FMMV, Box *box, Box **list)
{
        int i,j,k;
        Box *nb;

        memset(list, 0, 4*8*sizeof(Box*));

        k = 0;
        for (i=0; i<8; i++) {
                nb = box->neighbor[i];
                if (isTarget(nb) && hasTargetChilds(nb)) {
                        if (nb->noOfXin == -1) {
                                nb->noOfXin = noOfSourceNeighborsWithSourceChilds(nb);
                                for (j=0; j<4; j++) {
                                        if (isTarget(nb->child[j])) {
                                                list[k] = nb->child[j];
                                                ++(FMMV->noOfStoredXin);
                                                if (FMMV->noOfStoredXin>FMMV->maxNoOfStoredXin) {
                                                        FMMV->maxNoOfStoredXin = FMMV->noOfStoredXin;
                                                }
                                        }, &
                                        k++;
                                }
                        }
                        else {
                                for (j=0; j<4; j++) {
                                        if (isTarget(nb->child[j])) {
                                                list[k] = nb->child[j];
                                        }
                                        k++;
                                }
                        }
                }
                else {
                        k += 4;
                }
        }
}
*/

#define FREE_X(FMMV, x) {FMMV_FREE(FMMV, x, 2*FMMV->s_exp*sizeof(_FLOAT_)); x=0;}
/* static */ void prepare_X2L(FmmvHandle *FMMV, Box *box)
{
       int i,j,k;
       Box *nb;

       for (i=0; i<8; i++) {
                nb = box->neighbor[i];

                if (isTarget(nb) && hasTargetChilds(nb)) {
                        --(nb->noOfXin);
                        --(FMMV->noOfStoredXin);
                        if (nb->noOfXin == 0) {
                                for (j=0; j<4; j++) if (nb->child[j]) {
                                	X2L(FMMV, nb->child[j]);
					for (k=0; k<4; k++) {
                                            FREE_X(FMMV, nb->child[j]->X[k]);
                                            nb->child[j]->X[k] = 0;
					}    
                                }
                        }
                }
        }
}

