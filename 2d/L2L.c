#include"_fmmv.h"

void init_L2L(FmmvHandle *FMMV, int level) 
{
    if ((FMMV->beta!=0) && (level>=0)) {
        int p = FMMV->pL+1;  // Note pL+1 !!
	_FLOAT_ beta = FMMV->beta;
        _FLOAT_ r;
        _FLOAT_ II[4*FMM_P_MAX+4];
        _FLOAT_ R[2*FMM_P_MAX+3];
        _FLOAT_ B[2*FMM_P_MAX+3];
	double *FF = FMMV->FF;
	double *FF_inv = FMMV->FF_INV;
        int l, j;

        r = ldexp(.35355339059327376220/FMMV->scale, -level);  /* sqrt(2)/4 */
        bessel_i_scaled(r*beta, II, 2*p+2, 4*p+4);
        R[0] = 1.0;
        B[0] = 1.0;
        R[1] = r*beta;
        B[1] = beta;
        for (j=2; j<=2*p+2; j++) {
            R[j] = r*R[j-1];
            B[j] = beta*B[j-1];
        }
        
        for (l=1; l<=p; l++) {
            for (j=1; j<=p; j++) {
                FMMV->A [l][j] = ldexp(II[abs(l-j)]*FF[l-1]*FF_inv[j-1]*R[abs(l-j)+1]*B[j-l+abs(l-j)], level*(l-j)-j);  
                FMMV->Ac[l][j] = ldexp(II[l+j-2]*FF[l-1]*FF_inv[j-1]*R[l+j+1-2]*B[2*j-2], level*(l-j)-j);  
            }
         }
    }   
 
}

void finish_L2L(FmmvHandle *FMMV) 
{ 
	int p;

	if (FMMV->pL+1>FMMV->pM) {
   	    p = FMMV->pL+1;
	}
	else {
	    p = FMMV->pM;
	}
        if (FMMV->beta==0) {	
   	    FMMV_FREE(FMMV, FMMV->A[0], ((p+2)*(p+3)/2+p+1)*sizeof(_FLOAT_));
        }    
        else {
   	    FMMV_FREE(FMMV, FMMV->A[0],  ((p+1)*(p+1))*sizeof(_FLOAT_));
   	    FMMV_FREE(FMMV, FMMV->Ac[0], ((p+1)*(p+1))*sizeof(_FLOAT_));
        }
}


void L2L(FmmvHandle *FMMV, Box *box)
{
	int p = FMMV->pL;
	_FLOAT_ *L;
	int  l, j, i, j_mod_8 ; 
	_FLOAT_ h;
	_FLOAT_ s0[2*FMM_P_MAX];
	_FLOAT_ s1[2*FMM_P_MAX];
	_FLOAT_ s2[2*FMM_P_MAX];
	_FLOAT_ s3[2*FMM_P_MAX];
        _FLOAT_ c[8][2*FMM_P_MAX];

	const _FLOAT_ IS2 = 0.70710678118654752440;


	if (!isTarget(box)||!box->L) return;

	L = box->L;

        memset(c[0], 0, 2*(p+1)*sizeof(_FLOAT_));
        memset(c[1], 0, 2*(p+1)*sizeof(_FLOAT_));
        memset(c[2], 0, 2*(p+1)*sizeof(_FLOAT_));
        memset(c[3], 0, 2*(p+1)*sizeof(_FLOAT_));
        memset(c[4], 0, 2*(p+1)*sizeof(_FLOAT_));
        memset(c[5], 0, 2*(p+1)*sizeof(_FLOAT_));
        memset(c[6], 0, 2*(p+1)*sizeof(_FLOAT_));
        memset(c[7], 0, 2*(p+1)*sizeof(_FLOAT_));

	if (FMMV->beta==0) {
	    for (l=0; l<=p; l++) {
		//j_mod_8 = l%8;
		j_mod_8 = l&7;
                for (j=0; j<=l; j++)  {
			h = FMMV->A[l+1][j+1];
                        c[j_mod_8][2*j]   += h*L[2*l];
                        c[j_mod_8][2*j+1] += h*L[2*l+1];
                         
			j_mod_8--;
			if (j_mod_8<0) {
				j_mod_8 = 7;
			}	
		}	
            } 
	}
        else  {
	    for (l=0; l<=p; l++) {
		//j_mod_8 = l%8;
		j_mod_8 = l&7;
                for (j=0; j<=p; j++)  {
			h = FMMV->A[l+1][j+1];
                        c[j_mod_8][2*j]   += h*L[2*l];
                        c[j_mod_8][2*j+1] += h*L[2*l+1];

			j_mod_8--;
			if (j_mod_8<0) {
				j_mod_8 = 7;
			}	
		}	
		//j_mod_8 = 7-(l%8);
		j_mod_8 = 7-(l&7);
		for (j=1; j<=p; j++) {
			h = FMMV->Ac[l+1][j+1];
                        c[j_mod_8][2*j]   += h*L[2*l];
                        c[j_mod_8][2*j+1] -= h*L[2*l+1];  // "-" ... conjugate!!

			j_mod_8--;
			if (j_mod_8<0) {
				j_mod_8 += 8;
			}	
		}
            } 
        

        }

	for (i=0; i<=2*p; i+=2) {
		s0[i]   = 2.0*(c[0][i]-c[4][i]);
		s0[i+1] = 2.0*(c[0][i+1]-c[4][i+1]);
		s1[i]   = 2.0*IS2*(-c[1][i]+c[1][i+1] + c[5][i]-c[5][i+1]);
		s1[i+1] = 2.0*IS2*(-c[1][i]-c[1][i+1] + c[5][i]+c[5][i+1]);
		s2[i]   = 2.0*(-c[2][i+1]+c[6][i+1]);
		s2[i+1] = 2.0*(c[2][i]-c[6][i]);
		s3[i]   = 2.0*IS2*(c[3][i]+c[3][i+1] - c[7][i]-c[7][i+1]);
		s3[i+1] = 2.0*IS2*(-c[3][i]+c[3][i+1] + c[7][i]-c[7][i+1]);
        }

	if (isTarget(box->child[SW])) { /* k=0 */
		if (!box->child[SW]->L) {
		    box->child[SW]->L = (_FLOAT_*) FMMV_MALLOC(FMMV, 2*(p+1)*sizeof(_FLOAT_));
		    memset(box->child[SW]->L, 0, 2*(p+1)*sizeof(_FLOAT_));
		}	
		L = box->child[SW]->L;
		for (i=0; i<=2*p; i+=2) {
			L[i  ] += s0[i  ] + s1[i  ] + s2[i  ] + s3[i];
			L[i+1] += s0[i+1] + s1[i+1] + s2[i+1] + s3[i+1];
		}
	}	

	if (isTarget(box->child[SE])) { /* k=1 */
		if (!box->child[SE]->L) {
		    box->child[SE]->L = (_FLOAT_*) FMMV_MALLOC(FMMV, 2*(p+1)*sizeof(_FLOAT_));
		    memset(box->child[SE]->L, 0, 2*(p+1)*sizeof(_FLOAT_));
		}	
		L = box->child[SE]->L;
		for (i=0; i<=2*p; i+=2) {
			L[i  ] += s0[i  ] - s1[i+1] - s2[i  ] + s3[i+1];
			L[i+1] += s0[i+1] + s1[i  ] - s2[i+1] - s3[i  ];
		}
	}	

	if (isTarget(box->child[NE])) { /* k=2 */
		if (!box->child[NE]->L) {
		    box->child[NE]->L = (_FLOAT_*) FMMV_MALLOC(FMMV, 2*(p+1)*sizeof(_FLOAT_));
		    memset(box->child[NE]->L, 0, 2*(p+1)*sizeof(_FLOAT_));
		}	
		L = box->child[NE]->L;
		for (i=0; i<=2*p; i+=2) {
			L[i  ] += s0[i  ] - s1[i  ] + s2[i  ] - s3[i];
			L[i+1] += s0[i+1] - s1[i+1] + s2[i+1] - s3[i+1];
		}
	}	
	
	if (isTarget(box->child[NW])) { /* k=3 */
		if (!box->child[NW]->L) {
		    box->child[NW]->L = (_FLOAT_*) FMMV_MALLOC(FMMV, 2*(p+1)*sizeof(_FLOAT_));
		    memset(box->child[NW]->L, 0, 2*(p+1)*sizeof(_FLOAT_));
		}	
		L = box->child[NW]->L;
		for (i=0; i<=2*p; i+=2) {
			L[i  ] += s0[i  ] + s1[i+1] - s2[i  ] - s3[i+1];
			L[i+1] += s0[i+1] - s1[i  ] - s2[i+1] + s3[i  ];
		}
	}	
}

