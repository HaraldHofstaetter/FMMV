#include"_fmmv.h"



void init_M2M(FmmvHandle *FMMV, int level) 
{
    _FLOAT_ beta = FMMV->beta;
    
    if (beta==0.0) {
        // for all levels the coeffs are the same => calculate them once (for level<0)
        // coeffs also valid for L2L
        if (level<0) {
            int p, l, j;
            _FLOAT_ r, f;
	    _FLOAT_ F[FMM_P_MAX+3]; /* factorials */
	    _FLOAT_ R[FMM_P_MAX+3]; 
            
            if (FMMV->pL+1>FMMV->pM) {
	        p = FMMV->pL+1;
	    }
	    else {
	        p = FMMV->pM;
	    }

	    r = 2*0.3535533905932737622004222/FMMV->scale; /* sqrt(2)/4 */
	    F[0] = 1.0;
	    R[0] = 1.0;
	    for (l=1; l<=FMM_P_MAX+2; l++) {
		    F[l] = l*F[l-1];
		    R[l] = r*R[l-1];
	    }
    
	    FMMV->A[0] = FMMV_MALLOC(FMMV, ((p+2)*(p+3)/2+p+1)*sizeof(_FLOAT_));
	    
	    FMMV->A[0][0] = 1.0;
	    for (l=1;l<=p;l++) {
		    FMMV->A[0][l] = 1.0/(_FLOAT_) l;
	    }
	    FMMV->A[1] = FMMV->A[0]+p+1;
	    FMMV->A[1][0] = R[1]*0.5;
	    FMMV->A[1][1] = R[0]*0.5;
	    for (l=1; l<=p; l++) {
		    FMMV->A[l+1] = FMMV->A[l]+l+1;
		    f = ldexp(1.0, -(l+1));
		    FMMV->A[l+1][0] = f*R[l+1]/((_FLOAT_) (l+1)); 
		    for (j=0; j<=l; j++) {
			    FMMV->A[l+1][j+1] = f*R[l-j]*F[l]/(F[j]*F[l-j]);
		    }
	    }
        }
    }
    else { // beta!=0
        if (level<0)  {  // init
            int p, l;
        
            if (FMMV->pL+1>FMMV->pM) {
		    p = FMMV->pL+1;
	    }
	    else {
		    p = FMMV->pM;
	    }
	    FMMV->A[0] = FMMV_MALLOC(FMMV, ((p+1)*(p+1))*sizeof(_FLOAT_));
	    FMMV->Ac[0] = FMMV_MALLOC(FMMV, ((p+1)*(p+1))*sizeof(_FLOAT_));
	    for (l=0;l<=p;l++) {
		FMMV->A[l] = FMMV->A[0]+l*(p+1);
		FMMV->Ac[l] = FMMV->Ac[0]+l*(p+1);
            }    

        }
        else {
            int p = FMMV->pM;
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
            FMMV->A[0][0] = II[0]*r*beta;
            for (j=1; j<=p; j++) {
                FMMV->A [0][j] = ldexp(II[j]*FF_inv[j-1]*0.5*R[j+1]*B[2*j], -level*j-j);      
                FMMV->Ac[0][j] = ldexp(II[j]*FF_inv[j-1]*0.5*R[j+1]*B[2*j], -level*j-j);     
            }
            
            for (l=1; l<=p; l++) {
                FMMV->A[l][0] = ldexp(II[l]* FF[l-1]*2.0*R[l+1], level*l);  
                for (j=1; j<=p; j++) {
                    FMMV->A [l][j] = ldexp(II[abs(l-j)]*FF[l-1]*FF_inv[j-1]*R[abs(l-j)+1]*B[j-l+abs(l-j)], level*(l-j)-j);  
                    FMMV->Ac[l][j] = ldexp(II[l+j]*FF[l-1]*FF_inv[j-1]*R[l+j+1]*B[2*j], level*(l-j)-j);  
                }
             }
        }
    }    
}

void finish_M2M(FmmvHandle *FMMV) { /* dummy */ }

void M2M(FmmvHandle *FMMV, Box *box)
{
	int p = FMMV->pM;
	_FLOAT_ *M;
	int i, l, j, j_mod_8;
	_FLOAT_ s_re, s_im, h;
	_FLOAT_ s0[2*FMM_P_MAX];
	_FLOAT_ s1[2*FMM_P_MAX];
	_FLOAT_ s2[2*FMM_P_MAX];
	_FLOAT_ s3[2*FMM_P_MAX];
	_FLOAT_ c[8][2*FMM_P_MAX];
	const _FLOAT_ IS2 = 0.70710678118654752440;
	
	if (!isSource(box)) return;

	memset(s0, 0, 2*(p+1)*sizeof(_FLOAT_));
	memset(s1, 0, 2*(p+1)*sizeof(_FLOAT_));
	memset(s2, 0, 2*(p+1)*sizeof(_FLOAT_));
	memset(s3, 0, 2*(p+1)*sizeof(_FLOAT_));

	if (isSource(box->child[NE])&&box->child[NE]->M) { /* NE: k=0 */
		M = box->child[NE]->M;
		for (i=0; i<=2*p; i+=2) {
			s0[i  ] += M[i  ];
			s0[i+1] += M[i+1];
			s1[i  ] += M[i  ];
			s1[i+1] += M[i+1];
			s2[i  ] += M[i  ];
			s2[i+1] += M[i+1];
			s3[i  ] += M[i  ];
			s3[i+1] += M[i+1];
		}	
	}	

	if (isSource(box->child[NW])&&box->child[NW]->M) { /* NW: k=1 */
		M = box->child[NW]->M;
		for (i=0; i<=2*p; i+=2) {
			s0[i  ] += M[i  ];
			s0[i+1] += M[i+1];
			s1[i  ] -= M[i+1];
			s1[i+1] += M[i  ];
			s2[i  ] -= M[i  ];
			s2[i+1] -= M[i+1];
			s3[i  ] += M[i+1];
			s3[i+1] -= M[i  ];
		}	
	}	
	if (isSource(box->child[SW])&&box->child[SW]->M) { /* SW: k=2 */
		M = box->child[SW]->M;
		for (i=0; i<=2*p; i+=2) {
			s0[i  ] += M[i  ];
			s0[i+1] += M[i+1];
			s1[i  ] -= M[i  ];
			s1[i+1] -= M[i+1];
			s2[i  ] += M[i  ];
			s2[i+1] += M[i+1];
			s3[i  ] -= M[i  ];
			s3[i+1] -= M[i+1];
		}	
	}	
	if (isSource(box->child[SE])&&box->child[SE]->M) { /* SE: k=3 */
		M = box->child[SE]->M;
		for (i=0; i<=2*p; i+=2) {
			s0[i  ] += M[i  ];
			s0[i+1] += M[i+1];
			s1[i  ] += M[i+1];
			s1[i+1] -= M[i  ];
			s2[i  ] -= M[i  ];
			s2[i+1] -= M[i+1];
			s3[i  ] -= M[i+1];
			s3[i+1] += M[i  ];
		}	
	}	

	for (i=0; i<=2*p; i+=2) {
		c[0][i  ] =  s0[i  ];
		c[0][i+1] =  s0[i+1];
		c[1][i  ] =  IS2*(s1[i  ] - s1[i+1]);
		c[1][i+1] =  IS2*(s1[i  ] + s1[i+1]);
		c[2][i  ] = -s2[i+1];
		c[2][i+1] =  s2[i  ];
		c[3][i  ] = -IS2*(s3[i  ] + s3[i+1]);
		c[3][i+1] =  IS2*(s3[i  ] - s3[i+1]);
		c[4][i  ] = -s0[i  ];
		c[4][i+1] = -s0[i+1];
		c[5][i  ] =  IS2*(s1[i+1] - s1[i  ]);
		c[5][i+1] = -IS2*(s1[i  ] + s1[i+1]);
		c[6][i  ] =  s2[i+1];
		c[6][i+1] = -s2[i  ];
		c[7][i  ] =  IS2*(s3[i  ] + s3[i+1]);
		c[7][i+1] =  IS2*(s3[i+1] - s3[i  ]);
	}	

	if (!box->M) {
	    box->M = (_FLOAT_*) FMMV_MALLOC(FMMV, 2*(p+1)*sizeof(_FLOAT_));
	    memset(box->M, 0, 2*(p+1)*sizeof(_FLOAT_));
	}	

	M = box->M;

	if (FMMV->beta==0) {
	    for (l=0; l<=p; l++) {
                s_re = 0.0;
                s_im = 0.0;

		//j_mod_8 = l%8;
		j_mod_8 = l&7;
                for (j=0; j<=l; j++)  {
			h = FMMV->A[l][j];
			s_re += h * c[j_mod_8][2*j];
			s_im += h * c[j_mod_8][2*j+1];

			j_mod_8--;
			if (j_mod_8<0) {
				j_mod_8 = 7;
			}	
		}	
                
		
		M[2*l]   += s_re;
		M[2*l+1] += s_im;
            } 
	}
        else  {
	    for (l=0; l<=p; l++) {
                s_re = 0.0;
                s_im = 0.0;
		//j_mod_8 = l%8;
		j_mod_8 = l&7;
                for (j=0; j<=p; j++)  {
			h = FMMV->A[l][j];
			s_re += h * c[j_mod_8][2*j];
			s_im += h * c[j_mod_8][2*j+1];

			j_mod_8--;
			if (j_mod_8<0) {
				j_mod_8 = 7;
			}	
		}	
		//j_mod_8 = 7-(l%8);
		j_mod_8 = 7-(l&7);
		for (j=1; j<=p; j++) {
			h = FMMV->Ac[l][j];
			s_re += h * c[j_mod_8][2*j];
			s_im -= h * c[j_mod_8][2*j+1];  // "-" ... conjugate!!

			j_mod_8--;
			if (j_mod_8<0) {
				j_mod_8 += 8;
			}	
		}
		M[2*l]   += s_re;
		M[2*l+1] += s_im;
            } 


        }
}	
		

	
