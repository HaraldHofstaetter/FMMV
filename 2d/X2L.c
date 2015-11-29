#include "_fmmv.h"

_FLOAT_ zeros[2*FMM_S_EPS_MAX]; // should be initialized to zero


void X2L(FmmvHandle *FMMV,  Box *box)
{	
	int p  = FMMV->pM;
	int q = FMMV->s_eps;
	_FLOAT_ v;

	_FLOAT_ s0[2*FMM_S_EPS_MAX+4];
	_FLOAT_ s1[2*FMM_S_EPS_MAX+4];
	_FLOAT_ s2[2*FMM_S_EPS_MAX+4];
	_FLOAT_ s3[2*FMM_S_EPS_MAX+4];

	_FLOAT_ *X_N = box->X[XN];
	_FLOAT_ *X_S = box->X[XS];
	_FLOAT_ *X_W = box->X[XW];
	_FLOAT_ *X_E = box->X[XE];

	int n, j;
	
	if (!isTarget(box)) return;

	if (!X_N) X_N = zeros;
	if (!X_S) X_S = zeros;
	if (!X_W) X_W = zeros;
	if (!X_E) X_E = zeros;
	
	_FLOAT_ *L = box->L;

	for (n=0; n<q; n++) {
		s0[2*n  ] =  X_E[2*n  ] + X_N[2*n  ] + X_W[2*n  ] + X_S[2*n ];
		s0[2*n+1] =  X_E[2*n+1] + X_N[2*n+1] + X_W[2*n+1] + X_S[2*n+1];

		s1[2*n  ] = -X_E[2*n  ] - X_N[2*n+1] + X_W[2*n  ] + X_S[2*n+1];
		s1[2*n+1] = -X_E[2*n+1] + X_N[2*n  ] + X_W[2*n+1] - X_S[2*n  ];
		
		s2[2*n  ] =  X_E[2*n  ] - X_N[2*n  ] + X_W[2*n  ] - X_S[2*n  ];
		s2[2*n+1] =  X_E[2*n+1] - X_N[2*n+1] + X_W[2*n+1] - X_S[2*n+1];

		s3[2*n  ] = -X_E[2*n  ] + X_N[2*n+1] + X_W[2*n  ] - X_S[2*n+1];
		s3[2*n+1] = -X_E[2*n+1] - X_N[2*n  ] + X_W[2*n+1] + X_S[2*n  ];
	}

	if (!box->L) {
	    box->L = (_FLOAT_*) FMMV_MALLOC(FMMV, 2*(p+1)*sizeof(_FLOAT_));
	    memset(box->L, 0, 2*(p+1)*sizeof(_FLOAT_));
	}	
	L = box->L;

	
	if (FMMV->s_exp>q) {
		L[0] += X_E[2*q  ] + X_N[2*q  ] + X_W[2*q  ] + X_S[2*q  ];
		L[1] += X_E[2*q+1] + X_N[2*q+1] + X_W[2*q+1] + X_S[2*q+1];
	}
	for (j=0; j<=p; j+=4) {
		for (n=0; n<q; n++) {
			v = FMMV->vander_over_fact[n][j+1];
			L[2*j  ] += v*s0[2*n  ];
			L[2*j+1] += v*s0[2*n+1];
		}	
	}
	for (j=1; j<=p; j+=4) {
		for (n=0; n<q; n++) {
			v = FMMV->vander_over_fact[n][j+1];
			L[2*j  ] += v*s1[2*n  ];
			L[2*j+1] += v*s1[2*n+1];
		}	
	}
	for (j=2; j<=p; j+=4) {
		for (n=0; n<q; n++) {
			v = FMMV->vander_over_fact[n][j+1];
			L[2*j  ] += v*s2[2*n  ];
			L[2*j+1] += v*s2[2*n+1];
		}	
	}
	for (j=3; j<=p; j+=4) {
		for (n=0; n<q; n++) {
			v = FMMV->vander_over_fact[n][j+1];
			L[2*j  ] += v*s3[2*n  ];
			L[2*j+1] += v*s3[2*n+1];
		}	
	}
}	









	

