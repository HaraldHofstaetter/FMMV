#include "_fmmv.h"
#include<math.h> // ldexp

// Note: in the 3d version of M2X, all child boxes of the given box are handled
// simultanously, 
//

void M2X(FmmvHandle *FMMV,  Box *box, _FLOAT_ *X_N, _FLOAT_ *X_S, _FLOAT_ *X_W, _FLOAT_ *X_E)	
{
	//_FLOAT_ *x = FMMV->lambda;
	//_FLOAT_ *w = FMMV->w;
	int p = FMMV->pM;
	int q = FMMV->s_eps;
	_FLOAT_ *vander_over_fact_n; // vector x_n^k/k! 
	_FLOAT_ v; // vector x_n^k/k! 
	_FLOAT_ *M;

	_FLOAT_ s0[2*FMM_S_EPS_MAX+4];
	_FLOAT_ s1[2*FMM_S_EPS_MAX+4];
	_FLOAT_ s2[2*FMM_S_EPS_MAX+4];
	_FLOAT_ s3[2*FMM_S_EPS_MAX+4];
	int j, n;


	M = box->M;

	for ( n=0; n<q; n++) {
		vander_over_fact_n = FMMV->vander_over_fact[n];
		s0[2*n  ] = 0.0;
		s0[2*n+1] = 0.0;
		for (j=0; j<=p; j+=4) {
			v = vander_over_fact_n[j];
			s0[2*n  ] += v*M[2*j  ];
			s0[2*n+1] += v*M[2*j+1];
		}	
		s1[2*n  ] = 0.0;
		s1[2*n+1] = 0.0;
		for (j=1; j<=p; j+=4) {
			v = vander_over_fact_n[j];
			s1[2*n  ] += v*M[2*j  ];
			s1[2*n+1] += v*M[2*j+1];
		}	
		s2[2*n  ] = 0.0;
		s2[2*n+1] = 0.0;
		for (j=2; j<=p; j+=4) {
			v = vander_over_fact_n[j];
			s2[2*n  ] += v*M[2*j  ];
			s2[2*n+1] += v*M[2*j+1];
		}	
		s3[2*n  ] = 0.0;
		s3[2*n+1] = 0.0;
		for (j=3; j<=p; j+=4) {
			v = vander_over_fact_n[j];
			s3[2*n  ] += v*M[2*j  ];
			s3[2*n+1] += v*M[2*j+1];
		}	
	}
	
	for(n=0; n<q; n++) { //TODO: X_RRII ...
		// k-1==0
		X_E[2*n  ] = s0[2*n  ] + s1[2*n  ] + s2[2*n  ] + s3[2*n  ];
		X_E[2*n+1] = s0[2*n+1] + s1[2*n+1] + s2[2*n+1] + s3[2*n+1];
		
		// k-1==1
		X_N[2*n  ] = s0[2*n  ] + s1[2*n+1] - s2[2*n  ] - s3[2*n+1];
		X_N[2*n+1] = s0[2*n+1] - s1[2*n  ] - s2[2*n+1] + s3[2*n  ];

		// k-1==2
		X_W[2*n  ] = s0[2*n  ] - s1[2*n  ] + s2[2*n  ] - s3[2*n  ];
		X_W[2*n+1] = s0[2*n+1] - s1[2*n+1] + s2[2*n+1] - s3[2*n+1];
		
		// k-1==3
		X_S[2*n  ] = s0[2*n  ] - s1[2*n+1] - s2[2*n  ] + s3[2*n+1];
		X_S[2*n+1] = s0[2*n+1] + s1[2*n  ] - s2[2*n+1] - s3[2*n  ];
	}

	if (FMMV->s_exp>q) {
		_FLOAT_ M0C = -M[0]*(FMMV->C - box->level*0.69314718055994530942); // log(2)
		_FLOAT_ M1C = 0.0; // -M[1]*FMMV->C[box->level];
		X_E[2*q  ] = M0C; 
		X_E[2*q+1] = M1C; 
		X_N[2*q  ] = M0C;
		X_N[2*q+1] = M1C;
		X_W[2*q  ] = M0C;
		X_W[2*q+1] = M1C;
		X_S[2*q  ] = M0C;
		X_S[2*q+1] = M1C;
	}
}	

	

