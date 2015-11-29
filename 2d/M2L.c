enum {
D_N_1_1, D_N_m1_1, D_N_1_m1, D_N_m1_m1,
D_S_1_1, D_S_m1_1, D_S_1_m1, D_S_m1_m1,
D_W_1_1, D_W_m1_1, D_W_1_m1, D_W_m1_m1,
D_E_1_1, D_E_m1_1, D_E_1_m1, D_E_m1_m1,
D_N_m1_5, D_N_1_5, D_S_m1_m5, D_S_1_m5,
D_W_m5_m1, D_W_m5_1, D_E_5_m1, D_E_5_1,
D_W_m5_m3, D_S_m3_m5, D_S_m5_m5, D_W_m5_3,
D_N_m3_5, D_W_m5_5, D_N_3_5, D_E_5_3,
D_N_5_5, D_S_3_m5, D_E_5_m3, D_E_5_m5,
D_N_m2_4, D_N_0_4, D_N_2_4, D_N_4_4,
D_E_4_2, D_E_4_0, D_E_4_m2, D_E_4_m4,
D_S_2_m4, D_S_0_m4, D_S_m2_m4, D_S_m4_m4,
D_W_m4_m2, D_W_m4_0, D_W_m4_2, D_W_m4_4,
};
#include "_fmmv.h"
#include <string.h> /* memset */
_FLOAT_ *VEC_addmul_c(FmmvHandle *FMMV, int n, _FLOAT_ *a, _FLOAT_ *out, _FLOAT_ *in);
void gen_M2L_interaction_list(FmmvHandle *FMMV, Box *box, Box **dummy);
void prepare_X2L(FmmvHandle *FMMV, Box *box);


void M2L(FmmvHandle *FMMV, Box *box)
{
	_FLOAT_ *D_X2X = FMMV->D_X2X;
	Box *nc;
	int s_exp = FMMV->s_exp;
	int s_exp2 = 2*FMMV->s_exp;
	Box *IL[32]; /*6^2-2^2*/
	_FLOAT_ outE[2*FMM_S_EPS_MAX+4];
	_FLOAT_ outN[2*FMM_S_EPS_MAX+4];
	_FLOAT_ outW[2*FMM_S_EPS_MAX+4];
	_FLOAT_ outS[2*FMM_S_EPS_MAX+4];
	_FLOAT_ s_E[2*FMM_S_EPS_MAX+4];
	_FLOAT_ s_N[2*FMM_S_EPS_MAX+4];
	_FLOAT_ s_W[2*FMM_S_EPS_MAX+4];
	_FLOAT_ s_S[2*FMM_S_EPS_MAX+4];

	if (!isSource(box) || !hasSourceChilds(box)) return;

	gen_M2L_interaction_list(FMMV, box, IL);

	memset(s_E, 0, s_exp2*sizeof(_FLOAT_));
	memset(s_N, 0, s_exp2*sizeof(_FLOAT_));
	memset(s_W, 0, s_exp2*sizeof(_FLOAT_));
	memset(s_S, 0, s_exp2*sizeof(_FLOAT_));


	if (isSource(box->child[SW])) {
		M2X(FMMV, box->child[SW], outN, outS, outW, outE);
		VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_E_1_1, outE, s_E);
		VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_N_1_1, outN, s_N);
		VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_W_1_1, outW, s_W);
		VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_S_1_1, outS, s_S);
		if (box->neighbor[NW_]) {
			nc = box->neighbor[NW_]->child[SE];
			if (nc) {
				nc->X[XN] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_N_m2_4, outN, nc->X[XN]);
			}
		}
		if (box->neighbor[NM_]) {
			nc = box->neighbor[NM_]->child[SW];
			if (nc) {
				nc->X[XN] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_N_0_4, outN, nc->X[XN]);
			}
			nc = box->neighbor[NM_]->child[SE];
			if (nc) {
				nc->X[XN] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_N_2_4, outN, nc->X[XN]);
			}
		}
		if (box->neighbor[NE_]) {
			nc = box->neighbor[NE_]->child[SW];
			if (nc) {
				nc->X[XN] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_N_4_4, outN, nc->X[XN]);
			}
		}
		if (box->neighbor[ME_]) {
			nc = box->neighbor[ME_]->child[NW];
			if (nc) {
				nc->X[XE] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_E_4_2, outE, nc->X[XE]);
			}
			nc = box->neighbor[ME_]->child[SW];
			if (nc) {
				nc->X[XE] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_E_4_0, outE, nc->X[XE]);
			}
		}
		if (box->neighbor[SE_]) {
			nc = box->neighbor[SE_]->child[NW];
			if (nc) {
				nc->X[XE] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_E_4_m2, outE, nc->X[XE]);
			}
		}
	}

	if (isSource(box->child[NW])) {
		M2X(FMMV, box->child[NW], outN, outS, outW, outE);
		VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_E_1_m1, outE, s_E);
		VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_N_1_m1, outN, s_N);
		VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_W_1_m1, outW, s_W);
		VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_S_1_m1, outS, s_S);
		if (box->neighbor[NE_]) {
			nc = box->neighbor[NE_]->child[SW];
			if (nc) {
				nc->X[XE] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_E_4_2, outE, nc->X[XE]);
			}
		}
		if (box->neighbor[ME_]) {
			nc = box->neighbor[ME_]->child[NW];
			if (nc) {
				nc->X[XE] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_E_4_0, outE, nc->X[XE]);
			}
			nc = box->neighbor[ME_]->child[SW];
			if (nc) {
				nc->X[XE] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_E_4_m2, outE, nc->X[XE]);
			}
		}
		if (box->neighbor[SE_]) {
			nc = box->neighbor[SE_]->child[NW];
			if (nc) {
				nc->X[XE] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_E_4_m4, outE, nc->X[XE]);
			}
		}
		if (box->neighbor[SM_]) {
			nc = box->neighbor[SM_]->child[NE];
			if (nc) {
				nc->X[XS] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_S_2_m4, outS, nc->X[XS]);
			}
			nc = box->neighbor[SM_]->child[NW];
			if (nc) {
				nc->X[XS] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_S_0_m4, outS, nc->X[XS]);
			}
		}
		if (box->neighbor[SW_]) {
			nc = box->neighbor[SW_]->child[NE];
			if (nc) {
				nc->X[XS] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_S_m2_m4, outS, nc->X[XS]);
			}
		}
	}

	if (isSource(box->child[NE])) {
		M2X(FMMV, box->child[NE], outN, outS, outW, outE);
		VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_E_m1_m1, outE, s_E);
		VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_N_m1_m1, outN, s_N);
		VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_W_m1_m1, outW, s_W);
		VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_S_m1_m1, outS, s_S);
		if (box->neighbor[SE_]) {
			nc = box->neighbor[SE_]->child[NW];
			if (nc) {
				nc->X[XS] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_S_2_m4, outS, nc->X[XS]);
			}
		}
		if (box->neighbor[SM_]) {
			nc = box->neighbor[SM_]->child[NE];
			if (nc) {
				nc->X[XS] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_S_0_m4, outS, nc->X[XS]);
			}
			nc = box->neighbor[SM_]->child[NW];
			if (nc) {
				nc->X[XS] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_S_m2_m4, outS, nc->X[XS]);
			}
		}
		if (box->neighbor[SW_]) {
			nc = box->neighbor[SW_]->child[NE];
			if (nc) {
				nc->X[XS] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_S_m4_m4, outS, nc->X[XS]);
			}
		}
		if (box->neighbor[MW_]) {
			nc = box->neighbor[MW_]->child[SE];
			if (nc) {
				nc->X[XW] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_W_m4_m2, outW, nc->X[XW]);
			}
			nc = box->neighbor[MW_]->child[NE];
			if (nc) {
				nc->X[XW] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_W_m4_0, outW, nc->X[XW]);
			}
		}
		if (box->neighbor[NW_]) {
			nc = box->neighbor[NW_]->child[SE];
			if (nc) {
				nc->X[XW] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_W_m4_2, outW, nc->X[XW]);
			}
		}
	}

	if (isSource(box->child[SE])) {
		M2X(FMMV, box->child[SE], outN, outS, outW, outE);
		VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_E_m1_1, outE, s_E);
		VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_N_m1_1, outN, s_N);
		VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_W_m1_1, outW, s_W);
		VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_S_m1_1, outS, s_S);
		if (box->neighbor[SW_]) {
			nc = box->neighbor[SW_]->child[NE];
			if (nc) {
				nc->X[XW] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_W_m4_m2, outW, nc->X[XW]);
			}
		}
		if (box->neighbor[MW_]) {
			nc = box->neighbor[MW_]->child[SE];
			if (nc) {
				nc->X[XW] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_W_m4_0, outW, nc->X[XW]);
			}
			nc = box->neighbor[MW_]->child[NE];
			if (nc) {
				nc->X[XW] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_W_m4_2, outW, nc->X[XW]);
			}
		}
		if (box->neighbor[NW_]) {
			nc = box->neighbor[NW_]->child[SE];
			if (nc) {
				nc->X[XW] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_W_m4_4, outW, nc->X[XW]);
			}
		}
		if (box->neighbor[NM_]) {
			nc = box->neighbor[NM_]->child[SW];
			if (nc) {
				nc->X[XN] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_N_m2_4, outN, nc->X[XN]);
			}
			nc = box->neighbor[NM_]->child[SE];
			if (nc) {
				nc->X[XN] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_N_0_4, outN, nc->X[XN]);
			}
		}
		if (box->neighbor[NE_]) {
			nc = box->neighbor[NE_]->child[SW];
			if (nc) {
				nc->X[XN] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_N_2_4, outN, nc->X[XN]);
			}
		}
	}

	if (box->neighbor[NM_]) {
		nc = box->neighbor[NM_]->child[NW];
		if (nc) {
			nc->X[XN] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_N_m1_5, s_N, nc->X[XN]);
		}
		nc = box->neighbor[NM_]->child[NE];
		if (nc) {
			nc->X[XN] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_N_1_5, s_N, nc->X[XN]);
		}
	}

	if (box->neighbor[SM_]) {
		nc = box->neighbor[SM_]->child[SW];
		if (nc) {
			nc->X[XS] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_S_m1_m5, s_S, nc->X[XS]);
		}
		nc = box->neighbor[SM_]->child[SE];
		if (nc) {
			nc->X[XS] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_S_1_m5, s_S, nc->X[XS]);
		}
	}

	if (box->neighbor[MW_]) {
		nc = box->neighbor[MW_]->child[SW];
		if (nc) {
			nc->X[XW] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_W_m5_m1, s_W, nc->X[XW]);
		}
		nc = box->neighbor[MW_]->child[NW];
		if (nc) {
			nc->X[XW] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_W_m5_1, s_W, nc->X[XW]);
		}
	}

	if (box->neighbor[ME_]) {
		nc = box->neighbor[ME_]->child[SE];
		if (nc) {
			nc->X[XE] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_E_5_m1, s_E, nc->X[XE]);
		}
		nc = box->neighbor[ME_]->child[NE];
		if (nc) {
			nc->X[XE] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_E_5_1, s_E, nc->X[XE]);
		}
	}

	if (box->neighbor[SW_]) {
		nc = box->neighbor[SW_]->child[NW];
		if (nc) {
			nc->X[XW] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_W_m5_m3, s_W, nc->X[XW]);
		}
		nc = box->neighbor[SW_]->child[SE];
		if (nc) {
			nc->X[XS] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_S_m3_m5, s_S, nc->X[XS]);
		}
		nc = box->neighbor[SW_]->child[SW];
		if (nc) {
			nc->X[XS] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_S_m5_m5, s_S, nc->X[XS]);
		}
	}

	if (box->neighbor[NW_]) {
		nc = box->neighbor[NW_]->child[SW];
		if (nc) {
			nc->X[XW] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_W_m5_3, s_W, nc->X[XW]);
		}
		nc = box->neighbor[NW_]->child[NE];
		if (nc) {
			nc->X[XN] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_N_m3_5, s_N, nc->X[XN]);
		}
		nc = box->neighbor[NW_]->child[NW];
		if (nc) {
			nc->X[XW] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_W_m5_5, s_W, nc->X[XW]);
		}
	}

	if (box->neighbor[NE_]) {
		nc = box->neighbor[NE_]->child[NW];
		if (nc) {
			nc->X[XN] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_N_3_5, s_N, nc->X[XN]);
		}
		nc = box->neighbor[NE_]->child[SE];
		if (nc) {
			nc->X[XE] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_E_5_3, s_E, nc->X[XE]);
		}
		nc = box->neighbor[NE_]->child[NE];
		if (nc) {
			nc->X[XN] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_N_5_5, s_N, nc->X[XN]);
		}
	}

	if (box->neighbor[SE_]) {
		nc = box->neighbor[SE_]->child[SW];
		if (nc) {
			nc->X[XS] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_S_3_m5, s_S, nc->X[XS]);
		}
		nc = box->neighbor[SE_]->child[NE];
		if (nc) {
			nc->X[XE] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_E_5_m3, s_E, nc->X[XE]);
		}
		nc = box->neighbor[SE_]->child[SE];
		if (nc) {
			nc->X[XE] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_E_5_m5, s_E, nc->X[XE]);
		}
	}
	prepare_X2L(FMMV, box);
}

/*** dummy routines yet... ***/

#include<assert.h>

void M2L_ws2(FmmvHandle *FMMV, Box *box)
{
	assert(0);
}

void M2L_ws2_reduced(FmmvHandle *FMMV, Box *box)
{
	assert(0);
}

