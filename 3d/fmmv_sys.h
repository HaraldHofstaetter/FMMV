/* This file is automatically generated by gen_ida.py */
/* DO NOT EDIT! */

#ifndef _FMMV_H_
#define _FMMV_H_

#define FMM_DIM  3

#if (FMM_PRECISION==0)
	typedef float _FLOAT_;
	#define FMM_P_MAX 20
	#define FMM_S_EPS_MAX 20
	#define FMM_S_EXP_MAX 1200
	#define FMM_SIMD_VECTOR_LENGTH 4
#elif (FMM_PRECISION==1)
	typedef double _FLOAT_;
	#define FMM_P_MAX 30
	#define FMM_S_EPS_MAX 30
	#define FMM_S_EXP_MAX 2400
	#define FMM_SIMD_VECTOR_LENGTH 4
#endif

#define FMM_STANDARD 0

#if (FMM_PRECISION==0)

typedef struct {
	float x[4];
	float y[4];
	float z[4];
	float q[4];
	float pot[4];
} _DATA_standard_V0_t;

typedef struct {
	_DATA_standard_V0_t *_0; /* x , y , z , q , pot  */
} DATA_standard_t;

#elif (FMM_PRECISION==1)

typedef struct {
	double x[4];
	double y[4];
	double z[4];
	double q[4];
	double pot[4];
} _DATA_standard_V0_t;

typedef struct {
	_DATA_standard_V0_t *_0; /* x , y , z , q , pot  */
} DATA_standard_t;
#endif



#define FMM_DIPOLE 1

#if (FMM_PRECISION==0)

typedef struct {
	float x[4];
	float y[4];
	float z[4];
	float q[4];
	float mx[4];
	float my[4];
	float mz[4];
	float pot[4];
} _DATA_dipole_V0_t;

typedef struct {
	_DATA_dipole_V0_t *_0; /* x , y , z , q , mx , my , mz , pot  */
} DATA_dipole_t;

#elif (FMM_PRECISION==1)

typedef struct {
	double x[4];
	double y[4];
	double z[4];
	double q[4];
	double mx[4];
	double my[4];
	double mz[4];
	double pot[4];
} _DATA_dipole_V0_t;

typedef struct {
	_DATA_dipole_V0_t *_0; /* x , y , z , q , mx , my , mz , pot  */
} DATA_dipole_t;
#endif



#define FMM_GRAD 2

#if (FMM_PRECISION==0)

typedef struct {
	float x[4];
	float y[4];
	float z[4];
	float q[4];
	float pot[4];
	float gradx[4];
	float grady[4];
	float gradz[4];
} _DATA_grad_V0_t;

typedef struct {
	_DATA_grad_V0_t *_0; /* x , y , z , q , pot , gradx , grady , gradz  */
} DATA_grad_t;

#elif (FMM_PRECISION==1)

typedef struct {
	double x[4];
	double y[4];
	double z[4];
	double q[4];
	double pot[4];
	double gradx[4];
	double grady[4];
	double gradz[4];
} _DATA_grad_V0_t;

typedef struct {
	_DATA_grad_V0_t *_0; /* x , y , z , q , pot , gradx , grady , gradz  */
} DATA_grad_t;
#endif



#define FMM_DIPOLE_GRAD 3

#if (FMM_PRECISION==0)

typedef struct {
	float x[4];
	float y[4];
	float z[4];
	float q[4];
	float mx[4];
	float my[4];
	float mz[4];
	float pot[4];
	float gradx[4];
	float grady[4];
	float gradz[4];
} _DATA_dipole_grad_V0_t;

typedef struct {
	_DATA_dipole_grad_V0_t *_0; /* x , y , z , q , mx , my , mz , pot , gradx , grady , gradz  */
} DATA_dipole_grad_t;

#elif (FMM_PRECISION==1)

typedef struct {
	double x[4];
	double y[4];
	double z[4];
	double q[4];
	double mx[4];
	double my[4];
	double mz[4];
	double pot[4];
	double gradx[4];
	double grady[4];
	double gradz[4];
} _DATA_dipole_grad_V0_t;

typedef struct {
	_DATA_dipole_grad_V0_t *_0; /* x , y , z , q , mx , my , mz , pot , gradx , grady , gradz  */
} DATA_dipole_grad_t;
#endif



#define FMM_ST_STANDARD 4

#if (FMM_PRECISION==0)

typedef struct {
	float x[4];
	float y[4];
	float z[4];
	float q[4];
} _DATA_standard_ST_V0_t;

typedef struct {
	float tx[4];
	float ty[4];
	float tz[4];
	float pot[4];
} _DATA_standard_ST_V1_t;

typedef struct {
	_DATA_standard_ST_V0_t *_0; /* x , y , z , q  */
	_DATA_standard_ST_V1_t *_1; /* tx , ty , tz , pot  */
} DATA_standard_ST_t;

#elif (FMM_PRECISION==1)

typedef struct {
	double x[4];
	double y[4];
	double z[4];
	double q[4];
} _DATA_standard_ST_V0_t;

typedef struct {
	double tx[4];
	double ty[4];
	double tz[4];
	double pot[4];
} _DATA_standard_ST_V1_t;

typedef struct {
	_DATA_standard_ST_V0_t *_0; /* x , y , z , q  */
	_DATA_standard_ST_V1_t *_1; /* tx , ty , tz , pot  */
} DATA_standard_ST_t;
#endif



#define FMM_ST_DIPOLE 5

#if (FMM_PRECISION==0)

typedef struct {
	float x[4];
	float y[4];
	float z[4];
	float q[4];
	float mx[4];
	float my[4];
	float mz[4];
} _DATA_dipole_ST_V0_t;

typedef struct {
	float tx[4];
	float ty[4];
	float tz[4];
	float pot[4];
} _DATA_dipole_ST_V1_t;

typedef struct {
	_DATA_dipole_ST_V0_t *_0; /* x , y , z , q , mx , my , mz  */
	_DATA_dipole_ST_V1_t *_1; /* tx , ty , tz , pot  */
} DATA_dipole_ST_t;

#elif (FMM_PRECISION==1)

typedef struct {
	double x[4];
	double y[4];
	double z[4];
	double q[4];
	double mx[4];
	double my[4];
	double mz[4];
} _DATA_dipole_ST_V0_t;

typedef struct {
	double tx[4];
	double ty[4];
	double tz[4];
	double pot[4];
} _DATA_dipole_ST_V1_t;

typedef struct {
	_DATA_dipole_ST_V0_t *_0; /* x , y , z , q , mx , my , mz  */
	_DATA_dipole_ST_V1_t *_1; /* tx , ty , tz , pot  */
} DATA_dipole_ST_t;
#endif



#define FMM_ST_GRAD 6

#if (FMM_PRECISION==0)

typedef struct {
	float x[4];
	float y[4];
	float z[4];
	float q[4];
} _DATA_grad_ST_V0_t;

typedef struct {
	float tx[4];
	float ty[4];
	float tz[4];
	float pot[4];
	float gradx[4];
	float grady[4];
	float gradz[4];
} _DATA_grad_ST_V1_t;

typedef struct {
	_DATA_grad_ST_V0_t *_0; /* x , y , z , q  */
	_DATA_grad_ST_V1_t *_1; /* tx , ty , tz , pot , gradx , grady , gradz  */
} DATA_grad_ST_t;

#elif (FMM_PRECISION==1)

typedef struct {
	double x[4];
	double y[4];
	double z[4];
	double q[4];
} _DATA_grad_ST_V0_t;

typedef struct {
	double tx[4];
	double ty[4];
	double tz[4];
	double pot[4];
	double gradx[4];
	double grady[4];
	double gradz[4];
} _DATA_grad_ST_V1_t;

typedef struct {
	_DATA_grad_ST_V0_t *_0; /* x , y , z , q  */
	_DATA_grad_ST_V1_t *_1; /* tx , ty , tz , pot , gradx , grady , gradz  */
} DATA_grad_ST_t;
#endif



#define FMM_ST_DIPOLE_GRAD 7

#if (FMM_PRECISION==0)

typedef struct {
	float x[4];
	float y[4];
	float z[4];
	float q[4];
	float mx[4];
	float my[4];
	float mz[4];
} _DATA_dipole_grad_ST_V0_t;

typedef struct {
	float tx[4];
	float ty[4];
	float tz[4];
	float pot[4];
	float gradx[4];
	float grady[4];
	float gradz[4];
} _DATA_dipole_grad_ST_V1_t;

typedef struct {
	_DATA_dipole_grad_ST_V0_t *_0; /* x , y , z , q , mx , my , mz  */
	_DATA_dipole_grad_ST_V1_t *_1; /* tx , ty , tz , pot , gradx , grady , gradz  */
} DATA_dipole_grad_ST_t;

#elif (FMM_PRECISION==1)

typedef struct {
	double x[4];
	double y[4];
	double z[4];
	double q[4];
	double mx[4];
	double my[4];
	double mz[4];
} _DATA_dipole_grad_ST_V0_t;

typedef struct {
	double tx[4];
	double ty[4];
	double tz[4];
	double pot[4];
	double gradx[4];
	double grady[4];
	double gradz[4];
} _DATA_dipole_grad_ST_V1_t;

typedef struct {
	_DATA_dipole_grad_ST_V0_t *_0; /* x , y , z , q , mx , my , mz  */
	_DATA_dipole_grad_ST_V1_t *_1; /* tx , ty , tz , pot , gradx , grady , gradz  */
} DATA_dipole_grad_ST_t;
#endif



#define FMM_ADDITIONAL_DATA_IN_BOX\

#define FMM_ADDITIONAL_GLOBAL_DATA\
	_FLOAT_ (*dipoleMoments)[3];\
	_FLOAT_ *charges;\
	_FLOAT_ *potentials;\
	_FLOAT_ (*gradients)[3];\
	_FLOAT_ beta;	\
	_FLOAT_ *lambda;	\
	_FLOAT_ *w;	\
	int *M;	\
	int *FFT_rep;	\
	int len_F;	\
	_FLOAT_ **Tz_M2M;	\
	_FLOAT_ **Tz_L2L;	\
	_FLOAT_ *Tz_M2M_coeffs;	\
	_FLOAT_ *Tz_L2L_coeffs;	\
	_FLOAT_ *CMX_coeffs;	\
	_FLOAT_ *CXL_coeffs;	\
	_FLOAT_ *CXL_reduced_coeffs;	\
	_FLOAT_ **Ry_theta;	\
	_FLOAT_ **Ry_minus_theta;	\
	_FLOAT_ **Ry_pi_minus_theta;	\
	_FLOAT_ **Ry_minus_pi_minus_theta;	\
	_FLOAT_ **Ry_pi2;	\
	_FLOAT_ **Ry_minus_pi2;	\
	_FLOAT_ *Ry_theta_coeffs;	\
	_FLOAT_ *Ry_minus_theta_coeffs;	\
	_FLOAT_ *Ry_pi_minus_theta_coeffs;	\
	_FLOAT_ *Ry_minus_pi_minus_theta_coeffs;	\
	_FLOAT_ *Ry_pi2_coeffs;	\
	_FLOAT_ *Ry_minus_pi2_coeffs;	\
	int *P_MRT;	\
	int *P_MRT_plus;	\
	int *P_MRT_minus;	\
	int *P_Mriri2rrii;	\
	int *P_LRT;	\
	int *P_LRT_plus;	\
	int *P_LRT_minus;	\
	int *P_Lriri2rrii;	\
	int *P_X_riri2rrii;	\
	int *P_VF;	\
	int *neg_F;	\
	int len_neg_F;	\
	int *P_FV;	\
	int *neg_V;	\
	int len_neg_V;	\

#include"fmmv3d.h"
#include"fmmv_common.h"
#if (FMM_PRECISION==0)
  #include"simd.h"
#elif (FMM_PRECISION==1)
  #include"simd.h"
#endif /* FMM_PRECISION */

void gen_M_standard(FmmvHandle *FMMV, Box *box);
void gen_M_dipole(FmmvHandle *FMMV, Box *box);
void gen_M_grad(FmmvHandle *FMMV, Box *box);
void gen_M_dipole_grad(FmmvHandle *FMMV, Box *box);
void gen_M_ST_standard(FmmvHandle *FMMV, Box *box);
void gen_M_ST_dipole(FmmvHandle *FMMV, Box *box);
void gen_M_ST_grad(FmmvHandle *FMMV, Box *box);
void gen_M_ST_dipole_grad(FmmvHandle *FMMV, Box *box);
void eval_L_standard(FmmvHandle *FMMV, Box *box);
void eval_L_dipole(FmmvHandle *FMMV, Box *box);
void eval_L_grad(FmmvHandle *FMMV, Box *box);
void eval_L_dipole_grad(FmmvHandle *FMMV, Box *box);
void eval_L_ST_standard(FmmvHandle *FMMV, Box *box);
void eval_L_ST_dipole(FmmvHandle *FMMV, Box *box);
void eval_L_ST_grad(FmmvHandle *FMMV, Box *box);
void eval_L_ST_dipole_grad(FmmvHandle *FMMV, Box *box);
void eval_M_standard(FmmvHandle *FMMV, Box *target, Box *source);
void eval_M_dipole(FmmvHandle *FMMV, Box *target, Box *source);
void eval_M_grad(FmmvHandle *FMMV, Box *target, Box *source);
void eval_M_dipole_grad(FmmvHandle *FMMV, Box *target, Box *source);
void eval_M_ST_standard(FmmvHandle *FMMV, Box *target, Box *source);
void eval_M_ST_dipole(FmmvHandle *FMMV, Box *target, Box *source);
void eval_M_ST_grad(FmmvHandle *FMMV, Box *target, Box *source);
void eval_M_ST_dipole_grad(FmmvHandle *FMMV, Box *target, Box *source);
void gen_L_standard(FmmvHandle *FMMV, Box *target, Box *source);
void gen_L_dipole(FmmvHandle *FMMV, Box *target, Box *source);
void gen_L_grad(FmmvHandle *FMMV, Box *target, Box *source);
void gen_L_dipole_grad(FmmvHandle *FMMV, Box *target, Box *source);
void gen_L_ST_standard(FmmvHandle *FMMV, Box *target, Box *source);
void gen_L_ST_dipole(FmmvHandle *FMMV, Box *target, Box *source);
void gen_L_ST_grad(FmmvHandle *FMMV, Box *target, Box *source);
void gen_L_ST_dipole_grad(FmmvHandle *FMMV, Box *target, Box *source);
void gen_L_eval_M_standard(FmmvHandle *FMMV, Box *list3, Box *list4);
void gen_L_eval_M_dipole(FmmvHandle *FMMV, Box *list3, Box *list4);
void gen_L_eval_M_grad(FmmvHandle *FMMV, Box *list3, Box *list4);
void gen_L_eval_M_dipole_grad(FmmvHandle *FMMV, Box *list3, Box *list4);
void extrinsic_correction_standard(FmmvHandle *FMMV);
void extrinsic_correction_dipole(FmmvHandle *FMMV);
void extrinsic_correction_grad(FmmvHandle *FMMV);
void extrinsic_correction_dipole_grad(FmmvHandle *FMMV);
void extrinsic_correction_ST_standard(FmmvHandle *FMMV);
void extrinsic_correction_ST_dipole(FmmvHandle *FMMV);
void extrinsic_correction_ST_grad(FmmvHandle *FMMV);
void extrinsic_correction_ST_dipole_grad(FmmvHandle *FMMV);
void eval_direct_standard_acc0(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_standard_acc1(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_standard_acc2(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_dipole_acc0(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_dipole_acc1(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_dipole_acc2(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_grad_acc0(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_grad_acc1(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_grad_acc2(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_dipole_grad_acc0(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_dipole_grad_acc1(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_dipole_grad_acc2(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_ST_standard_acc0(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_ST_standard_acc1(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_ST_standard_acc2(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_ST_dipole_acc0(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_ST_dipole_acc1(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_ST_dipole_acc2(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_ST_grad_acc0(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_ST_grad_acc1(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_ST_grad_acc2(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_ST_dipole_grad_acc0(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_ST_dipole_grad_acc1(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_ST_dipole_grad_acc2(FmmvHandle *FMMV, Box *target, Box *source);
void eval_direct_periodic_standard_acc0(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_standard_acc1(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_standard_acc2(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_dipole_acc0(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_dipole_acc1(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_dipole_acc2(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_grad_acc0(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_grad_acc1(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_grad_acc2(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_dipole_grad_acc0(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_dipole_grad_acc1(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_dipole_grad_acc2(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_ST_standard_acc0(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_ST_standard_acc1(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_ST_standard_acc2(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_ST_dipole_acc0(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_ST_dipole_acc1(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_ST_dipole_acc2(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_ST_grad_acc0(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_ST_grad_acc1(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_ST_grad_acc2(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_ST_dipole_grad_acc0(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_ST_dipole_grad_acc1(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
void eval_direct_periodic_ST_dipole_grad_acc2(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);
#endif /* _FMMV_H_ */
