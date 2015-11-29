#ifndef _FMMV3D_H_
#define _FMMV3D_H_

#define X_RRII
#include "fmmv_sys.h"

/* macros for indexing: */
#define I(n,m) ((n)*(n)+(n)+(m))
#define J(n,m) ((n)*((n)+1)/2 + (m))

int init_all(FmmvHandle *fh);
void init_coeffs(FmmvHandle *FMMV);
void init_quad_coeffs(FmmvHandle *FMMV3D);
void init_F(FmmvHandle *FMMV);
void init_perm(FmmvHandle *FMMV);
void init_Ry(FmmvHandle *FMMV);
void double2single(int n, double *source, float *target);
int finish_all(FmmvHandle *fh);
void finish_coeffs(FmmvHandle *FMMV);
void finish_quad_coeffs(FmmvHandle *FMMV);
void finish_F(FmmvHandle *FMMV);
void finish_perm(FmmvHandle *FMMV);
void finish_Ry(FmmvHandle *FMMV);

void M2X(FmmvHandle *FMMV, int dir, Box *box, _FLOAT_ *X1, _FLOAT_ *X2);
void X2L(FmmvHandle *FMMV, int dir, Box *box, int reduced);

void Ry(int p, _FLOAT_ **bocks, _FLOAT_ *x, _FLOAT_ *y);
void Ry_pi(int p, _FLOAT_ *y);

void Rz_pi2(int p, _FLOAT_ *x, _FLOAT_ *y);
void Rz_minus_pi2(int p, _FLOAT_ *x, _FLOAT_ *y);
void Rz_pi4(int p, _FLOAT_ *x, _FLOAT_ *y);
void Rz_minus_pi4(int p, _FLOAT_ *x, _FLOAT_ *y);
void Rz_3pi4(int p, _FLOAT_ *x, _FLOAT_ *y);
void Rz_minus_3pi4(int p, _FLOAT_ *x, _FLOAT_ *y);


void F_M2X(FmmvHandle *fh,_FLOAT_ *x, _FLOAT_ *y);
void F_X2L(FmmvHandle *fh,_FLOAT_ *x, _FLOAT_ *y);

/* SIMD4 ****************/
void Ry_simd4(int p, _FLOAT_ **bocks, _FLOAT_ *x, _FLOAT_ *y);
void Ry_pi_simd4(int p, _FLOAT_ *y);

void Rz_simd4(int p, _FLOAT_ *x, _FLOAT_ *y, int back);
void Rz1_simd4(int p, _FLOAT_ *x, _FLOAT_ *y, int back);
void Rz_pi2_rrii_simd4(int p, _FLOAT_ *x, _FLOAT_ *y);
void Rz_minus_pi2_rrii_simd4(int p, _FLOAT_ *x, _FLOAT_ *y);



void F_M2X_simd4(FmmvHandle *fh,_FLOAT_ *x, _FLOAT_ *y);
void F_X2L_simd4(FmmvHandle *fh,_FLOAT_ *x, _FLOAT_ *y);

/* SIMD2 ****************/
void Ry_simd2(int p, _FLOAT_ **bocks, _FLOAT_ *x, _FLOAT_ *y);
void Ry_pi_simd2(int p, _FLOAT_ *y);

void Rz_simd2(int p, _FLOAT_ *x, _FLOAT_ *y, int k1, int k2, int back);
void Rz1_simd2(int p, _FLOAT_ *x, _FLOAT_ *y, int k1, int k2, int back);
void Rz_pi2_rrii_simd2(int p, _FLOAT_ *x, _FLOAT_ *y);
void Rz_minus_pi2_rrii_simd2(int p, _FLOAT_ *x, _FLOAT_ *y);


void F_M2X_simd2(FmmvHandle *fh,_FLOAT_ *x, _FLOAT_ *y);
void F_X2L_simd2(FmmvHandle *fh,_FLOAT_ *x, _FLOAT_ *y);

/********************/


void gemv_simd2(int m, int n, _FLOAT_ *A, int lda, _FLOAT_ *x, _FLOAT_ *y);
void gemv_trans_simd2(int m, int n, _FLOAT_ *A, int lda, _FLOAT_ *x, _FLOAT_ *y);
void tpmv_upper_simd2(int n, _FLOAT_ *A, _FLOAT_ *x);
void tpmv_lower_simd2(int n, _FLOAT_ *A, _FLOAT_ *x);
_FLOAT_ *vec2_add2_simd2(FmmvHandle *fh, const int n, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_2add2_simd2(FmmvHandle *fh, int n, const _FLOAT_ *x, _FLOAT_ **y1, _FLOAT_ **y2);
void vec2_copy_simd2(const int n, const _FLOAT_ *x1, const _FLOAT_ *x2, _FLOAT_ *y);
void P_2riri2rrii_simd2(int p, _FLOAT_ *x1, _FLOAT_ *x2, _FLOAT_ *y);
void P_rrii2riri2_simd2(FmmvHandle *fh, int p, const _FLOAT_ *x, _FLOAT_ **y1, _FLOAT_ **y2, _FLOAT_ *mem);
void P_X_2rrii2riri_simd2(int s, _FLOAT_ *x1, _FLOAT_ *x2, _FLOAT_ *y);
void P_X_riri2rrii2_simd2(int s, const _FLOAT_ *x, _FLOAT_ **y1, _FLOAT_ **y2);
void scale_X_simd2(int p, _FLOAT_ *xx, int level);
void perm(const int n, const int *vec, const _FLOAT_ *x, _FLOAT_ *y);
void perm_inv(const int n, const int *vec, const _FLOAT_ *x, _FLOAT_ *y);
void perm_simd2(const int n, const int *vec, const _FLOAT_ *x, _FLOAT_ *y);
void perm_inv_simd2(const int n, const int *vec, const _FLOAT_ *x, _FLOAT_ *y);

void gemv_simd4(int m, int n, _FLOAT_ *A, int lda, _FLOAT_ *x, _FLOAT_ *y);
void gemv_trans_simd4(int m, int n, _FLOAT_ *A, int lda, _FLOAT_ *x, _FLOAT_ *y);
void tpmv_upper_simd4(int n, _FLOAT_ *A, _FLOAT_ *x);
void tpmv_lower_simd4(int n, _FLOAT_ *A, _FLOAT_ *x);
_FLOAT_ *vec4_add2_simd4(FmmvHandle *fh, const int n, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_4add2_simd4(FmmvHandle *fh, int n, const _FLOAT_ *x, _FLOAT_ **y1, _FLOAT_ **y2, _FLOAT_ **y3, _FLOAT_ **y4);
void vec4_copy_simd4(const int n, const _FLOAT_ *x1, const _FLOAT_ *x2, const _FLOAT_ *x3, const _FLOAT_ *x4, _FLOAT_ *y);
void P_4riri2rrii_simd4(int p, _FLOAT_ *x1, _FLOAT_ *x2, _FLOAT_ *x3, _FLOAT_ *x4, _FLOAT_ *y);
void P_rrii2riri4_simd4(FmmvHandle *fh, int p, const _FLOAT_ *x, _FLOAT_ **y1, _FLOAT_ **y2, _FLOAT_ **y3, _FLOAT_ **y4, _FLOAT_ *mem);
void P_X_4rrii2riri_simd4(int s, _FLOAT_ *x1, _FLOAT_ *x2, _FLOAT_ *x3, _FLOAT_ *x4, _FLOAT_ *y);
void P_X_riri2rrii4_simd4(int s, const _FLOAT_ *x, _FLOAT_ **y1, _FLOAT_ **y2, _FLOAT_ **y3, _FLOAT_ **y4);
void scale_X_simd4(int p, _FLOAT_ *xx, int level);
void perm_simd4(const int n, const int *vec, const _FLOAT_ *x, _FLOAT_ *y);
void perm_inv_simd4(const int n, const int *vec, const _FLOAT_ *x, _FLOAT_ *y);

void compute_spherical_harmonics_generic(int p, _FLOAT_ *Y, _FLOAT_ sin_phi, _FLOAT_ cos_phi, _FLOAT_ cos_theta);

void core_gen_M_L(int p, _FLOAT_ *rho, _FLOAT_ q, _FLOAT_ *dML, _FLOAT_ *ML);
_FLOAT_ core_eval_L_M(int p, _FLOAT_ *rho, _FLOAT_ *LM, _FLOAT_ *Y);
void core_eval_L_M_grad_minus(int p, _FLOAT_ *scale, _FLOAT_ *L, _FLOAT_ *Y, _FLOAT_ *x, _FLOAT_ *y, _FLOAT_ *z);
void core_eval_L_M_grad_plus(int p, _FLOAT_ *scale, _FLOAT_ *M, _FLOAT_ *Y, _FLOAT_ *x, _FLOAT_ *y, _FLOAT_ *z);
void core_gen_M_L_dipole_minus(int p, _FLOAT_ *scale, _FLOAT_ mx, _FLOAT_ my, _FLOAT_ mz,  _FLOAT_ *Y, _FLOAT_ *M);
void core_gen_M_L_dipole_plus(int p, _FLOAT_ *scale, _FLOAT_ mx, _FLOAT_ my, _FLOAT_ mz,  _FLOAT_ *Y, _FLOAT_ *M);


extern _FLOAT_ _ONE_F_[];
extern _FLOAT_ _ZERO_F_[];
extern int _ZERO_[];

void gemv(int m, int n, _FLOAT_ *A, int lda, _FLOAT_ *x, _FLOAT_ *y);
void gemv_trans(int m, int n, _FLOAT_ *A, int lda, _FLOAT_ *x, _FLOAT_ *y);
void tpmv_upper(int n, _FLOAT_ *A, _FLOAT_ *x);
void tpmv_lower(int n, _FLOAT_ *A, _FLOAT_ *x);

void direct_method_xp( unsigned int NParticles, _FLOAT_ particles[][3], _FLOAT_ charges[], _FLOAT_ dipoleMoments[][3],
        unsigned int NTargets, _FLOAT_ targets[][3], _FLOAT_ pot[], _FLOAT_ grad[][3], 
        double beta, double *time, char **errorMessage);


#include<math.h> /* sqrt, sqrtf */
#if (FMM_PRECISION==0)
  #define SQRT(x) sqrtf(x)
  #define RECIP_SQRT0(x) (1.0/sqrtf(x))
  #define RECIP_SQRT1(x) (1.0/sqrtf(x))
  #define RECIP_SQRT2(x) (1.0/sqrtf(x))
  #define RECIP_SQRT(x) RECIP_SQRT1(x)
  #define EXP0(x) (expf(x))
  #define EXP1(x) (expf(x))
  #define EXP2(x) (expf(x))
  #define EXP(x) EXP1(x) 

/*
#ifdef USE_CBLAS
 #include "cblas.h"

  #define gemv(m,n,a,lda,x,y) \
	cblas_sgemv(CblasRowMajor,CblasNoTrans,m,n,1.0,a,lda,x,1,0.0,y,1)

  #define gemv_trans(m,n,a,lda,x,y) \
	cblas_sgemv(CblasRowMajor,CblasTrans,m,n,1.0,a,lda,x,1,0.0,y,1)

  #define tpmv_upper(n,a,x) \
	cblas_stpmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,n,a,x,1)
  
  #define tpmv_lower(n,a,x) \
	cblas_stpmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,n,a,x,1)

#else
  extern int sgemv_(char *, int *, int *, float *, 
	            float *, int *, float *, int *, float *, float *, int *);	
  extern int stpmv_(char *, char *, char *, int *, float *, float *, int *);

  #define gemv(m,n,a,lda,x,y) \
	sgemv_("T",_ZERO_+n, _ZERO_+m, _ONE_F_, a, _ZERO_+lda,x,_ZERO_+1,_ZERO_F_,y,_ZERO_+1)	

  #define gemv_trans(m,n,a,lda,x,y) \
	sgemv_("N",_ZERO_+n, _ZERO_+m, _ONE_F_, a, _ZERO_+lda,x,_ZERO_+1,_ZERO_F_,y,_ZERO_+1)	

  #define tpmv_upper(n,a,x) \
	stpmv_("L","T","N",_ZERO_+n,a,x,_ZERO_+1)
  
  #define tpmv_lower(n,a,x) \
	stpmv_("U","T","N",_ZERO_+n,a,x,_ZERO_+1)
  
#endif
*/

#else  
  #define SQRT(x) sqrt(x)
  #define RECIP_SQRT0(x) (1.0/sqrt(x))
  #define RECIP_SQRT1(x) (1.0/sqrt(x))
  #define RECIP_SQRT2(x) (1.0/sqrt(x))
  #define RECIP_SQRT(x) RECIP_SQRT2(x)
  #define EXP0(x) (exp(x))
  #define EXP1(x) (exp(x))
  #define EXP2(x) (exp(x))
  #define EXP(x) EXP2(x) 

/*
#ifdef USE_CBLAS
 #include "cblas.h"

  #define gemv(m,n,a,lda,x,y) \
	cblas_dgemv(CblasRowMajor,CblasNoTrans,m,n,1.0,a,lda,x,1,0.0,y,1)

  #define gemv_trans(m,n,a,lda,x,y) \
	cblas_dgemv(CblasRowMajor,CblasTrans,m,n,1.0,a,lda,x,1,0.0,y,1)

  #define tpmv_upper(n,a,x) \
	cblas_dtpmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,n,a,x,1)
  
  #define tpmv_lower(n,a,x) \
	cblas_dtpmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit,n,a,x,1)

#else
  extern int dgemv_(char *, int *, int *, double *, 
	            double *, int *, double *, int *, double *, double *, int *);	
  extern int dtpmv_(char *, char *, char *, int *, double *, double *, int *);

  #define gemv(m,n,a,lda,x,y) \
	dgemv_("T",_ZERO_+n, _ZERO_+m, _ONE_F_, a, _ZERO_+lda,x,_ZERO_+1,_ZERO_F_,y,_ZERO_+1)	

  #define gemv_trans(m,n,a,lda,x,y) \
	dgemv_("N",_ZERO_+n, _ZERO_+m, _ONE_F_, a, _ZERO_+lda,x,_ZERO_+1,_ZERO_F_,y,_ZERO_+1)	

  #define tpmv_upper(n,a,x) \
	dtpmv_("L","T","N",_ZERO_+n,a,x,_ZERO_+1)
  
  #define tpmv_lower(n,a,x) \
	dtpmv_("U","T","N",_ZERO_+n,a,x,_ZERO_+1)

#endif
*/
#endif

#endif /* _FMMV3D_H_ */

