#ifndef __FMMV2D_H_
#define __FMMV2D_H_
#include "fmmv_sys.h"

int init_all(FmmvHandle *FMMV);
void init_coeffs(FmmvHandle *FMMV);
void init_quad_coeffs(FmmvHandle *FMMV);
void double2single(int n, double *source, float *target);
void finish_coeffs(FmmvHandle *FMMV);
void finish_quad_coeffs(FmmvHandle *FMMV);

void M2X(FmmvHandle *FMMV,  Box *box, _FLOAT_ *X_N, _FLOAT_ *X_S, _FLOAT_ *X_W, _FLOAT_ *X_E);
void X2L(FmmvHandle *FMMV,  Box *box);

void direct_method_xp( unsigned int NParticles, _FLOAT_ particles[][2], _FLOAT_ charges[], _FLOAT_ dipoleMoments[][2],
        unsigned int NTargets, _FLOAT_ targets[][2], _FLOAT_ pot[], _FLOAT_ grad[][2], /* result of length NTarget */
        double beta, double *time, char **errorMessage);

void direct_method_complex_xp( unsigned int NParticles, _FLOAT_ particles[][2],
        _FLOAT_ complexCharges[][2], _FLOAT_ complexDipoleMoments[][2],
        unsigned int NTargets, _FLOAT_ targets[][2], 
        _FLOAT_ complexPotentials[][2],   _FLOAT_ complexGradients[][2],
        double beta, double *time, char **errorMessage);

float k0f(float x);
float k1f(float x);
float i0f(float x);
float i1f(float x);
double k0(double x);
double k1(double x);
double i0(double x);
double i1(double x);

#include<math.h> /* sqrt, sqrtf, etc. */
#if (FMM_PRECISION==0)
  #define RECIP0(x) (1.0/(x))
  #define RECIP1(x) (1.0/(x))
  #define RECIP2(x) (1.0/(x))
  #define RECIP(x) RECIP2(x)
  #define LOG0(x) (logf(x))
  #define LOG1(x) (logf(x))
  #define LOG2(x) (logf(x))
  #define LOG(x) LOG2(x)
  #define K00(x) (k0f(x))
  #define K01(x) (k0f(x))
  #define K02(x) (k0f(x))
  #define K0(x) K02(x)
  #define K10(x) (k1f(x))
  #define K11(x) (k1f(x))
  #define K12(x) (k1f(x))
  #define K1(x) K12(x)
  #define I00(x) (i0f(x))
  #define I01(x) (i0f(x))
  #define I02(x) (i0f(x))
  #define I0(x) I02(x)
  #define I10(x) (i1f(x))
  #define I11(x) (i1f(x))
  #define I12(x) (i1f(x))
  #define I1(x) I12(x)
  #define SQRT(x) sqrtf(x)
  #define ATAN20(y, x) atan2f(y, x)
  #define ATAN21(y, x) atan2f(y, x)
  #define ATAN22(y, x) atan2f(y, x)
  #define ATAN2(y, x) ATAN22(y, x)
#else  
  #define RECIP0(x) (1.0/(x))
  #define RECIP1(x) (1.0/(x))
  #define RECIP2(x) (1.0/(x))
  #define RECIP(x) RECIP2(x)
  #define LOG0(x) (log(x))
  #define LOG1(x) (log(x))
  #define LOG2(x) (log(x))
  #define LOG(x) LOG2(x)
  #define K00(x) (k0(x))
  #define K01(x) (k0(x))
  #define K02(x) (k0(x))
  #define K0(x) K02(x)
  #define K10(x) (k1(x))
  #define K11(x) (k1(x))
  #define K12(x) (k1(x))
  #define K1(x) K12(x)
  #define I00(x) (i0(x))
  #define I01(x) (i0(x))
  #define I02(x) (i0(x))
  #define I0(x) I02(x)
  #define I10(x) (i1(x))
  #define I11(x) (i1(x))
  #define I12(x) (i1(x))
  #define I1(x) I12(x)
  #define SQRT(x) sqrt(x)
  #define ATAN20(y, x) atan2(y, x)
  #define ATAN21(y, x) atan2(y, x)
  #define ATAN22(y, x) atan2(y, x)
  #define ATAN2(y, x) ATAN22(y, x)
#endif


#endif /* __FMMV2D_H_ */


