#include"_fmmv.h"
#include<math.h>

#if (FMM_PRECISION==0)
   #define EXP(x) expf(x)
   #define SINH(x) sinhf(x)
#else
   #define EXP(x) exp(x)
   #define SINH(x) sinh(x)
#endif

void bessel_i(_FLOAT_ x, _FLOAT_ i[], int n, int n1)
{
     int l;
     _FLOAT_ one_over_x, s;

     one_over_x = 1.0/x;
     i[n1-1] = 0.0;
     i[n1-2] = 1.0;

     for (l=n1-3; l>=0; l--) {
         i[l] = i[l+2] + (2*l+3)*i[l+1]*one_over_x;
     }

     s = SINH(x)*one_over_x/i[0];
     for (l=0;l<n; l++) {
         i[l] *= s;
     }    
}

void bessel_i_scaled(_FLOAT_ x, _FLOAT_ i[], int n, int n1)
/* computes i_n(x)/x^(n+1) */
{
     int l;
     _FLOAT_ s, x2;

     x2 = x*x;
     i[n1-1] = 0.0;
     i[n1-2] = 1.0;

     for (l=n1-3; l>=0; l--) {
#if (FMM_PRECISION==0)
         if (l==n-3) {
              i[l+2]/=i[l+1];
              i[l+1]=1.0;
         }     
#endif         
         i[l] = x2*i[l+2] + (2*l+3)*i[l+1];
     }

     s = SINH(x)/(x2*i[0]);
     for (l=0;l<n; l++) {
         i[l] *= s;
     }    
}


void bessel_k(_FLOAT_ x, _FLOAT_ k[], int n)
{
     int l;
     _FLOAT_ one_over_x;

     one_over_x = 1.0/x;
     k[0] = EXP(-x)*one_over_x; // *1.5707963267948966192; /* pi/2 */
     k[1] = k[0]*one_over_x*(x+1.0);
     for (l=2; l<n; l++) {
         k[l] = (2*l-1)*k[l-1]*one_over_x + k[l-2];
     }
}     

void bessel_k_scaled(_FLOAT_ x, _FLOAT_ k[], int n)
/* computes k_n(x)*x^(n+1) */
{
     int l;
     _FLOAT_ x2;

     x2 = x*x;
     k[0] = EXP(-x);
     k[1] = k[0]*(x+1.0);
     for (l=2; l<n; l++) {
         k[l] = (2*l-1)*k[l-1] + x2*k[l-2];
     }
}     

#define Re(n,m) (2*((n)*((n)+1)/2 + (m)))
#define Im(n,m) (2*((n)*((n)+1)/2 + (m))+1)

void core_eval_L_grad_yukawa(int p, _FLOAT_ *rho, _FLOAT_ *scaleL, _FLOAT_ *L, _FLOAT_ *Y, _FLOAT_ *x, _FLOAT_ *y, _FLOAT_ *z)
{

    int n,m;
    _FLOAT_ hx, hy ,hz;
    _FLOAT_ cx_m,cy_m,cz_m;
    _FLOAT_ cx_p,cy_p,cz_p;
    _FLOAT_ dxY_re_m, dxY_im_m, dyY_re_m, dyY_im_m, dzY_re_m, dzY_im_m;
    _FLOAT_ dxY_re_p, dxY_im_p, dyY_re_p, dyY_im_p, dzY_re_p, dzY_im_p;

    hx = 0.0;
    hy = 0.0;
    hz = 0.0;

    for (n=0;n<=p;n++) {
      cx_m = 0.0;
      cy_m = 0.0;
      cz_m = 0.0;
      cx_p = 0.0;
      cy_p = 0.0;
      cz_p = 0.0;
      for (m=0;m<=n;m++) {
        if (n==0) {
                dxY_re_m = 0.0;
                dxY_re_p = -sqrt(2)*Y[Re(n+1,1)];
                dxY_im_m = 0.0;
                dxY_im_p = 0;
                dyY_re_m = 0.0;
                dyY_re_p = -sqrt(2)*Y[Im(n+1,1)];
                dyY_im_m = 0.0;
                dyY_im_p = 0.0;
                dzY_re_m = 0.0;
                dzY_re_p = Y[Re(n+1,0)];
                dzY_im_m = 0.0;
                dzY_im_p = 0.0; 
        }
        else {
             if (m==0) {        
		dxY_re_m =  sqrt(n*(n-1))*Y[Re(n-1,1)];
		dxY_re_p = -sqrt((n+1)*(n+2))*Y[Re(n+1,1)];
                dxY_im_m =  0;
                dxY_im_p =  0;
		dyY_re_m =  sqrt(n*(n-1))*Y[Im(n-1,1)];
		dyY_re_p = -sqrt((n+1)*(n+2))*Y[Im(n+1,1)];
                dyY_im_m =  0;
                dyY_im_p =  0;
                dzY_re_m =  n*Y[Re(n-1,0)];
                dzY_re_p =  (n+1)*Y[Re(n+1,0)];
                dzY_im_m = 0.0; 
                dzY_im_p = 0.0; 
             }
             else {
		dxY_re_m = +sqrt((n-m)*(n-m-1))*Y[Re(n-1,m+1)] - sqrt((n+m)*(n+m-1))*Y[Re(n-1,m-1)]; 
		dxY_re_p = -sqrt((n+m+1)*(n+m+2))*Y[Re(n+1,m+1)] + sqrt((n-m+1)*(n-m+2))*Y[Re(n+1,m-1)];
		dxY_im_m = +sqrt((n-m)*(n-m-1))*Y[Im(n-1,m+1)] - sqrt((n+m)*(n+m-1))*Y[Im(n-1,m-1)];
		dxY_im_p = -sqrt((n+m+1)*(n+m+2))*Y[Im(n+1,m+1)] + sqrt((n-m+1)*(n-m+2))*Y[Im(n+1,m-1)]; 
		dyY_re_m = +sqrt((n-m)*(n-m-1))*Y[Im(n-1,m+1)] + sqrt((n+m)*(n+m-1))*Y[Im(n-1,m-1)];
		dyY_re_p = -sqrt((n+m+1)*(n+m+2))*Y[Im(n+1,m+1)] - sqrt((n-m+1)*(n-m+2))*Y[Im(n+1,m-1)]; 
		dyY_im_m = -sqrt((n-m)*(n-m-1))*Y[Re(n-1,m+1)] -sqrt((n+m)*(n+m-1))*Y[Re(n-1,m-1)]; 
		dyY_im_p = + sqrt((n+m+1)*(n+m+2))*Y[Re(n+1,m+1)] + sqrt((n-m+1)*(n-m+2))*Y[Re(n+1,m-1)]; 
                dzY_re_m = 2.0*sqrt((n-m)*(n+m))*Y[Re(n-1,m)];
                dzY_re_p = 2.0*sqrt((n+1-m)*(n+1+m))*Y[Re(n+1,m)];
                dzY_im_m = 2.0*sqrt((n-m)*(n+m))*Y[Im(n-1,m)];
                dzY_im_p = 2.0*sqrt((n+1-m)*(n+1+m))*Y[Im(n+1,m)];
             }
       }
       cx_m += dxY_re_m*L[Re(n,m)] - dxY_im_m*L[Im(n,m)]; 
       cx_p += dxY_re_p*L[Re(n,m)] - dxY_im_p*L[Im(n,m)]; 
       cy_m += dyY_re_m*L[Re(n,m)] - dyY_im_m*L[Im(n,m)]; 
       cy_p += dyY_re_p*L[Re(n,m)] - dyY_im_p*L[Im(n,m)]; 
       cz_m += dzY_re_m*L[Re(n,m)] - dzY_im_m*L[Im(n,m)]; 
       cz_p += dzY_re_p*L[Re(n,m)] - dzY_im_p*L[Im(n,m)]; 
     }
     if (n==0) {
	     hx += scaleL[n]*rho[n+1]*cx_p;
	     hy += scaleL[n]*rho[n+1]*cy_p;
	     hz += scaleL[n]*rho[n+1]*cz_p;
     }
     else{
	     hx += scaleL[n]*(rho[n-1]*cx_m + rho[n+1]*cx_p);///((double) (2*n+1));
	     hy += scaleL[n]*(rho[n-1]*cy_m + rho[n+1]*cy_p);///((double) (2*n+1));
	     hz += scaleL[n]*(rho[n-1]*cz_m + rho[n+1]*cz_p);///((double) (2*n+1));
     }	     
   }
   *x = hx;
   *y = hy;
   *z = hz;
}


void core_eval_L_grad_minus(int p, _FLOAT_ *scale, _FLOAT_ *L, _FLOAT_ *Y, _FLOAT_ *x, _FLOAT_ *y, _FLOAT_ *z)
{

    int n,m;
    _FLOAT_ hx, hy ,hz;
    _FLOAT_ cx,cy,cz;
    _FLOAT_ dx_re, dx_im, dy_re, dy_im, dz_re, dz_im;

    hx = 0.0;
    hy = 0.0;
    hz = 0.0;

    //for (n=0;n<=p;n++) {
    for (n=1;n<=p;n++) {
      cx = 0.0;
      cy = 0.0;
      cz = 0.0;
      for (m=0;m<=n;m++) {
             if (m==0) {        
		dx_re =  sqrt(n*(n-1))*Y[Re(n-1,1)];
                dx_im =  0;
		dy_re =  sqrt(n*(n-1))*Y[Im(n-1,1)];
                dy_im =  0;
                dz_re =  n*Y[Re(n-1,0)];
                dz_im = 0.0; 
             }
             else {
		dx_re = +sqrt((n-m)*(n-m-1))*Y[Re(n-1,m+1)] - sqrt((n+m)*(n+m-1))*Y[Re(n-1,m-1)]; 
		dx_im = +sqrt((n-m)*(n-m-1))*Y[Im(n-1,m+1)] - sqrt((n+m)*(n+m-1))*Y[Im(n-1,m-1)];
		dy_re = +sqrt((n-m)*(n-m-1))*Y[Im(n-1,m+1)] + sqrt((n+m)*(n+m-1))*Y[Im(n-1,m-1)];
		dy_im = -sqrt((n-m)*(n-m-1))*Y[Re(n-1,m+1)] - sqrt((n+m)*(n+m-1))*Y[Re(n-1,m-1)]; 
                dz_re = 2.0*sqrt((n-m)*(n+m))*Y[Re(n-1,m)];
                dz_im = 2.0*sqrt((n-m)*(n+m))*Y[Im(n-1,m)];
             }
             cx += dx_re*L[Re(n,m)] - dx_im*L[Im(n,m)]; 
             cy += dy_re*L[Re(n,m)] - dy_im*L[Im(n,m)]; 
             cz += dz_re*L[Re(n,m)] - dz_im*L[Im(n,m)]; 
     }
     hx += scale[n]*cx;
     hy += scale[n]*cy;
     hz += scale[n]*cz;
   }
   *x = hx;
   *y = hy;
   *z = hz;
}



void core_eval_L_grad_plus(int p, _FLOAT_ *scale, _FLOAT_ *L, _FLOAT_ *Y, _FLOAT_ *x, _FLOAT_ *y, _FLOAT_ *z)
{

    int n,m;
    _FLOAT_ hx, hy ,hz;
    _FLOAT_ cx,cy,cz;
    _FLOAT_ dx_re, dx_im, dy_re, dy_im, dz_re, dz_im;

    hx = 0.0;
    hy = 0.0;
    hz = 0.0;

    for (n=0;n<=p;n++) {
      cx = 0.0;
      cy = 0.0;
      cz = 0.0;
      for (m=0;m<=n;m++) {
             if (m==0) {        
		dx_re = -sqrt((n+1)*(n+2))*Y[Re(n+1,1)];
                dx_im =  0;
		dy_re = -sqrt((n+1)*(n+2))*Y[Im(n+1,1)];
                dy_im =  0;
                dz_re =  (n+1)*Y[Re(n+1,0)];
                dz_im = 0.0; 
             }
             else {
		dx_re = -sqrt((n+m+1)*(n+m+2))*Y[Re(n+1,m+1)] + sqrt((n-m+1)*(n-m+2))*Y[Re(n+1,m-1)];
		dx_im = -sqrt((n+m+1)*(n+m+2))*Y[Im(n+1,m+1)] + sqrt((n-m+1)*(n-m+2))*Y[Im(n+1,m-1)]; 
		dy_re = -sqrt((n+m+1)*(n+m+2))*Y[Im(n+1,m+1)] - sqrt((n-m+1)*(n-m+2))*Y[Im(n+1,m-1)]; 
		dy_im = +sqrt((n+m+1)*(n+m+2))*Y[Re(n+1,m+1)] + sqrt((n-m+1)*(n-m+2))*Y[Re(n+1,m-1)]; 
                dz_re = 2.0*sqrt((n+1-m)*(n+1+m))*Y[Re(n+1,m)];
                dz_im = 2.0*sqrt((n+1-m)*(n+1+m))*Y[Im(n+1,m)];
             }
             cx += dx_re*L[Re(n,m)] - dx_im*L[Im(n,m)]; 
             cy += dy_re*L[Re(n,m)] - dy_im*L[Im(n,m)]; 
             cz += dz_re*L[Re(n,m)] - dz_im*L[Im(n,m)]; 
     }
     hx += scale[n]*cx;
     hy += scale[n]*cy;
     hz += scale[n]*cz;
   }
   *x = hx;
   *y = hy;
   *z = hz;
}
