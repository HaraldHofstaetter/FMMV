/*
 * FMMV - the Fastest Multipole Method of Vienna
 * Copyright (c) 2006-2015 Harald Hofstaetter
 * http://www.harald-hofstaetter.at
 * 
 * This file is part of FMMV.
 * 
 * FMMV is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * FMMV is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with FMMV; if not, write to the Free Software  Foundation, Inc.,
 * 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 * 
 */

#include "_fmmv.h"
#include<math.h>
#include "simd.h"
#if (FMM_PRECISION==0)
   #include"simd4s.h"
#else
   #include"simd4d.h"
#endif


// NOTE: order of vector components!
/*
#if (FMM_PRECISION==0)
#define V4_EXP(x) V4_SET(expf(((V4_BASETYPE*) &(x))[3]), \
                         expf(((V4_BASETYPE*) &(x))[2]), \
                         expf(((V4_BASETYPE*) &(x))[1]), \
                         expf(((V4_BASETYPE*) &(x))[0]))

#define V4_SINH(x) V4_SET(sinhf(((V4_BASETYPE*) &(x))[3]), \
                          sinhf(((V4_BASETYPE*) &(x))[2]), \
                          sinhf(((V4_BASETYPE*) &(x))[1]), \
                          sinhf(((V4_BASETYPE*) &(x))[0]))
#else
#define V4_EXP(x) V4_SET(exp(((V4_BASETYPE*) &(x))[3]), \
                         exp(((V4_BASETYPE*) &(x))[2]), \
                         exp(((V4_BASETYPE*) &(x))[1]), \
                         exp(((V4_BASETYPE*) &(x))[0]))

#define V4_SINH(x) V4_SET(sinh(((V4_BASETYPE*) &(x))[3]), \
                          sinh(((V4_BASETYPE*) &(x))[2]), \
                          sinh(((V4_BASETYPE*) &(x))[1]), \
                          sinh(((V4_BASETYPE*) &(x))[0]))
#endif
*/

//#include"sse_mathfun.h"
__m128 exp_ps(__m128 x); 
#define V4_EXP(x) exp_ps(x)

V4_TYPE V4_SINH(V4_TYPE x)
{
   V4_TYPE ex = V4_EXP(x);
   return V4_MUL(V4_SET1(0.5),V4_SUB(ex, V4_DIV(V4_SET1(1.0),ex)));
}   

void bessel_i_simd4(V4_TYPE x, V4_BASETYPE *i, int n, int n1)
{
     int l;
     V4_TYPE one_over_x, s;

     //one_over_x = 1.0/x;
     one_over_x = V4_DIV(V4_SET1(1.0) ,x);
     //i[n1-1] = 0.0;
     V4_SETITEM(i, 4*(n1-1), V4_SET1(0.0));
     //i[n1-2] = 1.0;
     V4_SETITEM(i, 4*(n1-2), V4_SET1(1.0));

     for (l=n1-3; l>=0; l--) {
         //i[l] = i[l+2] + (2*l+3)*i[l+1]*one_over_x;
         V4_SETITEM(i, 4*l, V4_ADD(V4_GETITEM(i, 4*(l+2)),  
                    V4_MUL(V4_MUL(V4_SET1(2*l+3), V4_GETITEM(i, 4*(l+1))), one_over_x)));
     }

     //s = sinh(x)*one_over_x/i[0];
     s = V4_DIV(V4_MUL(V4_SINH(x), one_over_x), V4_GETITEM(i, 0));
     for (l=0;l<n; l++) {
         //i[l] *= s;
         V4_SETITEM(i, 4*l, V4_MUL(V4_GETITEM(i, 4*l), s));
     }    
}

void bessel_i_scaled_simd4(V4_TYPE x, V4_BASETYPE *i, int n, int n1)
/* computes i_n(x)/x^(n+1) */
{
     int l;
     V4_TYPE s, x2;
    
     //x2 = x*x;
     x2 = V4_MUL(x, x);
     //i[n1-1] = 0.0;
     V4_SETITEM(i, 4*(n1-1), V4_SET1(0.0));
     //i[n1-2] = 1.0;
     V4_SETITEM(i, 4*(n1-2), V4_SET1(1.0));

     for (l=n1-3; l>=0; l--) {
#if (FMM_PRECISION==0)
         if (l==n-3) {
              V4_SETITEM(i, 4*(l+2), V4_DIV(V4_GETITEM(i, 4*(l+2)), V4_GETITEM(i, 4*(l+1))));
              V4_SETITEM(i, 4*(l+1), V4_SET1(1.0));
         }     
#endif     
         //i[l] = x2*i[l+2] + (2*l+3)*i[l+1];
         V4_SETITEM(i, 4*l, V4_ADD(V4_MUL(x2,V4_GETITEM(i, 4*(l+2))),  
                    V4_MUL(V4_SET1(2*l+3), V4_GETITEM(i, 4*(l+1)))));
     }

     //s = sinh(x)/(x2*i[0]);
     s = V4_DIV(V4_SINH(x), V4_MUL(x2, V4_GETITEM(i, 0)));
     for (l=0;l<n; l++) {
         //i[l] *= s;
         V4_SETITEM(i, 4*l, V4_MUL(V4_GETITEM(i, 4*l), s));
     }    
}


void bessel_k_simd4(V4_TYPE x, V4_BASETYPE *k, int n)
{
     int l;
     V4_TYPE mx, one_over_x;

     //one_over_x = 1.0/x;
     one_over_x = V4_DIV(V4_SET1(1.0) ,x);
     //k[0] = exp(-x)*one_over_x; 
     mx = V4_NEG(x);
     V4_SETITEM(k, 0, V4_MUL(V4_EXP(mx),one_over_x));
     //k[1] = k[0]*one_over_x*(x+1.0);
     V4_SETITEM(k, 4, V4_MUL(V4_MUL(V4_GETITEM(k, 0), one_over_x), V4_ADD(x, V4_SET1(1.0))));
     for (l=2; l<n; l++) {
         //k[l] = + k[l-2] + (2*l-1)*k[l-1]*one_over_x;
         V4_SETITEM(k, 4*l, V4_ADD(V4_GETITEM(k, 4*(l-2)),  
                    V4_MUL(V4_MUL(V4_SET1(2*l-1), V4_GETITEM(k, 4*(l-1))), one_over_x)));
     }
}     

void bessel_k_scaled_simd4(V4_TYPE x, V4_BASETYPE *k, int n)
/* computes k_n(x)*x^(n+1) */
{
     int l;
     V4_TYPE mx, x2;

     //x2 = x*x;
     x2 = V4_MUL(x, x);
     //k[0] = exp(-x);
     mx = V4_NEG(x);
     V4_SETITEM(k, 0, V4_EXP(mx));
     //k[1] = k[0]*(x+1.0);
     V4_SETITEM(k, 4, V4_MUL(V4_GETITEM(k, 0), V4_ADD(x, V4_SET1(1.0))));
     for (l=2; l<n; l++) {
         //k[l] = x2*k[l-2] + (2*l-1)*k[l-1];
         V4_SETITEM(k, 4*l, V4_ADD(V4_MUL(x2,V4_GETITEM(k, 4*(l-2))),  
                    V4_MUL(V4_SET1(2*l-1), V4_GETITEM(k, 4*(l-1)))));
     }
}
