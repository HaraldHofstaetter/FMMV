/*
 * FMMV - the Fastest Multipole Method of Vienna
 * Copyright (c) 2006-2010 Harald Hofstaetter
 * Institute Mathematics
 * University of Vienna, Austria
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
   #include"simd2s.h"
#else
   #include"simd2d.h"
#endif

// NOTE: order of vector components!
#if (FMM_PRECISION==0)
#define V2_EXP(x) V2_SET(expf(((V2_BASETYPE*) &(x))[1]), \
                         expf(((V2_BASETYPE*) &(x))[0]))

#define V2_SINH(x) V2_SET(sinhf(((V2_BASETYPE*) &(x))[1]), \
                          sinhf(((V2_BASETYPE*) &(x))[0]))
#else
#define V2_EXP(x) V2_SET(exp(((V2_BASETYPE*) &(x))[1]), \
                         exp(((V2_BASETYPE*) &(x))[0]))

#define V2_SINH(x) V2_SET(sinh(((V2_BASETYPE*) &(x))[1]), \
                          sinh(((V2_BASETYPE*) &(x))[0]))
#endif

void bessel_i_simd2(V2_TYPE x, V2_BASETYPE *i, int n, int n1)
{
     int l;
     V2_TYPE one_over_x, s;

     //one_over_x = 1.0/x;
     one_over_x = V2_DIV(V2_SET1(1.0) ,x);
     //i[n1-1] = 0.0;
     V2_SETITEM(i, 2*(n1-1), V2_SET1(0.0));
     //i[n1-2] = 1.0;
     V2_SETITEM(i, 2*(n1-2), V2_SET1(1.0));

     for (l=n1-3; l>=0; l--) {
         //i[l] = i[l+2] + (2*l+3)*i[l+1]*one_over_x;
         V2_SETITEM(i, 2*l, V2_ADD(V2_GETITEM(i, 2*(l+2)),  
                    V2_MUL(V2_MUL(V2_SET1(2*l+3), V2_GETITEM(i, 2*(l+1))), one_over_x)));
     }

     //s = sinh(x)*one_over_x/i[0];
     s = V2_DIV(V2_MUL(V2_SINH(x), one_over_x), V2_GETITEM(i, 0));
     for (l=0;l<n; l++) {
         //i[l] *= s;
         V2_SETITEM(i, 2*l, V2_MUL(V2_GETITEM(i, 2*l), s));
     }    
}

void bessel_i_scaled_simd2(V2_TYPE x, V2_BASETYPE *i, int n, int n1)
/* computes i_n(x)/x^(n+1) */
{
     int l;
     V2_TYPE s, x2;
    
     //x2 = x*x;
     x2 = V2_MUL(x, x);
     //i[n1-1] = 0.0;
     V2_SETITEM(i, 2*(n1-1), V2_SET1(0.0));
     //i[n1-2] = 1.0;
     V2_SETITEM(i, 2*(n1-2), V2_SET1(1.0));

     for (l=n1-3; l>=0; l--) {
         //i[l] = x2*i[l+2] + (2*l+3)*i[l+1];
         V2_SETITEM(i, 2*l, V2_ADD(V2_MUL(x2,V2_GETITEM(i, 2*(l+2))),  
                    V2_MUL(V2_SET1(2*l+3), V2_GETITEM(i, 2*(l+1)))));
     }

     //s = sinh(x)/(x2*i[0]);
     s = V2_DIV(V2_SINH(x), V2_MUL(x2, V2_GETITEM(i, 0)));
     for (l=0;l<n; l++) {
         //i[l] *= s;
         V2_SETITEM(i, 2*l, V2_MUL(V2_GETITEM(i, 2*l), s));
     }    
}


void bessel_k_simd2(V2_TYPE x, V2_BASETYPE *k, int n)
{
     int l;
     V2_TYPE mx, one_over_x;

     //one_over_x = 1.0/x;
     one_over_x = V2_DIV(V2_SET1(1.0) ,x);
     //k[0] = exp(-x)*one_over_x; 
     mx = V2_NEG(x);
     V2_SETITEM(k, 0, V2_MUL(V2_EXP(mx),one_over_x));
     //k[1] = k[0]*one_over_x*(x+1.0);
     V2_SETITEM(k, 2, V2_MUL(V2_MUL(V2_GETITEM(k, 0), one_over_x), V2_ADD(x, V2_SET1(1.0))));
     for (l=2; l<n; l++) {
         //k[l] = + k[l-2] + (2*l-1)*k[l-1]*one_over_x;
         V2_SETITEM(k, 2*l, V2_ADD(V2_GETITEM(k, 2*(l-2)),  
                    V2_MUL(V2_MUL(V2_SET1(2*l-1), V2_GETITEM(k, 2*(l-1))), one_over_x)));
     }
}     

void bessel_k_scaled_simd2(V2_TYPE x, V2_BASETYPE *k, int n)
/* computes k_n(x)*x^(n+1) */
{
     int l;
     V2_TYPE mx, x2;

     //x2 = x*x;
     x2 = V2_MUL(x, x);
     //k[0] = exp(-x);
     mx = V2_NEG(x);
     V2_SETITEM(k, 0, V2_EXP(mx));
     //k[1] = k[0]*(x+1.0);
     V2_SETITEM(k, 2, V2_MUL(V2_GETITEM(k, 0), V2_ADD(x, V2_SET1(1.0))));
     for (l=2; l<n; l++) {
         //k[l] = x2*k[l-2] + (2*l-1)*k[l-1];
         V2_SETITEM(k, 2*l, V2_ADD(V2_MUL(x2,V2_GETITEM(k, 2*(l-2))),  
                    V2_MUL(V2_SET1(2*l-1), V2_GETITEM(k, 2*(l-1)))));
     }
}
