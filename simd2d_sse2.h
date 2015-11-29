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

#ifndef _SIMD2D_H_
#define _SIMD2D_H_

#include<emmintrin.h>
typedef __m128d V2_TYPE;      
//__m128d vmldInvSqrt2(__m128d a);

#define V2_INIT /* NONE */

typedef double V2_BASETYPE;

#define V2_ZERO      _mm_setzero_pd()
#define V2_NEGZERO   _mm_set1_pd(-0.0)

#define V2_ADD(x, y) _mm_add_pd(x, y)
#define V2_SUB(x, y) _mm_sub_pd(x, y)
#define V2_MUL(x, y) _mm_mul_pd(x, y)
#define V2_DIV(x, y) _mm_div_pd(x, y)
#define V2_XOR(x, y) _mm_xor_pd(x, y)
#define V2_NEG(x)    _mm_xor_pd(x, V2_NEGZERO)
#define V2_FMA(x, y, z) _mm_add_pd(_mm_mul_pd(x, y), z)
#define V2_FMS(x, y, z) _mm_sub_pd(_mm_mul_pd(x, y), z)
#define V2_FNMS(x, y, z) _mm_sub_pd(z, _mm_mul_pd(x, y))


#define V2_SET(x, y)   _mm_set_pd(x, y)
#define V2_SET1(x)           _mm_set1_pd(x)
#define V2_LOAD(x)           _mm_load_pd(x)
#define V2_LOAD1(x)          _mm_load1_pd(x)
#define V2_LOAD2(x, y)	     _mm_loadh_pd(_mm_load_sd(x), y)
#define V2_STORE(x, y)       _mm_store_pd(x, y)
#define V2_LOAD0(x)          _mm_load_sd(x)
#define V2_STORE0(x, y)      _mm_store_sd(x, y)


#define V2_GETITEM(x, i)     _mm_load_pd((x)+(i))
#define V2_GETBASEITEM(x, i) _mm_load1_pd((x)+(i))
#define V2_SETITEM(x, i, y)  _mm_store_pd((x)+(i), y)

static __inline void V2_STORE2(double *x, double *y, __m128d z) 
{
	_mm_store_sd(x, z); 
	_mm_storeh_pd(y, z);
}	
	
static __inline __m128d V2_HORIZADD(__m128d v0, __m128d v1) 
{
	return _mm_add_pd(_mm_shuffle_pd(v0, v1, _MM_SHUFFLE2(1, 1)), 
			  _mm_shuffle_pd(v0, v1, _MM_SHUFFLE2(0, 0)));

}

static __inline __m128d V2_HORIZADD1(__m128d x)
{
	return _mm_add_sd(x, _mm_shuffle_pd(x, x, _MM_SHUFFLE2(1, 1)));
}	

#define V2_SQRT(x)   _mm_sqrt_pd(x)

//#define V2_RECIP_SQRT(x) vmldInvSqrt2(x)

#define V2_RECIP_SQRT0(x) _mm_cvtps_pd(_mm_rsqrt_ps(_mm_cvtpd_ps(x)))


static __inline V2_TYPE V2_RECIP_SQRT1(V2_TYPE x)
{
        V2_TYPE est = V2_RECIP_SQRT0(x);
        V2_TYPE est2 = V2_MUL(est, est);
        V2_TYPE esth = V2_MUL(est, V2_SET1(0.5));
        return V2_FMA(V2_FNMS(x, est2, V2_SET1(1.0)), esth, est);
}


static __inline V2_TYPE V2_RECIP_SQRT2(V2_TYPE x)
{
	V2_TYPE one_half = V2_SET1(0.5);
	V2_TYPE three = V2_SET1(3.0);
        V2_TYPE est = V2_RECIP_SQRT0(x);

	est = V2_MUL(one_half, V2_MUL(est, V2_FNMS(x, V2_MUL(est, est), three)));
	return V2_MUL(one_half, V2_MUL(est, V2_FNMS(x, V2_MUL(est, est), three)));
}

#define V2_RECIP_SQRT(x) V2_RECIP_SQRT2(x)

#endif /* _SIMD2D_ */
