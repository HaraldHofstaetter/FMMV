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

#ifndef _SIMD4D_H_
#define _SIMD4D_H_

#include<immintrin.h>

#define V4_LITTLE_ENDIAN
#define V4_INIT /* NONE */

typedef __m256d V4_TYPE;      
typedef double V4_BASETYPE;

#define V4_ZERO      _mm256_setzero_pd()
#define V4_NEGZERO   _mm256_set1_pd(-0.0)

#define V4_ADD(x, y) _mm256_add_pd(x, y)
#define V4_SUB(x, y) _mm256_sub_pd(x, y)
#define V4_MUL(x, y) _mm256_mul_pd(x, y)
#define V4_DIV(x, y) _mm256_div_pd(x, y)
#define V4_XOR(x, y) _mm256_xor_pd(x, y)
#define V4_NEG(x)    _mm256_xor_pd(x, V4_NEGZERO)
#define V4_FMA(x, y, z) _mm256_add_pd(_mm256_mul_pd(x, y), z)
#define V4_FMS(x, y, z) _mm256_sub_pd(_mm256_mul_pd(x, y), z)
#define V4_FNMS(x, y, z) _mm256_sub_pd(z, _mm256_mul_pd(x, y))

#define V4_SET(x, y, z, w)   _mm256_set_pd(x, y, z, w)
#define V4_SET1(x)           _mm256_set1_pd(x)
#define V4_LOAD(x)           _mm256_load_pd(x)
#define V4_LOAD1(x)          _mm256_broadcast_sd(x)
#define V4_STORE(x, y)       _mm256_store_pd(x, y)
#define V4_LOAD0(x)          _mm256_castpd128_pd256(_mm_load_sd(x))  /* CHECK THIS */
#define V4_STORE0(x, y)      _mm_store_sd(x, _mm256_castpd256_pd128(y))

#define V4_GETITEM(x, i)     _mm256_load_pd((x)+(i))
#define V4_GETBASEITEM(x, i) _mm256_broadcast_sd((x)+(i))
#define V4_SETITEM(x, i, y)  _mm256_store_pd((x)+(i), y)

#define V4_SQRT(x)   _mm256_sqrt_pd(x)
#define V4_RECIP_SQRT0(x) _mm256_cvtps_pd(_mm_rsqrt_ps(_mm256_cvtpd_ps(x)))

//#include<xmmintrin.h>
//__m128 vmlsInvSqrt4(__m128 a)
//#define V4_RECIP_SQRT(x) vmlsInvSqrt4(x)

static __inline V4_TYPE V4_RECIP_SQRT1(V4_TYPE x)
{
        /* one Newton iteration */
        V4_TYPE est = V4_RECIP_SQRT0(x);
        V4_TYPE est2 = V4_MUL(est, est);
        V4_TYPE esth = V4_MUL(est, V4_SET1(0.5));
        return V4_FMA(V4_FNMS(x, est2, V4_SET1(1.0)), esth, est);
}

static __inline V4_TYPE V4_RECIP_SQRT2(V4_TYPE x)
{
        /* two Newton iterations */
        V4_TYPE one_half = V4_SET1(0.5);
        V4_TYPE three = V4_SET1(3.0);
        V4_TYPE est = V4_RECIP_SQRT0(x);

        est = V4_MUL(one_half, V4_MUL(est, V4_FNMS(x, V4_MUL(est, est), three)));
        return V4_MUL(one_half, V4_MUL(est, V4_FNMS(x, V4_MUL(est, est), three)));
}

//#define V4_RECIP_SQRT(x) V4_RECIP_SQRT2(x)
#define V4_RECIP_SQRT(x) V4_DIV(V4_SET1(1.0),V4_SQRT(x))


#define V4_RECIP0(x) _mm256_cvtps_pd(_mm_rcp_ps(_mm256_cvtpd_ps(x)))

static __inline V4_TYPE V4_RECIP1(V4_TYPE x)
{
        /* one Newton iteration */
        V4_TYPE est = V4_RECIP0(x);
        return V4_MUL(est, V4_FNMS(est, x, V4_SET1(2.0)));
}

static __inline V4_TYPE V4_RECIP2(V4_TYPE x)
{
        /* two Newton iterations */
        V4_TYPE est = V4_RECIP0(x);
        V4_TYPE est1 = V4_MUL(est, V4_FNMS(est, x, V4_SET1(2.0)));
        return V4_MUL(est1, V4_FNMS(est1, x, V4_SET1(2.0)));
}

//#define V4_RECIP(x) V4_RECIP2(x)
#define V4_RECIP(x) V4_DIV(V4_SET1(1.0),(x))


static __inline __m256d V4_HORIZADD(__m256d a, __m256d b, __m256d c, __m256d d) 
{
/* From: answer to question
 * http://stackoverflow.com/questions/10833234/4-horizontal-double-precision-sums-in-one-go-with-avx
 * by Norbert P.
 */   
   /* {a[0]+a[1], b[0]+b[1], a[2]+a[3], b[2]+b[3]} */
    __m256d sumab = _mm256_hadd_pd(a, b);
   /* {c[0]+c[1], d[0]+d[1], c[2]+c[3], d[2]+d[3]} */
    __m256d sumcd = _mm256_hadd_pd(c, d);
   /* {a[0]+a[1], b[0]+b[1], c[2]+c[3], d[2]+d[3]} */
    __m256d blend = _mm256_blend_pd(sumab, sumcd, 0b1100);
   /* {a[2]+a[3], b[2]+b[3], c[0]+c[1], d[0]+d[1]} */
    __m256d perm = _mm256_permute2f128_pd(sumab, sumcd, 0x21);

    return _mm256_add_pd(perm, blend); 
}

static __inline __m256d V4_HORIZADD1(__m256d x)
{
        __m128d x128 = _mm_add_pd(_mm256_extractf128_pd(x, 1), _mm256_castpd256_pd128(x)); 
        return _mm256_castpd128_pd256(_mm_add_pd(x128, _mm_shuffle_pd(x128, x128, 0x55)));   
        //return _mm256_castpd128_pd256(_mm_hadd_pd(x128, _mm_setzero_pd()));   //0x55));
}

#endif /* _SIMD4D_H_ */
