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

#ifndef _SIMD4S_H_
#define _SIMD4S_H_

#define V4_LITTLE_ENDIAN
#define V4_INIT /* NONE */

typedef __m128 V4_TYPE;      
typedef float V4_BASETYPE;

#define V4_ZERO      _mm_setzero_ps()
#define V4_NEGZERO   _mm_set1_ps(-0.0)

#define V4_ADD(x, y) _mm_add_ps(x, y)
#define V4_SUB(x, y) _mm_sub_ps(x, y)
#define V4_MUL(x, y) _mm_mul_ps(x, y)
#define V4_DIV(x, y) _mm_div_ps(x, y)
#define V4_XOR(x, y) _mm_xor_ps(x, y)
#define V4_NEG(x)    _mm_xor_ps(x, V4_NEGZERO)
#define V4_FMA(x, y, z) _mm_add_ps(_mm_mul_ps(x, y), z)
#define V4_FMS(x, y, z) _mm_sub_ps(_mm_mul_ps(x, y), z)
#define V4_FNMS(x, y, z) _mm_sub_ps(z, _mm_mul_ps(x, y))

#define V4_SET(x, y, z, w)   _mm_set_ps(x, y, z, w)
#define V4_SET1(x)           _mm_set1_ps(x)
#define V4_LOAD(x)           _mm_load_ps(x)
#define V4_LOAD1(x)          _mm_load_ps1(x)
#define V4_STORE(x, y)       _mm_store_ps(x, y)
#define V4_LOAD0(x)          _mm_load_ss(x)
#define V4_STORE0(x, y)      _mm_store_ss(x, y)
#define V4_GETITEM(x, i)     _mm_load_ps((x)+(i))
#define V4_GETBASEITEM(x, i) _mm_load_ps1((x)+(i))
#define V4_SETITEM(x, i, y)  _mm_store_ps((x)+(i), y)

#define V4_SQRT(x)   _mm_sqrt_ps(x)
#define V4_RECIP_SQRT0(x) _mm_rsqrt_ps(x)

//#include<xmmintrin.h>
//__m128 vmlsInvSqrt4(__m128 a)
//#define V4_RECIP_SQRT(x) vmlsInvSqrt4(x)

static __inline V4_TYPE V4_RECIP_SQRT1(V4_TYPE x)
{
        V4_TYPE est = V4_RECIP_SQRT0(x);
        V4_TYPE est2 = V4_MUL(est, est);
        V4_TYPE esth = V4_MUL(est, V4_SET1(0.5));
        return V4_FMA(V4_FNMS(x, est2, V4_SET1(1.0)), esth, est);
}

static __inline V4_TYPE V4_RECIP_SQRT2(V4_TYPE x)
{
        V4_TYPE one_half = V4_SET1(0.5);
        V4_TYPE three = V4_SET1(3.0);
        V4_TYPE est = V4_RECIP_SQRT0(x);

        est = V4_MUL(one_half, V4_MUL(est, V4_FNMS(x, V4_MUL(est, est), three)));
        return V4_MUL(one_half, V4_MUL(est, V4_FNMS(x, V4_MUL(est, est), three)));
}

#define V4_RECIP_SQRT(x) V4_RECIP_SQRT1(x)

#define V4_RECIP0(x) _mm_rcp_ps(x)
#define V4_RECIP1(x) V4_DIV(V4_SET1(1.0),(x))
#define V4_RECIP2(x) V4_DIV(V4_SET1(1.0),(x))
#define V4_RECIP(x) V4_RECIP(x)


#ifdef _MSC_VER
static __inline __m128 V4_HORIZADD0(__m128 *v0, __m128 *v1, __m128 *v2, __m128 *v3) 
{
__m128 h0=_mm_add_ps(_mm_movelh_ps(*v0, *v1), _mm_movehl_ps(*v1, *v0));
__m128 h1=_mm_add_ps(_mm_movelh_ps(*v2, *v3), _mm_movehl_ps(*v3, *v2)); 
return _mm_add_ps(_mm_shuffle_ps(h0, h1, 0xDD), _mm_shuffle_ps(h0, h1, 0x88));
}	

#define V4_HORIZADD(v0, v1, v2, v3) V4_HORIZADD0(&(v0), &(v1), &(v2), &(v3)) 

#else
static __inline __m128 V4_HORIZADD(__m128 v0, __m128 v1, __m128 v2, __m128 v3) 
{
__m128 h0=_mm_add_ps(_mm_movelh_ps(v0, v1), _mm_movehl_ps(v1, v0));
__m128 h1=_mm_add_ps(_mm_movelh_ps(v2, v3), _mm_movehl_ps(v3, v2)); 
return _mm_add_ps(_mm_shuffle_ps(h0, h1, 0xDD), _mm_shuffle_ps(h0, h1, 0x88));
}
#endif


static __inline __m128 V4_HORIZADD1(__m128 x)
{
	return _mm_add_ss(x, 
	_mm_add_ss(_mm_movehl_ps(x, x),
	_mm_add_ss(_mm_shuffle_ps(x, x, _MM_SHUFFLE(1,1,1,1)),	
        _mm_shuffle_ps(x, x, _MM_SHUFFLE(3,3,3,3)))));
   
}


static __inline __m128 V4_LOAD_LEFT(void *i, int mask)
{
    __m128 x = _mm_load_ps(i);
    switch (mask) {
        case 0:
            return x;
        case 1:
            return _mm_shuffle_ps(x, x, _MM_SHUFFLE(3,2,1,1));
        case 2:
            return _mm_shuffle_ps(x, x, _MM_SHUFFLE(3,2,2,2));
        case 3:
            return _mm_shuffle_ps(x, x, _MM_SHUFFLE(3,3,3,3));
    }    
}


static __inline __m128 V4_LOAD_LEFT0(void *i, int mask)
{
    __m128 x = _mm_load_ps(i);
    switch (mask) {
        case 0:
            return x;
        case 1: {
            __m128i m = {(u_int64_t) 0xFFFFFFFF00000000, (u_int64_t) 0xFFFFFFFFFFFFFFFF}; 
            return _mm_and_ps(x, (__m128) m); }
        case 2: {
            __m128i m = {(u_int64_t) 0x0, (u_int64_t) 0xFFFFFFFFFFFFFFFF}; 
            return _mm_and_ps(x, (__m128) m); }
        case 3: {
            __m128i m = {(u_int64_t) 0x0, (u_int64_t) 0xFFFFFFFF00000000}; 
            return _mm_and_ps(x, (__m128) m); }
    }    
}


static __inline __m128 V4_LOAD_RIGHT(void *i, int mask)
{
    __m128 x = _mm_load_ps(i);
    switch (mask) {
        case 0:
            return x;
        case 1:
            return _mm_shuffle_ps(x, x, _MM_SHUFFLE(0,0,0,0));
        case 2:
            return _mm_shuffle_ps(x, x, _MM_SHUFFLE(1,1,1,0));
        case 3:
            return _mm_shuffle_ps(x, x, _MM_SHUFFLE(2,2,1,0));
    }    
}


static __inline __m128 V4_LOAD_RIGHT0(void *i, int mask)
{
    __m128 x = _mm_load_ps(i);
    switch (mask) {
        case 0: 
            return x;
        case 1: {
            __m128i m = {(u_int64_t) 0x00000000FFFFFFFF, (u_int64_t) 0x0}; 
            return _mm_and_ps(x, (__m128) m); }
        case 2: {
            __m128i m = {(u_int64_t) 0xFFFFFFFFFFFFFFFF, (u_int64_t) 0x0}; 
            return _mm_and_ps(x, (__m128) m); }
        case 3: {
            __m128i m = {(u_int64_t) 0xFFFFFFFFFFFFFFFF, (u_int64_t) 0x00000000FFFFFFFF}; 
            return _mm_and_ps(x, (__m128) m); }
    }    
}



#endif /* _SIMD4S_H_ */
