/*
 * FMMV - the Fastest Multipole Method of Vienna
 * Copyright (c) 2006 Harald Hofstaetter
 * Institute for Analysis and Scientific Computing
 * Vienna University of Technology
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

//#include <Accelerate/Accelerate.h>

#define V4_INIT /* NONE */

#undef V4_LITTLE_ENDIAN /* big endian */

typedef vector float V4_TYPE;      
typedef float V4_BASETYPE;

#define V4_ZERO      ((vector float) (0.0f))
#define V4_NEGZERO      ((vector float) (-0.0f))

#define V4_ADD(x, y) vec_add(x, y)
#define V4_SUB(x, y) vec_sub(x, y)
#define V4_MUL(x, y) vec_madd(x, y, V4_NEGZERO)

static inline vector float V4_REC(vector float x)
{
    vector float est = vec_re(x);
    return vec_madd(vec_nmsub(est, x, (vector float) (1.0)), est, est);
}

#define V4_DIV(x, y) V4_MUL(x, V4_REC(y))
#define V4_XOR(x, y) vec_xor(x, y)
#define V4_NEG(x)    vec_xor(x, V4_NEGZERO)
#define V4_FMA(x, y, z) vec_madd(x, y, z)
//#define V4_FMS(x, y, z) V4_SUB(V4_MUL(x, y), z)
#define V4_FMS(x, y, z) V4_NEG(vec_nmsub(x, y, z))
#define V4_FNMS(x, y, z) vec_nmsub(x, y, z)

#define V4_SET(x, y, z, w)   ((vector float) (x, y, z, w))
#define V4_SET1(x)          vec_splat((vector float)(x), 0)
#define V4_LOAD(x)           vec_ld(0, x)
#define V4_LOAD1(x)          vec_splat((vector float)*(x), 0)
#define V4_STORE(x, y)       vec_st(y, 0, x)
#define V4_LOAD0(x)          vec_perm((vector float) 0.0f, (vector float)*(x), (vector unsigned char) ( 16, 17, 18, 19, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15))
#define V4_STORE0(x, y)      vec_ste(y, 0, x)
#define V4_GETITEM(x, i)     V4_LOAD((x)+(i))
#define V4_GETBASEITEM(x, i) V4_LOAD1((x)+(i))
#define V4_SETITEM(x, i, y)  V4_STORE((x)+(i), y)


#define V4_RECIP_SQRT0(x) vec_rsqrte(x)

static inline V4_TYPE V4_RECIP_SQRT1(V4_TYPE x)
{
        V4_TYPE est = V4_RECIP_SQRT0(x);
        V4_TYPE est2 = V4_MUL(est, est);
        V4_TYPE esth = V4_MUL(est, V4_SET1(0.5));

        return V4_FMA(V4_FNMS(x, est2, V4_SET1(1.0)), esth, est);
}

static inline V4_TYPE V4_RECIP_SQRT2(V4_TYPE x)
{
        V4_TYPE one_half = V4_SET1(0.5);
        V4_TYPE three = V4_SET1(3.0);
        V4_TYPE est = V4_RECIP_SQRT0(x);

        est = V4_MUL(one_half, V4_MUL(est, V4_FNMS(x, V4_MUL(est, est), three)));
        return V4_MUL(one_half, V4_MUL(est, V4_FNMS(x, V4_MUL(est, est), three)));
}

#define V4_RECIP_SQRT(x) V4_RECIP_SQRT1(x)

#define V4_SQRT(x) V4_MUL(x, V4_RECIP_SQRT(x))


static inline vector float V4_HORIZADD(vector float v0, vector float v1, vector float v2, vector float v3)
{
        vector unsigned char low = (vector unsigned char)
                         (0, 1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 20, 21, 22, 23);
        vector unsigned char high = (vector unsigned char)
                         (8, 9, 10, 11, 12, 13, 14, 15, 24, 25, 26, 27, 28, 29, 30, 31);
        vector unsigned char even = (vector unsigned char)
                         (0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19, 24, 25, 26, 27);
        vector unsigned char odd = (vector unsigned char)
                         (4, 5, 6, 7, 12, 13, 14, 15, 20, 21, 22, 23, 28, 29, 30, 31);
        vector float h0 = vec_add(vec_perm(v0, v1, low), vec_perm(v0, v1, high));
        vector float h1 = vec_add(vec_perm(v2, v3, low), vec_perm(v2, v3, high));
        return vec_add(vec_perm(h0, h1, even), vec_perm(h0, h1, odd));
}



static __inline__ V4_TYPE V4_HORIZADD1(V4_TYPE x)
{
	x = vec_add(x, vec_sld(x,x,8));
	x = vec_add(x, vec_sld(x,x,4));
	return x;
}	

#endif /* _SIMD4S_H_ */
