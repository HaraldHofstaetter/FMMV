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

#ifndef _SIMD_H_
#define _SIMD_H_

#include<xmmintrin.h>

#define FMM_SIMD_ALIGNMENT (sizeof(_FLOAT_)*FMM_SIMD_VECTOR_LENGTH)

/*
void* _mm_malloc(size_t size, size_t alignment); 
void _mm_free(void * ptr);
*/

#undef SIMD_ALIGN
#ifdef _MSC_VER
  #define SIMD_ALIGN __declspec(align(FMM_SIMD_ALIGNMENT))
#else
  #define SIMD_ALIGN __attribute__ ((aligned(FMM_SIMD_ALIGNMENT)))
#endif

#undef FMMV_MALLOC
#define FMMV_MALLOC(FMMV, n) _mm_malloc(n, FMM_SIMD_ALIGNMENT); \
        FMMV->allocatedMemory+=(n); \
        if (FMMV->allocatedMemory>FMMV->maxAllocatedMemory) \
                FMMV->maxAllocatedMemory=FMMV->allocatedMemory;

#undef FMMV_FREE
#define FMMV_FREE(FMMV, x, n) if (x) FMMV->allocatedMemory-=(n); \
        _mm_free(x);

#endif /* _SIMD_H_ */

