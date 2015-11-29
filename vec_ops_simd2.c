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

#include"_fmmv.h"

#include "simd.h"
#if (FMM_PRECISION==0)
   #include"simd2s.h"
#else
   #include"simd2d.h"
#endif


void VEC_mul_simd2(const int n, const V2_BASETYPE *d, const V2_BASETYPE *x, V2_BASETYPE *y)
{
	/* pointers d,x,y must be 16 byte aligned */
	int i;
	int n0 = (n>>1);
	int n00 = (n0>>2);
	int n01 = (n0>>2)<<2;
	V2_BASETYPE *yp;
	const V2_BASETYPE *xp, *dp;
	yp = y;
	xp = x;	
	dp = d;
	for (i=0; i<n00; i++) {
		V2_STORE(yp, V2_MUL(V2_LOAD(dp), V2_LOAD(xp)));
		V2_STORE(yp+2, V2_MUL(V2_LOAD(dp+2), V2_LOAD(xp+2)));
		V2_STORE(yp+4, V2_MUL(V2_LOAD(dp+4), V2_LOAD(xp+4)));
		V2_STORE(yp+6, V2_MUL(V2_LOAD(dp+6), V2_LOAD(xp+6)));
		yp += 8;
		xp += 8;
		dp += 8;
	}
	for (i=n01; i<n0; i++) {
		V2_STORE(yp, V2_MUL(V2_LOAD(dp), V2_LOAD(xp)));
		yp += 2;
		xp += 2;
		dp += 2;
	}	
	if (n&1) {
		*yp = *dp * *xp;
	}
	
}

void VEC_mul2_simd2(const int n, const V2_BASETYPE *d, const V2_BASETYPE *x, V2_BASETYPE *y)
{
	/* pointers x,y must be 16 byte aligned */
	int i;
	int n0 = (n>>2);
	int n1 = (n>>2)<<2;

	for (i=0; i<n0; i++) {
		V2_STORE(y, V2_MUL(V2_LOAD1(d), V2_LOAD(x)));
		V2_STORE(y+2, V2_MUL(V2_LOAD1(d+1), V2_LOAD(x+2)));
		V2_STORE(y+4, V2_MUL(V2_LOAD1(d+2), V2_LOAD(x+4)));
		V2_STORE(y+6, V2_MUL(V2_LOAD1(d+3), V2_LOAD(x+6)));
		y += 8;
		d += 4;
		x += 8;
	}
	for (i=n1; i<n; i++) {
		V2_STORE(y, V2_MUL(V2_LOAD1(d), V2_LOAD(x)));
		y += 2;
		d++;
		x += 2;
	}
}

void VEC_scale_simd2(const int n, const V2_BASETYPE d0, const V2_BASETYPE *x, V2_BASETYPE *y)
{
	/* pointers d,x,y must be 16 byte aligned */
	int i;
	V2_TYPE d = V2_SET1(d0);
	int n1 = (n>>1)<<1;
	for (i=0; i<n1; i+=2) {
		V2_STORE(y+i, V2_MUL(d, V2_LOAD(x+i)));
	}
	if (n&1) {
		y[n1] = d0*x[n1];
	}
}


void VEC_add_simd2(const int n, const V2_BASETYPE *d, const V2_BASETYPE *x, V2_BASETYPE *y)
{
	/* pointers d,x,y must be 16 byte aligned */
	int i;
	int n0 = (n>>1);
	int n00 = (n0>>2);
	int n01 = (n0>>2)<<2;
	V2_BASETYPE *yp;
        const V2_BASETYPE *xp, *dp;
	yp = y;
	xp = x;	
	dp = d;
	for (i=0; i<n00; i++) {
		V2_STORE(yp, V2_ADD(V2_LOAD(dp), V2_LOAD(xp)));
		V2_STORE(yp+2, V2_ADD(V2_LOAD(dp+2), V2_LOAD(xp+2)));
		V2_STORE(yp+4, V2_ADD(V2_LOAD(dp+4), V2_LOAD(xp+4)));
		V2_STORE(yp+6, V2_ADD(V2_LOAD(dp+6), V2_LOAD(xp+6)));
		yp += 8;
		xp += 8;
		dp += 8;
	}
	for (i=n01; i<n0; i++) {
		V2_STORE(yp, V2_ADD(V2_LOAD(dp), V2_LOAD(xp)));
		yp += 2;
		xp += 2;
		dp += 2;
	}	
	if (n&1) {
		*yp = *dp + *xp;
	}
}

#include <string.h> /* memcpy */

V2_BASETYPE* VEC_add2_simd2(FmmvHandle *FMMV, const int n, const V2_BASETYPE *x, V2_BASETYPE *y)
{
	if (y) {
		int i;
		int n0 = (n>>1);
		int n00 = (n0>>2);
		int n01 = (n0>>2)<<2;
		V2_BASETYPE *yp;
	        const V2_BASETYPE *xp;
		yp = y;
		xp = x;	
		for (i=0; i<n00; i++) {
			V2_STORE(yp, V2_ADD(V2_LOAD(yp), V2_LOAD(xp)));
			V2_STORE(yp+2, V2_ADD(V2_LOAD(yp+2), V2_LOAD(xp+2)));
			V2_STORE(yp+4, V2_ADD(V2_LOAD(yp+4), V2_LOAD(xp+4)));
			V2_STORE(yp+6, V2_ADD(V2_LOAD(yp+6), V2_LOAD(xp+6)));
			yp += 8;
			xp += 8;
		}
		for (i=n01; i<n0; i++) {
			V2_STORE(yp, V2_ADD(V2_LOAD(yp), V2_LOAD(xp)));
			yp += 2;
			xp += 2;
		}	
		if (n&1) {
			*yp += *xp;
		}
		
	}
	else {
		y = (V2_BASETYPE*) FMMV_MALLOC(FMMV, n*sizeof(V2_BASETYPE));
		memcpy((void*) y, (void*) x, n*sizeof(V2_BASETYPE));
	}
	return y;
}	

void VEC_sub_simd2(const int n, const V2_BASETYPE *d, const V2_BASETYPE *x, V2_BASETYPE *y)
{
	/* pointers d,x,y must be 16 byte aligned */
	int i;
	int n0 = (n>>1);
	int n00 = (n0>>2);
	int n01 = (n0>>2)<<2;
	V2_BASETYPE *yp;
	const V2_BASETYPE  *xp, *dp;
	yp = y;
	xp = x;	
	dp = d;
	for (i=0; i<n00; i++) {
		V2_STORE(yp, V2_SUB(V2_LOAD(dp), V2_LOAD(xp)));
		V2_STORE(yp+2, V2_SUB(V2_LOAD(dp+2), V2_LOAD(xp+2)));
		V2_STORE(yp+4, V2_SUB(V2_LOAD(dp+4), V2_LOAD(xp+4)));
		V2_STORE(yp+6, V2_SUB(V2_LOAD(dp+6), V2_LOAD(xp+6)));
		yp += 8;
		xp += 8;
		dp += 8;
	}
	for (i=n01; i<n0; i++) {
		V2_STORE(yp, V2_SUB(V2_LOAD(dp), V2_LOAD(xp)));
		yp += 2;
		xp += 2;
		dp += 2;
	}	
	if (n&1) {
		*yp = *dp - *xp;
		xp++;
		yp++;
		dp++;
	}
}

void VEC_mul_c_rrii_simd2(const int n, const V2_BASETYPE *d, const V2_BASETYPE *x, V2_BASETYPE *y)
{
	/* Note: n must be divisible by 2; pointers d,x,y must be 16 byte aligned */
	int i;
	const V2_BASETYPE *di=d+n;
	const V2_BASETYPE *xi=x+n;
	V2_BASETYPE *yi=y+n;
	V2_TYPE dir, dii, xir, xii;

	for (i=0; i<n; i+=2) {
		dir = V2_LOAD(d+i);
		dii = V2_LOAD(di+i);
		xir = V2_LOAD(x+i);
		xii = V2_LOAD(xi+i);
		V2_STORE(y+i,V2_FMS(dir, xir, V2_MUL(dii,xii)));
		V2_STORE(yi+i,V2_FMA(dii, xir, V2_MUL(dir,xii)));
	}
}


void VEC_mul_cj_rrii_simd2(const int n, const V2_BASETYPE *d, const V2_BASETYPE *x, V2_BASETYPE *y)
{
	/* Note: n must be divisible by 2; pointers d,x,y must be 16 byte aligned */
	int i;
	const V2_BASETYPE *di=d+n;
	const V2_BASETYPE *xi=x+n;
	V2_BASETYPE *yi=y+n;
	V2_TYPE dir, dii, xir, xii;

	for (i=0; i<n; i+=2) {
		dir = V2_LOAD(d+i);
		dii = V2_LOAD(di+i);
		xir = V2_LOAD(x+i);
		xii = V2_LOAD(xi+i);
		V2_STORE(y+i,V2_FMA(dir, xir, V2_MUL(dii,xii)));
		V2_STORE(yi+i,V2_FMS(dir, xii, V2_MUL(dii,xir)));
	}
}	

void VEC_mul_ccj_rrii_simd2(const int n, const V2_BASETYPE *d, const V2_BASETYPE *x, V2_BASETYPE *y, V2_BASETYPE *z)
{
	/* Note: n must be divisible by 2; pointers d,x,y must be 16 byte aligned */
	int i;
	const V2_BASETYPE *di=d+n;
	const V2_BASETYPE *xi=x+n;
	V2_BASETYPE *yi=y+n;
	V2_BASETYPE *zi=z+n;
	V2_TYPE dir, dii, xir, xii;
	V2_TYPE a0, b0, c0, d0;

	for (i=0; i<n; i+=2) {
		
		dir = V2_LOAD(d+i);
		dii = V2_LOAD(di+i);
		xir = V2_LOAD(x+i);
		xii = V2_LOAD(xi+i);

		a0 = V2_MUL(dir,xir);
		b0 = V2_MUL(dii,xii);
		c0 = V2_MUL(dir,xii);
		d0 = V2_MUL(dii,xir);
		
		V2_STORE(y+i, V2_SUB(a0, b0));
		V2_STORE(yi+i, V2_ADD(c0, d0));
		V2_STORE(z+i, V2_ADD(a0, b0));
		V2_STORE(zi+i, V2_SUB(c0, d0));
	}
}	
	

void gemv_simd2(int m, int n, V2_BASETYPE *A, int lda, V2_BASETYPE *x, V2_BASETYPE *y)
{	
	int i, j;
	V2_BASETYPE *col;
	V2_TYPE temp;

	if (n>0) {
		for (i=0; i<m; i++) {
			col = A+lda*i;
			temp = V2_MUL(V2_LOAD1(col), V2_LOAD(x));
			for (j=1; j<n; j++) {
				temp = V2_FMA(V2_LOAD1(col+j), V2_LOAD(x+(j<<1)), temp);
			}
			V2_STORE(y+(i<<1), temp);
		}
	}
}

void gemv_trans_simd2(int m, int n, V2_BASETYPE *A, int lda, V2_BASETYPE *x, V2_BASETYPE *y)
{
	int i, j;
	V2_BASETYPE *col;
	V2_TYPE temp;

	temp = V2_ZERO;
	for (j=0; j<n; j++) {
		V2_STORE(y+(j<<1), temp);
	}		
	for (i=0; i<m; i++) {
		temp = V2_LOAD(x+(i<<1));
		col = A+lda*i;
		for (j=0; j<n; j++) {
			V2_STORE(y+(j<<1), V2_FMA(temp, V2_LOAD1(col+j), V2_LOAD(y+(j<<1))));
		}
	}	
}	

void tpmv_upper_simd2(int n, V2_BASETYPE *A, V2_BASETYPE *x)
{	
	int i, j, kk, k;
	V2_TYPE temp;

	kk = 0;
	for (j=0; j<n; j++) {
		temp = V2_MUL(V2_LOAD(x+(j<<1)), V2_LOAD1(A+kk));
		k = kk + 1;
		for (i=j+1; i<n; i++) {
			temp = V2_FMA(V2_LOAD(x+(i<<1)), V2_LOAD1(A+k), temp);
			k += 1;
		}	
		V2_STORE(x+(j<<1), temp);	
		kk += (n - j);
	}
}	

void tpmv_lower_simd2(int n, V2_BASETYPE *A, V2_BASETYPE *x)
{	
	int i, j, kk, k;
	V2_TYPE temp;

	kk = (n*(n+1))/2-1;
	for (j=n-1; j>=0; j--) {
		temp = V2_MUL(V2_LOAD(x+(j<<1)), V2_LOAD1(A+kk));
		k = kk - 1;	
		for (i=j-1; i>=0; i--) {
			temp = V2_FMA(V2_LOAD(x+(i<<1)), V2_LOAD1(A+k), temp);
			k -= 1;
		}	
		V2_STORE(x+(j<<1), temp);	
		kk -= (j + 1) ;
	}
}	


V2_BASETYPE *vec2_add2_simd2(FmmvHandle *FMMV, const int n, const V2_BASETYPE *x, V2_BASETYPE *y)
{
	int i, ii;
	
	if (!y) {
		y = (V2_BASETYPE*) FMMV_MALLOC(FMMV, n*sizeof(V2_BASETYPE));
		for (ii=0; ii<(n>>1); ii++) {
			i = ii<<2;
			V2_STORE(y+(ii<<1), V2_HORIZADD( V2_LOAD(x+i), V2_LOAD(x+i+2)));
						
		}
		for (i = ii<<1; i<n; i++) {
			ii = i<<1;
			y[i] =  x[ii] + x[ii+1];
		}	

	}
	else {
		for (ii=0; ii<(n>>1); ii++) {
			i = ii<<2;
			V2_STORE(y+(ii<<1) ,V2_ADD(V2_LOAD(y+(ii<<1)), 
				V2_HORIZADD( V2_LOAD(x+i), V2_LOAD(x+i+2))));

		}	
		for (i = ii<<1; i<n; i++) {
			ii = i<<1;
			y[i] +=  x[ii] + x[ii+1];
		}	
	}

	return y;
}	
			
void vec2_copy_simd2(const int n, const V2_BASETYPE *x1, const V2_BASETYPE *x2, V2_BASETYPE *y)
{
	int i;
	for (i=0; i<n; i++) {
			V2_STORE(y+(i<<1), V2_SET(x1[i], x2[i]));
	}
}	

void VEC_2add2_simd2(FmmvHandle *FMMV, int n, const V2_BASETYPE *x, V2_BASETYPE **y1, V2_BASETYPE **y2)
{
	int i, ii;	
	if (y2) {
		if (!*y2) {
			*y2 = (V2_BASETYPE*) FMMV_MALLOC(FMMV, n*sizeof(V2_BASETYPE));
			for (i=0, ii=0; i<n; i++, ii+=2) {
				(*y2)[i] = x[ii];
			}
		}
		else {
			for (i=0, ii=0; i<n; i++, ii+=2) {
				(*y2)[i] += x[ii];
			}
		}	
	}
	if (y1) {
		if (!*y1) {
			*y1 = (V2_BASETYPE*) FMMV_MALLOC(FMMV, n*sizeof(V2_BASETYPE));
			for (i=0, ii=1; i<n; i++, ii+=2) {
				(*y1)[i] = x[ii];
			}
		}
		else {
			for (i=0, ii=1; i<n; i++, ii+=2) {
				(*y1)[i] += x[ii];
			}
		}	
	}
}	

void P_2riri2rrii_simd2(int p, V2_BASETYPE *x1, V2_BASETYPE *x2, V2_BASETYPE *y)
{
	int n,m,i,ii;

	ii=0;
	for (n=0; n<=p; n++) {
		i = n*(n+1);
		for (m=0; m<=n; m++) {
			V2_STORE(y+ii, V2_SET(x1[i], x2[i]));
			i += 2;
			ii += 2;
		}	
		i = n*(n+1)+3;
		for (m=1; m<=n; m++) {
			V2_STORE(y+ii, V2_SET(x1[i], x2[i]));
			i += 2;
			ii += 2;
		}
	}	
}

void P_rrii2riri2_simd2(FmmvHandle *FMMV, int p, const V2_BASETYPE *x, V2_BASETYPE **y1, V2_BASETYPE **y2, V2_BASETYPE *mem)
{
        int n = (p+1)*(p+2);
	int i, ii, j, m, mi;
	int n1;
	
        n1 = n;
        if (n%2) n1 = n + 2 -(n%2);

	ii = 0;
	for (i=0; i<=p; i++) {
		m = 2*i*i;
		mi = m+2*i;
		mem[ii] = x[m];
		mem[n1+ii] = x[m+1];
		ii++;
		mem[ii] = 0.0;
		mem[n1+ii] = 0.0;
		ii++;
		for (j=2; j<=2*i; j+=2) {
			mem[ii] = x[m+j];
			mem[n1+ii] = x[m+j+1];
			ii++;
			mem[ii] = x[mi+j];
			mem[n1+ii] = x[mi+j+1];
			ii++;
		}	
	}
	
        if (y1) {
                *y1 = VEC_add2_simd2(FMMV, n, mem+1*n1, *y1);
        }
        if (y2) {
                *y2 = VEC_add2_simd2(FMMV, n, mem+0*n1, *y2);
        }
}

void P_X_2rrii2riri_simd2(int s, V2_BASETYPE *x1, V2_BASETYPE *x2, V2_BASETYPE *y)
{
	int i;
	int s0 = (s>>2);
	int s1 = (s>>2)<<2;
	V2_BASETYPE *yp;
	V2_BASETYPE *x1p, *x2p;
	
	yp = y;
	x1p = x1;
	x2p = x2;
	for (i=0; i<s0; i++) {
		V2_STORE(yp   , V2_SET(*x1p, *x2p));
		V2_STORE(yp+4 , V2_SET(*(x1p+1), *(x2p+1)));
		V2_STORE(yp+8 , V2_SET(*(x1p+2), *(x2p+2)));
		V2_STORE(yp+12, V2_SET(*(x1p+3), *(x2p+3)));
		yp += 16;
		x1p+=4;
		x2p+=4;
	}
	for (i=s1; i<s; i++) {
		V2_STORE(yp, V2_SET(*x1p, *x2p));
		yp += 4;
		x1p++;
		x2p++;
	}
	yp = y+2;
	x1p = x1+s;
	x2p = x2+s;
	for (i=0; i<s0; i++) {
		V2_STORE(yp   , V2_SET(*x1p, *x2p));
		V2_STORE(yp+4 , V2_SET(*(x1p+1), *(x2p+1)));
		V2_STORE(yp+8 , V2_SET(*(x1p+2), *(x2p+2)));
		V2_STORE(yp+12, V2_SET(*(x1p+3), *(x2p+3)));
		yp += 16;
		x1p+=4;
		x2p+=4;
	}
	for (i=s1; i<s; i++) {
		V2_STORE(yp, V2_SET(*x1p, *x2p));
		yp += 4;
		x1p++;
		x2p++;
	}
}	


void P_X_riri2rrii2_simd2(int s, const V2_BASETYPE *x, V2_BASETYPE **y1,  V2_BASETYPE **y2)
{
        int i;
	
	int s0 = (s>>2);
	int s1 = (s>>2)<<2;
	const V2_BASETYPE *xp, *x1p;
	V2_BASETYPE *yp, *y1p;
	if (y2) {
		xp = x;
		x1p = x+2;
		yp = *y2;
		y1p = (*y2)+s;
		for (i=0; i<s0; i++) {
			*yp += *xp;
			*(yp+1) += *(xp+4);
			*(yp+2) += *(xp+8);
			*(yp+3) += *(xp+12);
			*y1p += *x1p;
			*(y1p+1) += *(x1p+4);
			*(y1p+2) += *(x1p+8);
			*(y1p+3) += *(x1p+12);
			xp+=16;
			x1p+=16;
			yp+=4;
			y1p+=4;
		}
		
		for (i=s1; i<s; i++) {
			*yp += *xp;
			*y1p += *x1p;
			xp+=4;
			x1p+=4;
			yp++;
			y1p++;
		}
		
	}
	if (y1) {
		xp = x+1;
		x1p = x+3;
		yp = *y1;
		y1p = (*y1)+s;
		for (i=0; i<s0; i++) {
			*yp += *xp;
			*(yp+1) += *(xp+4);
			*(yp+2) += *(xp+8);
			*(yp+3) += *(xp+12);
			*y1p += *x1p;
			*(y1p+1) += *(x1p+4);
			*(y1p+2) += *(x1p+8);
			*(y1p+3) += *(x1p+12);
			xp+=16;
			x1p+=16;
			yp+=4;
			y1p+=4;
		}
		
		for (i=s1; i<s; i++) {
			*yp += *xp;
			*y1p += *x1p;
			xp+=4;
			x1p+=4;
			yp++;
			y1p++;
		}
		
	}
}

void perm_simd2(const int n, const int *vec, const V2_BASETYPE *x, V2_BASETYPE *y)
{
	int i;
	for (i=0; i<n; i++) {
		V2_SETITEM(y, i<<1, V2_GETITEM(x, vec[i]<<1));
	}	
}	

void perm_inv_simd2(const int n, const int *vec, const V2_BASETYPE *x, V2_BASETYPE *y)
{
	int i;
	for (i=0; i<n; i++) {
		V2_SETITEM(y, vec[i]<<1, V2_GETITEM(x, i<<1));
	}	
}	

void neg_simd2(const int n, const int *vec, V2_BASETYPE *x)
{
        int i;
        for (i=0; i<n; i++) {
                V2_SETITEM(x, vec[i]<<1, V2_NEG(V2_GETITEM(x, vec[i]<<1)));
        }
}
void VEC_conj(const int n, _FLOAT_ *y)
{
	int i;
	for (i=1; i<2*n+1; i+=2) {
		y[i] = -y[i];
	}
}	

void VEC_conj2_simd2(const int n, V2_BASETYPE *y)
{
	int i;
	for (i=2; i<4*n+2; i+=4) {
		V2_SETITEM(y, i, V2_NEG(V2_GETITEM(y, i)));
	}
}	


#include<math.h> // ldexpf

void scale_X_simd2(int p, V2_BASETYPE *xx, int level)
{ 
    V2_TYPE scale_p = V2_SET1(ldexp(+1.0, level)); 
    V2_TYPE scale_m = V2_SET1(ldexp(-1.0, level)); 
    int n, j;
    int jj = 0;

    for (n=0; n<=p; n++) {
       if(n&1) { /* n odd */
           for(j=0;j<2*n+1;j++) { // NOTE: Check this for reduced_scheme !!!
               // xx[jj] = -ldexp(xx[jj], level);
               V2_SETITEM(xx, jj, V2_MUL(V2_GETITEM(xx, jj), scale_m));
               jj+=2;
           }
       }
       else {
           for(j=0;j<2*n+1;j++) {
               //xx[jj] = +ldexp(xx[jj], level);
               V2_SETITEM(xx, jj, V2_MUL(V2_GETITEM(xx, jj), scale_p));
               jj+=2;
           }
       }
    }
}
