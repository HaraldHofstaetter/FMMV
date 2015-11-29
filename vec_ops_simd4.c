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

#include"_fmmv.h"

#include"simd.h"
#if (FMM_PRECISION==0)
   #include"simd4s.h"
#else
   #include"simd4d.h"
#endif


void VEC_mul_simd4(const int n, const V4_BASETYPE *d, const V4_BASETYPE *x, V4_BASETYPE *y)
{
	/* pointers d,x,y must be 16 byte aligned */
	int i;
	int n0 = (n>>2);
	int n1 = (n>>2)<<2;
	int n00 = (n0>>2);
	int n01 = (n0>>2)<<2;
	V4_BASETYPE *yp;
	const V4_BASETYPE *xp, *dp;
	yp = y;
	xp = x;	
	dp = d;
	for (i=0; i<n00; i++) {
		V4_STORE(yp, V4_MUL(V4_LOAD(dp), V4_LOAD(xp)));
		V4_STORE(yp+4, V4_MUL(V4_LOAD(dp+4), V4_LOAD(xp+4)));
		V4_STORE(yp+8, V4_MUL(V4_LOAD(dp+8), V4_LOAD(xp+8)));
		V4_STORE(yp+12, V4_MUL(V4_LOAD(dp+12), V4_LOAD(xp+12)));
		yp += 16;
		xp += 16;
		dp += 16;
	}
	for (i=n01; i<n0; i++) {
		V4_STORE(yp, V4_MUL(V4_LOAD(dp), V4_LOAD(xp)));
		yp += 4;
		xp += 4;
		dp += 4;
	}	
	for ( i=n1; i<n; i++) {
		*yp = *dp * *xp;
		xp++;
		yp++;
		dp++;
	}
	
}

void VEC_mul4_simd4(const int n, const V4_BASETYPE *d, const V4_BASETYPE *x, V4_BASETYPE *y)
{
	/* pointers x,y must be 16 byte aligned */
	int i;
	int n0 = (n>>2);
	int n1 = (n>>2)<<2;

	for (i=0; i<n0; i++) {
		V4_STORE(y, V4_MUL(V4_LOAD1(d), V4_LOAD(x)));
		V4_STORE(y+4, V4_MUL(V4_LOAD1(d+1), V4_LOAD(x+4)));
		V4_STORE(y+8, V4_MUL(V4_LOAD1(d+2), V4_LOAD(x+8)));
		V4_STORE(y+12, V4_MUL(V4_LOAD1(d+3), V4_LOAD(x+12)));
		y += 16;
		d += 4;
		x += 16;
	}
	for (i=n1; i<n; i++) {
		V4_STORE(y, V4_MUL(V4_LOAD1(d), V4_LOAD(x)));
		y += 4;
		d++;
		x += 4;
	}
}

void VEC_scale_simd4(const int n, const V4_BASETYPE d0, const V4_BASETYPE *x, V4_BASETYPE *y)
{
	/* pointers d,x,y must be 16 byte aligned */
	int i;
	V4_TYPE d = V4_SET1(d0);
	int n1 = (n>>2)<<2;
	for (i=0; i<n1; i+=4) {
		V4_STORE(y+i, V4_MUL(d, V4_LOAD(x+i)));
	}
	for ( i=n1; i<n; i++) {
		y[i] = d0*x[i];
	}
}


void VEC_add_simd4(const int n, const V4_BASETYPE *d, const V4_BASETYPE *x, V4_BASETYPE *y)
{
	/* pointers d,x,y must be 16 byte aligned */
	int i;
	int n0 = (n>>2);
	int n1 = (n>>2)<<2;
	int n00 = (n0>>2);
	int n01 = (n0>>2)<<2;
	V4_BASETYPE *yp;
        const V4_BASETYPE *xp, *dp;
	yp = y;
	xp = x;	
	dp = d;
	for (i=0; i<n00; i++) {
		V4_STORE(yp, V4_ADD(V4_LOAD(dp), V4_LOAD(xp)));
		V4_STORE(yp+4, V4_ADD(V4_LOAD(dp+4), V4_LOAD(xp+4)));
		V4_STORE(yp+8, V4_ADD(V4_LOAD(dp+8), V4_LOAD(xp+8)));
		V4_STORE(yp+12, V4_ADD(V4_LOAD(dp+12), V4_LOAD(xp+12)));
		yp += 16;
		xp += 16;
		dp += 16;
	}
	for (i=n01; i<n0; i++) {
		V4_STORE(yp, V4_ADD(V4_LOAD(dp), V4_LOAD(xp)));
		yp += 4;
		xp += 4;
		dp += 4;
	}	
	for ( i=n1; i<n; i++) {
		*yp = *dp + *xp;
		xp++;
		yp++;
		dp++;
	}
}

#include <string.h> /* memcpy */

V4_BASETYPE* VEC_add2_simd4(FmmvHandle *FMMV, const int n, const V4_BASETYPE *x, V4_BASETYPE *y)
{
	if (y) {
		int i;
		int n0 = (n>>2);
		int n1 = (n>>2)<<2;
		int n00 = (n0>>2);
		int n01 = (n0>>2)<<2;
		V4_BASETYPE *yp;
	        const V4_BASETYPE *xp;
		yp = y;
		xp = x;	
		for (i=0; i<n00; i++) {
			V4_STORE(yp, V4_ADD(V4_LOAD(yp), V4_LOAD(xp)));
			V4_STORE(yp+4, V4_ADD(V4_LOAD(yp+4), V4_LOAD(xp+4)));
			V4_STORE(yp+8, V4_ADD(V4_LOAD(yp+8), V4_LOAD(xp+8)));
			V4_STORE(yp+12, V4_ADD(V4_LOAD(yp+12), V4_LOAD(xp+12)));
			yp += 16;
			xp += 16;
		}
		for (i=n01; i<n0; i++) {
			V4_STORE(yp, V4_ADD(V4_LOAD(yp), V4_LOAD(xp)));
			yp += 4;
			xp += 4;
		}	
		for ( i=n1; i<n; i++) {
			*yp += *xp;
			xp++;
			yp++;
		}
		
	}
	else {
		y = (V4_BASETYPE*) FMMV_MALLOC(FMMV, n*sizeof(V4_BASETYPE));
		memcpy((void*) y, (void*) x, n*sizeof(V4_BASETYPE));
	}
	return y;
}	
/*
V4_BASETYPE* VEC_add2_simd4(FmmvHandle *FMMV, const int n, const V4_BASETYPE *x, V4_BASETYPE *y)
{
	if (y) {
		int i;
		int n1 = (n>>2)<<2;
		for (i=0; i<n1; i+=4) {
			V4_STORE(y+i, V4_ADD(V4_LOAD(y+i), V4_LOAD(x+i)));
		}
		for ( i=n1; i<n; i++) {
			y[i] += +x[i];
		}
	}
	else {
		y = (V4_BASETYPE*) FMMV_MALLOC(FMMV, n*sizeof(V4_BASETYPE));
		memcpy((void*) y, (void*) x, n*sizeof(V4_BASETYPE));
	}
	return y;
}
*/


void VEC_sub_simd4(const int n, const V4_BASETYPE *d, const V4_BASETYPE *x, V4_BASETYPE *y)
{
	/* pointers d,x,y must be 16 byte aligned */
	int i;
	int n0 = (n>>2);
	int n1 = (n>>2)<<2;
	int n00 = (n0>>2);
	int n01 = (n0>>2)<<2;
	V4_BASETYPE *yp;
	const V4_BASETYPE  *xp, *dp;
	yp = y;
	xp = x;	
	dp = d;
	for (i=0; i<n00; i++) {
		V4_STORE(yp, V4_SUB(V4_LOAD(dp), V4_LOAD(xp)));
		V4_STORE(yp+4, V4_SUB(V4_LOAD(dp+4), V4_LOAD(xp+4)));
		V4_STORE(yp+8, V4_SUB(V4_LOAD(dp+8), V4_LOAD(xp+8)));
		V4_STORE(yp+12, V4_SUB(V4_LOAD(dp+12), V4_LOAD(xp+12)));
		yp += 16;
		xp += 16;
		dp += 16;
	}
	for (i=n01; i<n0; i++) {
		V4_STORE(yp, V4_SUB(V4_LOAD(dp), V4_LOAD(xp)));
		yp += 4;
		xp += 4;
		dp += 4;
	}	
	for ( i=n1; i<n; i++) {
		*yp = *dp - *xp;
		xp++;
		yp++;
		dp++;
	}
}

#if 0
void VEC_mul_c_rrii_simd4(const int n, const V4_BASETYPE *d, const V4_BASETYPE *x, V4_BASETYPE *y)
{
	/* Note: n must be divisible by 4; pointers d,x,y must be 16 byte aligned */
	int i;
	const V4_BASETYPE *dr=d;
	const V4_BASETYPE *xr=x;
	const V4_BASETYPE *di=d+n;
	const V4_BASETYPE *xi=x+n;
	V4_BASETYPE *yr=y;
	V4_BASETYPE *yi=y+n;
	V4_TYPE dr0, di0, xr0, xi0;
	/*
	V4_TYPE dr1, di1, xr1, xi1;
	V4_TYPE dr2, di2, xr2, xi2;
	V4_TYPE dr3, di3, xr3, xi3;
	int n0 = (n>>3);
	int n1 = (n>>3)<<2;
	for (i=0; i<n0; i++) {
		dr0 = V4_LOAD(dr);
		dr1 = V4_LOAD(dr+4);
		dr2 = V4_LOAD(dr+8);
		dr3 = V4_LOAD(dr+12);
		di0 = V4_LOAD(di);
		di1 = V4_LOAD(di+4);
		di2 = V4_LOAD(di+8);
		di3 = V4_LOAD(di+12);
		xr0 = V4_LOAD(xr);
		xr1 = V4_LOAD(xr+4);
		xr2 = V4_LOAD(xr+8);
		xr3 = V4_LOAD(xr+12);
		xi0 = V4_LOAD(xi);
		xi1 = V4_LOAD(xi+4);
		xi2 = V4_LOAD(xi+8);
		xi3 = V4_LOAD(xi+12);
		V4_STORE(yr,V4_SUB(V4_MUL(dr0,xr0),V4_MUL(di0,xi0)));
		V4_STORE(yr+4,V4_SUB(V4_MUL(dr1,xr1),V4_MUL(di1,xi1)));
		V4_STORE(yr+8,V4_SUB(V4_MUL(dr2,xr2),V4_MUL(di2,xi2)));
		V4_STORE(yr+12,V4_SUB(V4_MUL(dr3,xr3),V4_MUL(di3,xi3)));
		V4_STORE(yi,V4_ADD(V4_MUL(di0,xr0),V4_MUL(dr0,xi0)));
		V4_STORE(yi+4,V4_ADD(V4_MUL(di1,xr1),V4_MUL(dr1,xi1)));
		V4_STORE(yi+8,V4_ADD(V4_MUL(di2,xr2),V4_MUL(dr2,xi2)));
		V4_STORE(yi+12,V4_ADD(V4_MUL(di3,xr3),V4_MUL(dr3,xi3)));
		dr += 8;
		di += 8;
		xr += 8;
		xi += 8;
		yr += 8;
		yi += 8;
	}
	for (i=n1; i<(n>>2); i++) {
	*/
	for (i=0; i<(n>>2); i++) {
		dr0 = V4_LOAD(dr);
		di0 = V4_LOAD(di);
		xr0 = V4_LOAD(xr);
		xi0 = V4_LOAD(xi);
		V4_STORE(yr,V4_SUB(V4_MUL(dr0,xr0),V4_MUL(di0,xi0)));
		V4_STORE(yi,V4_ADD(V4_MUL(di0,xr0),V4_MUL(dr0,xi0)));
		dr += 4;
		di += 4;
		xr += 4;
		xi += 4;
		yr += 4;
		yi += 4;
	}
}
#endif

void VEC_mul_c_rrii_simd4(const int n, const V4_BASETYPE *d, const V4_BASETYPE *x, V4_BASETYPE *y)
{
	/* Note: n must be divisible by 4; pointers d,x,y must be 16 byte aligned */
	int i;
	const V4_BASETYPE *di=d+n;
	const V4_BASETYPE *xi=x+n;
	V4_BASETYPE *yi=y+n;
	V4_TYPE dir, dii, xir, xii;

	for (i=0; i<n; i+=4) {
		dir = V4_LOAD(d+i);
		dii = V4_LOAD(di+i);
		xir = V4_LOAD(x+i);
		xii = V4_LOAD(xi+i);
		V4_STORE(y+i,V4_FMS(dir, xir, V4_MUL(dii,xii)));
		V4_STORE(yi+i,V4_FMA(dii, xir, V4_MUL(dir,xii)));
	}
}


void VEC_mul_cj_rrii_simd4(const int n, const V4_BASETYPE *d, const V4_BASETYPE *x, V4_BASETYPE *y)
{
	/* Note: n must be divisible by 4; pointers d,x,y must be 16 byte aligned */
	int i;
	const V4_BASETYPE *di=d+n;
	const V4_BASETYPE *xi=x+n;
	V4_BASETYPE *yi=y+n;
	V4_TYPE dir, dii, xir, xii;

	for (i=0; i<n; i+=4) {
		dir = V4_LOAD(d+i);
		dii = V4_LOAD(di+i);
		xir = V4_LOAD(x+i);
		xii = V4_LOAD(xi+i);
		V4_STORE(y+i,V4_FMA(dir, xir, V4_MUL(dii,xii)));
		V4_STORE(yi+i,V4_FMS(dir, xii, V4_MUL(dii,xir)));
	}
}	

void VEC_mul_ccj_rrii_simd4(const int n, const V4_BASETYPE *d, const V4_BASETYPE *x, V4_BASETYPE *y, V4_BASETYPE *z)
{
	/* Note: n must be divisible by 4; pointers d,x,y must be 16 byte aligned */
	int i;
	const V4_BASETYPE *di=d+n;
	const V4_BASETYPE *xi=x+n;
	V4_BASETYPE *yi=y+n;
	V4_BASETYPE *zi=z+n;
	V4_TYPE dir, dii, xir, xii;
	V4_TYPE a0, b0, c0, d0;

	for (i=0; i<n; i+=4) {
		
		dir = V4_LOAD(d+i);
		dii = V4_LOAD(di+i);
		xir = V4_LOAD(x+i);
		xii = V4_LOAD(xi+i);

		a0 = V4_MUL(dir,xir);
		b0 = V4_MUL(dii,xii);
		c0 = V4_MUL(dir,xii);
		d0 = V4_MUL(dii,xir);
		
		V4_STORE(y+i, V4_SUB(a0, b0));
		V4_STORE(yi+i, V4_ADD(c0, d0));
		V4_STORE(z+i, V4_ADD(a0, b0));
		V4_STORE(zi+i, V4_SUB(c0, d0));
	}
}	
	

void gemv_simd4(int m, int n, V4_BASETYPE *A, int lda, V4_BASETYPE *x, V4_BASETYPE *y)
{	
	int i, j;
	V4_BASETYPE *col;
	V4_TYPE temp;

	if (n>0) {
		for (i=0; i<m; i++) {
			col = A+lda*i;
			temp = V4_MUL(V4_LOAD1(col), V4_LOAD(x));
			for (j=1; j<n; j++) {
				temp = V4_FMA(V4_LOAD1(col+j), V4_LOAD(x+(j<<2)), temp);
			}
			V4_STORE(y+(i<<2), temp);
		}
	}
}

void gemv_trans_simd4(int m, int n, V4_BASETYPE *A, int lda, V4_BASETYPE *x, V4_BASETYPE *y)
{
	int i, j;
	V4_BASETYPE *col;
	V4_TYPE temp;

	temp = V4_ZERO;
	for (j=0; j<n; j++) {
		V4_STORE(y+(j<<2), temp);
	}		
	for (i=0; i<m; i++) {
		temp = V4_LOAD(x+(i<<2));
		col = A+lda*i;
		for (j=0; j<n; j++) {
			V4_STORE(y+(j<<2), V4_FMA(temp, V4_LOAD1(col+j), V4_LOAD(y+(j<<2))));
		}
	}	
}	

void tpmv_upper_simd4(int n, V4_BASETYPE *A, V4_BASETYPE *x)
{	
	int i, j, kk, k;
	V4_TYPE temp;

	kk = 0;
	for (j=0; j<n; j++) {
		temp = V4_MUL(V4_LOAD(x+(j<<2)), V4_LOAD1(A+kk));
		k = kk + 1;
		for (i=j+1; i<n; i++) {
			temp = V4_FMA(V4_LOAD(x+(i<<2)), V4_LOAD1(A+k), temp);
			k += 1;
		}	
		V4_STORE(x+(j<<2), temp);	
		kk += (n - j);
	}
}	

void tpmv_lower_simd4(int n, V4_BASETYPE *A, V4_BASETYPE *x)
{	
	int i, j, kk, k;
	V4_TYPE temp;

	kk = (n*(n+1))/2-1;
	for (j=n-1; j>=0; j--) {
		temp = V4_MUL(V4_LOAD(x+(j<<2)), V4_LOAD1(A+kk));
		k = kk - 1;	
		for (i=j-1; i>=0; i--) {
			temp = V4_FMA(V4_LOAD(x+(i<<2)), V4_LOAD1(A+k), temp);
			k -= 1;
		}	
		V4_STORE(x+(j<<2), temp);	
		kk -= (j + 1) ;
	}
}	


V4_BASETYPE *vec4_add2_simd4(FmmvHandle *FMMV, const int n, const V4_BASETYPE *x, V4_BASETYPE *y)
{
	int i, ii;
	
	if (!y) {
		y = (V4_BASETYPE*) FMMV_MALLOC(FMMV, n*sizeof(V4_BASETYPE));
		for (ii=0; ii<(n>>2); ii++) {
			i = ii<<4;
			V4_STORE(y+(ii<<2), V4_HORIZADD( V4_LOAD(x+i), V4_LOAD(x+i+4), V4_LOAD(x+i+8), V4_LOAD(x+i+12)));
						
		}
		for (i = ii<<2; i<n; i++) {
			ii = i<<2;
			y[i] =  x[ii] + x[ii+1] + x[ii+2] + x[ii+3];
		}	

	}
	else {
		for (ii=0; ii<(n>>2); ii++) {
			i = ii<<4;
			V4_STORE(y+(ii<<2) ,V4_ADD(V4_LOAD(y+(ii<<2)), 
				V4_HORIZADD( V4_LOAD(x+i), V4_LOAD(x+i+4), V4_LOAD(x+i+8), V4_LOAD(x+i+12))));

		}	
		for (i = ii<<2; i<n; i++) {
			ii = i<<2;
			y[i] +=  x[ii] + x[ii+1] + x[ii+2] + x[ii+3];
		}	
	}

	return y;
}	
			
void vec4_copy_simd4(const int n, const V4_BASETYPE *x1, const V4_BASETYPE *x2, const V4_BASETYPE *x3, const V4_BASETYPE *x4, V4_BASETYPE *y)
{
	int i;
	for (i=0; i<n; i++) {
			V4_STORE(y+(i<<2), V4_SET(x1[i], x2[i], x3[i], x4[i]));
	}
}	

void VEC_4add2_simd4(FmmvHandle *FMMV, int n, const V4_BASETYPE *x, V4_BASETYPE **y0, V4_BASETYPE **y1, V4_BASETYPE **y2, V4_BASETYPE **y3 /*, V4_BASETYPE *mem*/)
{
#ifdef V4_LITTLE_ENDIAN
    #define Y0 y3
    #define Y1 y2
    #define Y2 y1
    #define Y3 y0
#else
    #define Y0 y0
    #define Y1 y1
    #define Y2 y2
    #define Y3 y3
#endif	
	int i, ii;	
/*	
	int n1;
	n1 = n;
	if (n%4) n1 = n + 4 -(n%4);

	for (i=0; i<n; i++) {
		ii = (i<<2);
		mem[i] = x[ii];
		mem[n1+i] = x[ii+1];
		mem[2*n1+i] = x[ii+2];
		mem[3*n1+i] = x[ii+3];
	}	
	if (Y3) {
		*Y3 = VEC_add2_simd4(FMMV, n, mem+3*n1, *Y3);
	}	
	if (Y2) {
		*Y2 = VEC_add2_simd4(FMMV, n, mem+2*n1, *Y2);
	}	
	if (Y1) {
		*Y1 = VEC_add2_simd4(FMMV, n, mem+1*n1, *Y1);
	}	
	if (Y0) {
		*Y0 = VEC_add2_simd4(FMMV, n, mem+0*n1, *Y0);
	}
*/	
	if (Y0) {
		if (!*Y0) {
			*Y0 = (V4_BASETYPE*) FMMV_MALLOC(FMMV, n*sizeof(V4_BASETYPE));
			for (i=0, ii=0; i<n; i++, ii+=4) {
				(*Y0)[i] = x[ii];
			}
		}
		else {
			for (i=0, ii=0; i<n; i++, ii+=4) {
				(*Y0)[i] += x[ii];
			}
		}	
	}
	if (Y1) {
		if (!*Y1) {
			*Y1 = (V4_BASETYPE*) FMMV_MALLOC(FMMV, n*sizeof(V4_BASETYPE));
			for (i=0, ii=1; i<n; i++, ii+=4) {
				(*Y1)[i] = x[ii];
			}
		}
		else {
			for (i=0, ii=1; i<n; i++, ii+=4) {
				(*Y1)[i] += x[ii];
			}
		}	
	}
	if (Y2) {
		if (!*Y2) {
			*Y2 = (V4_BASETYPE*) FMMV_MALLOC(FMMV, n*sizeof(V4_BASETYPE));
			for (i=0, ii=2; i<n; i++, ii+=4) {
				(*Y2)[i] = x[ii];
			}
		}
		else {
			for (i=0, ii=2; i<n; i++, ii+=4) {
				(*Y2)[i] += x[ii];
			}
		}	
	}
	if (Y3) {
		if (!*Y3) {
			*Y3 = (V4_BASETYPE*) FMMV_MALLOC(FMMV, n*sizeof(V4_BASETYPE));
			for (i=0, ii=3; i<n; i++, ii+=4) {
				(*Y3)[i] = x[ii];
			}
		}
		else {
			for (i=0, ii=3; i<n; i++, ii+=4) {
				(*Y3)[i] += x[ii];
			}
		}	
	}
#undef Y0	
#undef Y1	
#undef Y2	
#undef Y3	
}	

void P_4riri2rrii_simd4(int p, V4_BASETYPE *x1, V4_BASETYPE *x2, V4_BASETYPE *x3, V4_BASETYPE *x4, V4_BASETYPE *y)
{
	int n,m,i,ii;

	ii=0;
	for (n=0; n<=p; n++) {
		i = n*(n+1);
		for (m=0; m<=n; m++) {
			V4_STORE(y+ii, V4_SET(x1[i], x2[i], x3[i], x4[i]));
			i += 2;
			ii += 4;
		}	
		i = n*(n+1)+3;
		for (m=1; m<=n; m++) {
			V4_STORE(y+ii, V4_SET(x1[i], x2[i], x3[i], x4[i]));
			i += 2;
			ii += 4;
		}
	}	
}

void P_rrii2riri4_simd4(FmmvHandle *FMMV, int p, const V4_BASETYPE *x, V4_BASETYPE **y0, V4_BASETYPE **y1, V4_BASETYPE **y2, V4_BASETYPE **y3, V4_BASETYPE *mem)
{
#ifdef V4_LITTLE_ENDIAN
    #define Y0 y3
    #define Y1 y2
    #define Y2 y1
    #define Y3 y0
#else
    #define Y0 y0
    #define Y1 y1
    #define Y2 y2
    #define Y3 y3
#endif	

        int n = (p+1)*(p+2);
	int i, ii, j, m, mi;
	int n1;
	
        n1 = n;
        if (n%4) n1 = n + 4 -(n%4);

	ii = 0;
	for (i=0; i<=p; i++) {
		m = 4*i*i;
		mi = m+4*i;
		mem[ii] = x[m];
		mem[n1+ii] = x[m+1];
		mem[2*n1+ii] = x[m+2];
		mem[3*n1+ii] = x[m+3];
		ii++;
		mem[ii] = 0.0;
		mem[n1+ii] = 0.0;
		mem[2*n1+ii] = 0.0;
		mem[3*n1+ii] = 0.0;
		ii++;
		for (j=4; j<=4*i; j+=4) {
			mem[ii] = x[m+j];
			mem[n1+ii] = x[m+j+1];
			mem[2*n1+ii] = x[m+j+2];
			mem[3*n1+ii] = x[m+j+3];
			ii++;
			mem[ii] = x[mi+j];
			mem[n1+ii] = x[mi+j+1];
			mem[2*n1+ii] = x[mi+j+2];
			mem[3*n1+ii] = x[mi+j+3];
			ii++;
		}	
	}
	
        if (Y3) {
                *Y3 = VEC_add2_simd4(FMMV, n, mem+3*n1, *Y3);
        }
        if (Y2) {
                *Y2 = VEC_add2_simd4(FMMV, n, mem+2*n1, *Y2);
        }
        if (Y1) {
                *Y1 = VEC_add2_simd4(FMMV, n, mem+1*n1, *Y1);
        }
        if (Y0) {
                *Y0 = VEC_add2_simd4(FMMV, n, mem+0*n1, *Y0);
        }
/*	
	
	if (Y0) {
		if (!*Y0) {
			*Y0 = (V4_BASETYPE*) FMMV_MALLOC(FMMV, n*sizeof(V4_BASETYPE));
			ii = 0;
			for (i=0; i<=p; i++) {
				m = 4*i*i;
				mi = m+4*i;
				(*Y0)[ii] = x[m];
				(*Y0)[ii+1] = 0.0;
				ii+=2;
				for (j=4; j<=4*i; j+=4) {
					(*Y0)[ii] = x[m+j];
					(*Y0)[ii+1] = x[mi+j];
					ii+=2;
				}	
			}
		}
		else {
			ii = 0;
			for (i=0; i<=p; i++) {
				m = 4*i*i;
				mi = m+4*i;
				(*Y0)[ii] += x[m];
				ii+=2;
				for (j=4; j<=4*i; j+=4) {
					(*Y0)[ii] += x[m+j];
					(*Y0)[ii+1] += x[mi+j];
					ii+=2;
				}	
			}
		}	
	}
	if (Y1) {
		if (!*Y1) {
			*Y1 = (V4_BASETYPE*) FMMV_MALLOC(FMMV, n*sizeof(V4_BASETYPE));
			ii = 0;
			for (i=0; i<=p; i++) {
				m = 4*i*i+1;
				mi = m+4*i;
				(*Y1)[ii] = x[m];
				(*Y1)[ii+1] = 0.0;
				ii+=2;
				for (j=4; j<=4*i; j+=4) {
					(*Y1)[ii] = x[m+j];
					(*Y1)[ii+1] = x[mi+j];
					ii+=2;
				}	
			}
		}
		else {
			ii = 0;
			for (i=0; i<=p; i++) {
				m = 4*i*i+1;
				mi = m+4*i;
				(*Y1)[ii] += x[m];
				ii+=2;
				for (j=4; j<=4*i; j+=4) {
					(*Y1)[ii] += x[m+j];
					(*Y1)[ii+1] += x[mi+j];
					ii+=2;
				}	
			}
		}	
	}
	if (Y2) {
		if (!*Y2) {
			*Y2 = (V4_BASETYPE*) FMMV_MALLOC(FMMV, n*sizeof(V4_BASETYPE));
			ii = 0;
			for (i=0; i<=p; i++) {
				m = 4*i*i+2;
				mi = m+4*i;
				(*Y2)[ii] = x[m];
				(*Y2)[ii+1] = 0.0;
				ii+=2;
				for (j=4; j<=4*i; j+=4) {
					(*Y2)[ii] = x[m+j];
					(*Y2)[ii+1] = x[mi+j];
					ii+=2;
				}	
			}
		}
		else {
			ii = 0;
			for (i=0; i<=p; i++) {
				m = 4*i*i+2;
				mi = m+4*i;
				(*Y2)[ii] += x[m];
				ii+=2;
				for (j=4; j<=4*i; j+=4) {
					(*Y2)[ii] += x[m+j];
					(*Y2)[ii+1] += x[mi+j];
					ii+=2;
				}	
			}
		}	
	}
	if (Y3) {
		if (!*Y3) {
			*Y3 = (V4_BASETYPE*) FMMV_MALLOC(FMMV, n*sizeof(V4_BASETYPE));
			ii = 0;
			for (i=0; i<=p; i++) {
				m = 4*i*i+3;
				mi = m+4*i;
				(*Y3)[ii] = x[m];
				(*Y3)[ii+1] = 0.0;
				ii+=2;
				for (j=4; j<=4*i; j+=4) {
					(*Y3)[ii] = x[m+j];
					(*Y3)[ii+1] = x[mi+j];
					ii+=2;
				}	
			}
		}
		else {
			ii = 0;
			for (i=0; i<=p; i++) {
				m = 4*i*i+3;
				mi = m+4*i;
				(*Y3)[ii] += x[m];
				ii+=2;
				for (j=4; j<=4*i; j+=4) {
					(*Y3)[ii] += x[m+j];
					(*Y3)[ii+1] += x[mi+j];
					ii+=2;
				}	
			}
		}	
	}
*/	
#undef Y0	
#undef Y1	
#undef Y2	
#undef Y3	
}


void P_X_4rrii2riri_simd4(int s, V4_BASETYPE *x1, V4_BASETYPE *x2, V4_BASETYPE *x3, V4_BASETYPE *x4, V4_BASETYPE *y)
{
	int i;
	int s0 = (s>>2);
	int s1 = (s>>2)<<2;
	V4_BASETYPE *yp;
	V4_BASETYPE *x1p, *x2p, *x3p, *x4p;
	
	yp = y;
	x1p = x1;
	x2p = x2;
	x3p = x3;
	x4p = x4;
	for (i=0; i<s0; i++) {
		V4_STORE(yp   , V4_SET(*x1p, *x2p, *x3p, *x4p));
		V4_STORE(yp+8 , V4_SET(*(x1p+1), *(x2p+1), *(x3p+1), *(x4p+1)));
		V4_STORE(yp+16, V4_SET(*(x1p+2), *(x2p+2), *(x3p+2), *(x4p+2)));
		V4_STORE(yp+24, V4_SET(*(x1p+3), *(x2p+3), *(x3p+3), *(x4p+3)));
		yp += 32;
		x1p+=4;
		x2p+=4;
		x3p+=4;
		x4p+=4;
	}
	for (i=s1; i<s; i++) {
		V4_STORE(yp, V4_SET(*x1p, *x2p, *x3p, *x4p));
		yp += 8;
		x1p++;
		x2p++;
		x3p++;
		x4p++;
	}
	yp = y+4;
	x1p = x1+s;
	x2p = x2+s;
	x3p = x3+s;
	x4p = x4+s;
	for (i=0; i<s0; i++) {
		V4_STORE(yp   , V4_SET(*x1p, *x2p, *x3p, *x4p));
		V4_STORE(yp+8 , V4_SET(*(x1p+1), *(x2p+1), *(x3p+1), *(x4p+1)));
		V4_STORE(yp+16, V4_SET(*(x1p+2), *(x2p+2), *(x3p+2), *(x4p+2)));
		V4_STORE(yp+24, V4_SET(*(x1p+3), *(x2p+3), *(x3p+3), *(x4p+3)));
		yp += 32;
		x1p+=4;
		x2p+=4;
		x3p+=4;
		x4p+=4;
	}
	for (i=s1; i<s; i++) {
		V4_STORE(yp, V4_SET(*x1p, *x2p, *x3p, *x4p));
		yp += 8;
		x1p++;
		x2p++;
		x3p++;
		x4p++;
	}
/*
	for (i=0, ii=0; i<s; i++, ii+=8)  y[ii]=x4[i];
	for (i=s, ii=4; i<2*s; i++, ii+=8)  y[ii]=x4[i];

	for (i=0, ii=1; i<s; i++, ii+=8)  y[ii]=x3[i];
	for (i=s, ii=5; i<2*s; i++, ii+=8)  y[ii]=x3[i];

	for (i=0, ii=2; i<s; i++, ii+=8)  y[ii]=x2[i];
	for (i=s, ii=6; i<2*s; i++, ii+=8)  y[ii]=x2[i];

	for (i=0, ii=3; i<s; i++, ii+=8)  y[ii]=x1[i];
	for (i=s, ii=7; i<2*s; i++, ii+=8)  y[ii]=x1[i];
*/	
}	


void P_X_riri2rrii4_simd4(int s, const V4_BASETYPE *x, V4_BASETYPE **y0, V4_BASETYPE **y1, V4_BASETYPE **y2, V4_BASETYPE **y3)
{
#ifdef V4_LITTLE_ENDIAN
    #define Y0 y3
    #define Y1 y2
    #define Y2 y1
    #define Y3 y0
#else
    #define Y0 y0
    #define Y1 y1
    #define Y2 y2
    #define Y3 y3
#endif	
        int i;
	
	int s0 = (s>>2);
	int s1 = (s>>2)<<2;
	const V4_BASETYPE *xp, *x1p;
	V4_BASETYPE *yp, *y1p;
	if (Y0) {
		/*if (!*Y0) {
			*Y0 = (V4_BASETYPE*) FMMV_MALLOC(FMMV, n*sizeof(V4_BASETYPE));
			for (i=0, ii=0; i<s; i++, ii+=8) {
				(*Y0)[i] = x[ii];
				(*Y0)[i+s] = x[ii+4];
			}
		}
		else {
			for (i=0, ii=0; i<s; i++, ii+=8) {
				(*Y0)[i] += x[ii];
				(*Y0)[i+s] += x[ii+4];
			}
		}*/	
		xp = x;
		x1p = x+4;
		yp = *Y0;
		y1p = (*Y0)+s;
		for (i=0; i<s0; i++) {
			*yp += *xp;
			*(yp+1) += *(xp+8);
			*(yp+2) += *(xp+16);
			*(yp+3) += *(xp+24);
			*y1p += *x1p;
			*(y1p+1) += *(x1p+8);
			*(y1p+2) += *(x1p+16);
			*(y1p+3) += *(x1p+24);
			xp+=32;
			x1p+=32;
			yp+=4;
			y1p+=4;
		}
		for (i=s1; i<s; i++) {
			*yp += *xp;
			*y1p += *x1p;
			xp+=8;
			x1p+=8;
			yp++;
			y1p++;
		}
		
	}
	if (Y1) {
		/*if (!*Y1) {
			*Y1 = (V4_BASETYPE*) FMMV_MALLOC(FMMV, n*sizeof(V4_BASETYPE));
			for (i=0, ii=1; i<s; i++, ii+=8) {
				(*Y1)[i] = x[ii];
				(*Y1)[i+s] = x[ii+4];
			}
		}
		else {
			for (i=0, ii=1; i<s; i++, ii+=8) {
				(*Y1)[i] += x[ii];
				(*Y1)[i+s] += x[ii+4];
			}
		}*/	
		xp = x+1;
		x1p = x+5;
		yp = *Y1;
		y1p = (*Y1)+s;
		for (i=0; i<s0; i++) {
			*yp += *xp;
			*(yp+1) += *(xp+8);
			*(yp+2) += *(xp+16);
			*(yp+3) += *(xp+24);
			*y1p += *x1p;
			*(y1p+1) += *(x1p+8);
			*(y1p+2) += *(x1p+16);
			*(y1p+3) += *(x1p+24);
			xp+=32;
			x1p+=32;
			yp+=4;
			y1p+=4;
		}
		
		for (i=s1; i<s; i++) {
			*yp += *xp;
			*y1p += *x1p;
			xp+=8;
			x1p+=8;
			yp++;
			y1p++;
		}
		
	}
	if (Y2) {
		/*if (!*Y2) {
			*Y2 = (V4_BASETYPE*) FMMV_MALLOC(FMMV, n*sizeof(V4_BASETYPE));
			for (i=0, ii=2; i<s; i++, ii+=8) {
				(*Y2)[i] = x[ii];
				(*Y2)[i+s] = x[ii+4];
			}
		}
		else {
			for (i=0, ii=2; i<s; i++, ii+=8) {
				(*Y2)[i] += x[ii];
				(*Y2)[i+s] += x[ii+4];
			}
		}*/
		xp = x+2;
		x1p = x+6;
		yp = *Y2;
		y1p = (*Y2)+s;
		for (i=0; i<s0; i++) {
			*yp += *xp;
			*(yp+1) += *(xp+8);
			*(yp+2) += *(xp+16);
			*(yp+3) += *(xp+24);
			*y1p += *x1p;
			*(y1p+1) += *(x1p+8);
			*(y1p+2) += *(x1p+16);
			*(y1p+3) += *(x1p+24);
			xp+=32;
			x1p+=32;
			yp+=4;
			y1p+=4;
		}
		
		for (i=s1; i<s; i++) {
			*yp += *xp;
			*y1p += *x1p;
			xp+=8;
			x1p+=8;
			yp++;
			y1p++;
		}
		
	}
	if (Y3) {
		/*if (!*Y3) {
			*Y3 = (V4_BASETYPE*) FMMV_MALLOC(FMMV, n*sizeof(V4_BASETYPE));
			for (i=0, ii=3; i<s; i++, ii+=8) {
				(*Y3)[i] = x[ii];
				(*Y3)[i+s] = x[ii+4];
			}
		}
		else {
			for (i=0, ii=3; i<s; i++, ii+=8) {
				(*Y3)[i] += x[ii];
				(*Y3)[i+s] += x[ii+4];
			}
		}*/
		xp = x+3;
		x1p = x+7;
		yp = *Y3;
		y1p = (*Y3)+s;
		for (i=0; i<s0; i++) {
			*yp += *xp;
			*(yp+1) += *(xp+8);
			*(yp+2) += *(xp+16);
			*(yp+3) += *(xp+24);
			*y1p += *x1p;
			*(y1p+1) += *(x1p+8);
			*(y1p+2) += *(x1p+16);
			*(y1p+3) += *(x1p+24);
			xp+=32;
			x1p+=32;
			yp+=4;
			y1p+=4;
		}
		
		for (i=s1; i<s; i++) {
			*yp += *xp;
			*y1p += *x1p;
			xp+=8;
			x1p+=8;
			yp++;
			y1p++;
		}
		
	}
#undef Y0	
#undef Y1	
#undef Y2	
#undef Y3	
	
}


void perm_simd4(const int n, const int *vec, const V4_BASETYPE *x, V4_BASETYPE *y)
{
	int i;
	for (i=0; i<n; i++) {
		V4_SETITEM(y, i<<2, V4_GETITEM(x, vec[i]<<2));
	}	
}	

void perm_inv_simd4(const int n, const int *vec, const V4_BASETYPE *x, V4_BASETYPE *y)
{
	int i;
	for (i=0; i<n; i++) {
		V4_SETITEM(y, vec[i]<<2, V4_GETITEM(x, i<<2));
	}	
}	

void neg_simd4(const int n, const int *vec, V4_BASETYPE *x)
{
        int i;
        for (i=0; i<n; i++) {
                V4_SETITEM(x, vec[i]<<2, V4_NEG(V4_GETITEM(x, vec[i]<<2)));
        }
}

void VEC_conj(const int n, _FLOAT_ *y)
{
	int i;
	for (i=1; i<2*n+1; i+=2) {
		y[i] = -y[i];
	}
}	

void VEC_conj4_simd4(const int n, V4_BASETYPE *y)
{
	int i;
	for (i=4; i<8*n+4; i+=8) {
		V4_SETITEM(y, i, V4_NEG(V4_GETITEM(y, i)));
	}
}	


#include<math.h> // ldexpf

void scale_X_simd4(int p, V4_BASETYPE *xx, int level)
{ 
    V4_TYPE scale_p = V4_SET1(ldexpf(+1.0, level)); 
    V4_TYPE scale_m = V4_SET1(ldexpf(-1.0, level)); 
    int n, j;
    int jj = 0;

    for (n=0; n<=p; n++) {
       if(n&1) { /* n odd */
           for(j=0;j<2*n+1;j++) { // NOTE: Check this for reduced_scheme !!!
               // xx[jj] = -ldexp(xx[jj], level);
               V4_SETITEM(xx, jj, V4_MUL(V4_GETITEM(xx, jj), scale_m));
               jj+=4;
           }
       }
       else {
           for(j=0;j<2*n+1;j++) {
               //xx[jj] = +ldexp(xx[jj], level);
               V4_SETITEM(xx, jj, V4_MUL(V4_GETITEM(xx, jj), scale_p));
               jj+=4;
           }
       }
    }
}
