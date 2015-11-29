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


void VEC_mul(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y)
{
	int i;
	for (i=0; i<n; i++) {
		y[i] = d[i]*x[i];
	}
}	

void VEC_scale(const int n, const _FLOAT_ d, const _FLOAT_ *x, _FLOAT_ *y)
{
	int i;
	for (i=0; i<n; i++) {
		y[i] = d*x[i];
	}
}	

void VEC_add(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y)
{
	int i;
	for (i=0; i<n; i++) {
		y[i] = d[i]+x[i];
	}
}	

_FLOAT_* VEC_add2(FmmvHandle *FMMV, const int n, const _FLOAT_ *x, _FLOAT_ *y)
{
	int i;
	if (y) {
		for (i=0; i<n; i++) {
			y[i] += x[i];
		}
	}
	else {
		y = (_FLOAT_*) FMMV_MALLOC(FMMV, n*sizeof(_FLOAT_));
		memcpy(y, x, n*sizeof(_FLOAT_));
	}		
	return y;
}	


void VEC_sub(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y)
{
	int i;
	for (i=0; i<n; i++) {
		y[i] = d[i]-x[i];
	}
}	

		
void VEC_mul_c(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y)
{
	int i;
	for (i=0; i<2*n; i+=2) {
		y[i] = d[i]*x[i] - d[i+1]*x[i+1];
		y[i+1] = d[i+1]*x[i] + d[i]*x[i+1];
	}
}	

void VEC_mul_c_rrii(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y)
{
	int i;
	const _FLOAT_ *di=d+n;
	const _FLOAT_ *xi=x+n;
	_FLOAT_ *yi=y+n;

	for (i=0; i<n; i++) {
		y[i] = d[i]*x[i] - di[i]*xi[i];
		yi[i] = di[i]*x[i] + d[i]*xi[i];
	}
}


		
void VEC_mul_cj(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y)
{
	int i;
	for (i=0; i<n*2; i+=2) {
		y[i] = d[i]*x[i] + d[i+1]*x[i+1];
		y[i+1] = d[i]*x[i+1] - d[i+1]*x[i];
	}
}	

void VEC_mul_cj_rrii(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y)
{
	int i;
	const _FLOAT_ *di=d+n;
	const _FLOAT_ *xi=x+n;
	_FLOAT_ *yi=y+n;

	for (i=0; i<n; i++) {
		y[i] = d[i]*x[i] + di[i]*xi[i];
		yi[i] = d[i]*xi[i] - di[i]*x[i];
	}
}	

void VEC_mul_ccj(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y, _FLOAT_ *z)
{
	_FLOAT_ a0, b0, c0, d0;

	int i;

	for (i=0; i<n*2; i+=2) {
		a0 = d[i]*x[i];
		b0 = d[i+1]*x[i+1];
		c0 = d[i]*x[i+1];
		d0 = d[i+1]*x[i];
		y[i] = a0 - b0;
		y[i+1] = c0 + d0;	
		z[i] = a0 + b0;
		z[i+1] = c0 - d0;
	}
}	

void VEC_mul_ccj_rrii(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y, _FLOAT_ *z)
{
	int i;
	const _FLOAT_ *di=d+n;
	const _FLOAT_ *xi=x+n;
	_FLOAT_ *yi=y+n;
	_FLOAT_ *zi=z+n;
	_FLOAT_ a0, b0, c0, d0;

	for (i=0; i<n; i++) {
		a0 = d[i]*x[i];
		b0 = di[i]*xi[i];
		c0 = d[i]*xi[i];
		d0 = di[i]*x[i];
		y[i] = a0 - b0;
		yi[i] = c0 + d0;	
		z[i] = a0 + b0;
		zi[i] = c0 - d0;
	}
}	

void perm(const int n, const int *vec, const _FLOAT_ *x, _FLOAT_ *y)
{
	int i;
	for (i=0; i<n; i++) {
		y[i] = x[vec[i]];
	}	
}	

void perm_inv(const int n, const int *vec, const _FLOAT_ *x, _FLOAT_ *y)
{
	int i;
	for (i=0; i<n; i++) {
		y[vec[i]] = x[i];
	}	
}	

void neg(const int n, const int *vec, _FLOAT_ *x)
{
	int i;
	for (i=0; i<n; i++) {
		x[vec[i]] = -x[vec[i]];
	}	
}	

void VEC_conj(const int n, _FLOAT_ *y)
{
	int i;
	for (i=1; i<2*n+1; i+=2) {
		y[i] = -y[i];
	}
}	

/* y = y+a*x */	
_FLOAT_ *VEC_addmul_c(FmmvHandle *FMMV, int n, _FLOAT_ *d, _FLOAT_ *x, _FLOAT_ *y)
{
	int i;
	if (!y) {
		y = (_FLOAT_*) FMMV_MALLOC(FMMV, 2*n*sizeof(_FLOAT_));
		memset(y, 0, 2*n*sizeof(_FLOAT_));

	}	

	for (i=0; i<2*n; i+=2) {
		y[i] += d[i]*x[i] - d[i+1]*x[i+1];
		y[i+1] += d[i+1]*x[i] + d[i]*x[i+1];
	}
	
	return y;
}	


void gemv(int m, int n, _FLOAT_ *A, int lda, _FLOAT_ *x, _FLOAT_ *y)
{	
	int i, j;
	_FLOAT_ *col;
	_FLOAT_ temp;

	if (n>0) {
		for (i=0; i<m; i++) {
			col = A+lda*i;
			temp = col[0]*x[0];
			for (j=1; j<n; j++) {
				temp += col[j]*x[j];
			}
			y[i] = temp;
		}
	}
}

void gemv_trans(int m, int n, _FLOAT_ *A, int lda, _FLOAT_ *x, _FLOAT_ *y)
{
	int i, j;
	_FLOAT_ *col;
	_FLOAT_ temp;

	for (j=0; j<n; j++) {
		y[j] = 0.0;
	}		
	for (i=0; i<m; i++) {
		temp = x[i];
		col = A+lda*i;
		for (j=0; j<n; j++) {
                        y[j] += temp*col[j];
		}
	}	
}	

void tpmv_upper(int n, _FLOAT_ *A, _FLOAT_ *x)
{	
	int i, j, kk, k;
	_FLOAT_ temp;

	kk = 0;
	for (j=0; j<n; j++) {
		temp = x[j]*A[kk];
		k = kk + 1;
		for (i=j+1; i<n; i++) {
			temp += x[i]*A[k];
			k += 1;
		}	
                x[j] = temp;
		kk += (n - j);
	}
}	

void tpmv_lower(int n, _FLOAT_ *A, _FLOAT_ *x)
{	
	int i, j, kk, k;
	_FLOAT_ temp;

	kk = (n*(n+1))/2-1;
	for (j=n-1; j>=0; j--) {
                temp = x[j]*A[kk]; 
		k = kk - 1;	
		for (i=j-1; i>=0; i--) {
			temp += x[i]*A[k];
			k -= 1;
		}	
                x[j] = temp;
		kk -= (j + 1) ;
	}
}	

