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

void Ry(int p, _FLOAT_ **blocks, _FLOAT_ *x, _FLOAT_ *y)
{
	int k,j, j1;

	gemv(1, 1, blocks[0], 1, x, y);
	k = 1;
	for (j=1; j<=p; j++) {
		j1 = j+1;
		gemv(j1, j1, blocks[2*j], j1, x+k, y+k);
		k += j1;
		gemv(j, j, blocks[2*j+1], j, x+k, y+k);
		k += j;
	}
}	

		

void Ry_pi(int p, /*_FLOAT_ *x,*/ _FLOAT_ *y)
{
	int len = (p+1)*(p+1);
	int dk, k0, k1, i;
	/* VEC_COPY(len, x, y); */
	
	k0 = 1;
	k1 = 3;
	dk = 2;
	while (k1<=len) {
		for (i=k0; i<k1; i++) {
			/* y[i] = -x[i]; */
			y[i] = -y[i];
		}
		dk += 2;	
		k0 = k1 + dk;
		dk += 2;
		k1 = k0 + dk;
	}	
	for (i=k0; i<len; i++) {
		/* y[i] = -x[i]; */
		y[i] = -y[i];
	}
}	
	
	
	


	
