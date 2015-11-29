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

#include "_fmmv.h"

void Rz_pi4(int p, _FLOAT_ *x, _FLOAT_ *y)
{
	_FLOAT_ S2 = 0.70710678118654752440;
	_FLOAT_ mS2 = -0.70710678118654752440;
	_FLOAT_ *yp = y;
	_FLOAT_ *xp = x;
	int k, n;

	for (n=0; n<=p; n++) {
		for (k=8; k<=n; k+=8) {
			yp[0] = xp[0];
			yp[1] = xp[1];
			yp[2] = S2*(xp[2]-xp[3]);
			yp[3] = S2*(xp[2]+xp[3]);
			yp[4] = -xp[5];
			yp[5] = xp[4];
			yp[6] = mS2*(xp[6]+xp[7]);
			yp[7] = S2*(xp[6]-xp[7]);
			yp[8] = -xp[8];
			yp[9] = -xp[9];
			yp[10] = mS2*(xp[10]-xp[11]);
			yp[11] = mS2*(xp[10]+xp[11]);
			yp[12] = xp[13];
			yp[13] = -xp[12];
			yp[14] = S2*(xp[14]+xp[15]);
			yp[15] = mS2*(xp[14]-xp[15]);
			yp += 16;
			xp += 16;
		}	

		switch (n+8-k) {/* Note: no break statements! */
		case 7:	
			yp[14] = S2*(xp[14]+xp[15]);
			yp[15] = mS2*(xp[14]-xp[15]);
		case 6:	
			yp[12] = xp[13];
			yp[13] = -xp[12];
		case 5:	
			yp[10] = mS2*(xp[10]-xp[11]);
			yp[11] = mS2*(xp[10]+xp[11]);
		case 4:	
			yp[8] = -xp[8];
			yp[9] = -xp[9];
		case 3:	
			yp[6] = mS2*(xp[6]+xp[7]);
			yp[7] = S2*(xp[6]-xp[7]);
		case 2:	
			yp[4] = -xp[5];
			yp[5] = xp[4];
		case 1:	
			yp[2] = S2*(xp[2]-xp[3]);
			yp[3] = S2*(xp[2]+xp[3]);
		case 0:	
			yp[0] = xp[0];
			yp[1] = xp[1];
		}
		yp += 2*(n+9-k);
		xp += 2*(n+9-k);
	}
}	

void Rz_minus_pi4(int p, _FLOAT_ *x, _FLOAT_ *y)
{
	_FLOAT_ S2 = 0.70710678118654752440;
	_FLOAT_ mS2 = -0.70710678118654752440;
	_FLOAT_ *yp = y;
	_FLOAT_ *xp = x;
	int k,n;

	for (n=0; n<=p; n++) {
		for (k=8; k<=n; k+=8) {
			yp[0] = xp[0];
			yp[1] = xp[1];
			yp[2] = S2*(xp[2]+xp[3]);
			yp[3] = mS2*(xp[2]-xp[3]);
			yp[4] = xp[5];
			yp[5] = -xp[4];
			yp[6] = mS2*(xp[6]-xp[7]);
			yp[7] = mS2*(xp[6]+xp[7]);
			yp[8] = -xp[8];
			yp[9] = -xp[9];
			yp[10] = mS2*(xp[10]+xp[11]);
			yp[11] = S2*(xp[10]-xp[11]);
			yp[12] = -xp[13];
			yp[13] = xp[12];
			yp[14] = S2*(xp[14]-xp[15]);
			yp[15] = S2*(xp[14]+xp[15]);
			yp += 16;
			xp += 16;
		}	

		switch (n+8-k) {/* Note: no break statements! */
		case 7:	
			yp[14] = S2*(xp[14]-xp[15]);
			yp[15] = S2*(xp[14]+xp[15]);
		case 6:	
			yp[12] = -xp[13];
			yp[13] = xp[12];
		case 5:	
			yp[10] = mS2*(xp[10]+xp[11]);
			yp[11] = S2*(xp[10]-xp[11]);
		case 4:	
			yp[8] = -xp[8];
			yp[9] = -xp[9];
		case 3:	
			yp[6] = mS2*(xp[6]-xp[7]);
			yp[7] = mS2*(xp[6]+xp[7]);
		case 2:	
			yp[4] = xp[5];
			yp[5] = -xp[4];
		case 1:	
			yp[2] = S2*(xp[2]+xp[3]);
			yp[3] = mS2*(xp[2]-xp[3]);
		case 0:	
			yp[0] = xp[0];
			yp[1] = xp[1];
		}
		yp += 2*(n+9-k);
		xp += 2*(n+9-k);
	}
}	

void Rz_3pi4(int p, _FLOAT_ *x, _FLOAT_ *y)
{
	_FLOAT_ S2 = 0.70710678118654752440;
	_FLOAT_ mS2 = -0.70710678118654752440;
	_FLOAT_ *yp = y;
	_FLOAT_ *xp = x;
	int k,n;
	
	for (n=0; n<=p; n++) {
		for (k=8; k<=n; k+=8) {
			yp[0] = xp[0];
			yp[1] = xp[1];
			yp[2] = mS2*(xp[2]+xp[3]);
			yp[3] = S2*(xp[2]-xp[3]);
			yp[4] = xp[5];
			yp[5] = -xp[4];
			yp[6] = S2*(xp[6]-xp[7]);
			yp[7] = S2*(xp[6]+xp[7]);
			yp[8] = -xp[8];
			yp[9] = -xp[9];
			yp[10] = S2*(xp[10]+xp[11]);
			yp[11] = mS2*(xp[10]-xp[11]);
			yp[12] = (-xp[13]);
			yp[13] = xp[12];
			yp[14] = mS2*(xp[14]-xp[15]);
			yp[15] = mS2*(xp[14]+xp[15]);
			yp += 16;
			xp += 16;
		}	

		switch (n+8-k) {/* Note: no break statements! */
		case 7:	
			yp[14] = mS2*(xp[14]-xp[15]);
			yp[15] = mS2*(xp[14]+xp[15]);
		case 6:	
			yp[12] = (-xp[13]);
			yp[13] = xp[12];
		case 5:	
			yp[10] = S2*(xp[10]+xp[11]);
			yp[11] = mS2*(xp[10]-xp[11]);
		case 4:	
			yp[8] = -xp[8];
			yp[9] = -xp[9];
		case 3:	
			yp[6] = S2*(xp[6]-xp[7]);
			yp[7] = S2*(xp[6]+xp[7]);
		case 2:	
			yp[4] = xp[5];
			yp[5] = -xp[4];
		case 1:	
			yp[2] = mS2*(xp[2]+xp[3]);
			yp[3] = S2*(xp[2]-xp[3]);
		case 0:	
			yp[0] = xp[0];
			yp[1] = xp[1];
		}
		yp += 2*(n+9-k);
		xp += 2*(n+9-k);
	}
}	

/* minus_3pi4 */
void Rz_minus_3pi4(int p, _FLOAT_ *x, _FLOAT_ *y)
{
	_FLOAT_ S2 = 0.70710678118654752440;
	_FLOAT_ mS2 = -0.70710678118654752440;
	_FLOAT_ *yp = y;
	_FLOAT_ *xp = x;
	int k,n;
	
	for (n=0; n<=p; n++) {
		for (k=8; k<=n; k+=8) {
			yp[0] = xp[0];
			yp[1] = xp[1];
			yp[2] = mS2*(xp[2]-xp[3]);
			yp[3] = mS2*(xp[2]+xp[3]);
			yp[4] = -xp[5];
			yp[5] = xp[4];
			yp[6] = S2*(xp[6]+xp[7]);
			yp[7] = mS2*(xp[6]-xp[7]);
			yp[8] = -xp[8];
			yp[9] = -xp[9];
			yp[10] = S2*(xp[10]-xp[11]);
			yp[11] = S2*(xp[10]+xp[11]);
			yp[12] = xp[13];
			yp[13] = -xp[12];
			yp[14] = mS2*(xp[14]+xp[15]);
			yp[15] = S2*(xp[14]-xp[15]);
			yp += 16;
			xp += 16;
		}	

		switch (n+8-k) {/* Note: no break statements! */
		case 7:	
			yp[14] = mS2*(xp[14]+xp[15]);
			yp[15] = S2*(xp[14]-xp[15]);
		case 6:	
			yp[12] = xp[13];
			yp[13] = -xp[12];
		case 5:	
			yp[10] = S2*(xp[10]-xp[11]);
			yp[11] = S2*(xp[10]+xp[11]);
		case 4:	
			yp[8] = -xp[8];
			yp[9] = -xp[9];
		case 3:	
			yp[6] = S2*(xp[6]+xp[7]);
			yp[7] = mS2*(xp[6]-xp[7]);
		case 2:	
			yp[4] = -xp[5];
			yp[5] = xp[4];
		case 1:	
			yp[2] = mS2*(xp[2]-xp[3]);
			yp[3] = mS2*(xp[2]+xp[3]);
		case 0:	
			yp[0] = xp[0];
			yp[1] = xp[1];
		}	
		yp += 2*(n+9-k);
		xp += 2*(n+9-k);
	}
}

void Rz_pi2(int p, _FLOAT_ *x, _FLOAT_ *y)
{
	_FLOAT_ *yp = y;
	_FLOAT_ *xp = x;
	int k,n;
	
	for (n=0; n<=p; n++) {
		for (k=4; k<=n; k+=4) {
			yp[0] = xp[0];
			yp[1] = xp[1];
			yp[2] = -xp[3];
			yp[3] = xp[2];
			yp[4] = -xp[4];
			yp[5] = -xp[5];
			yp[6] = xp[7];
			yp[7] = -xp[6];
			yp += 8;
			xp += 8;
		}	

		switch (n+4-k) {/* Note: no break statements! */
		case 3:	
			yp[6] = xp[7];
			yp[7] = -xp[6];
		case 2:	
			yp[4] = -xp[4];
			yp[5] = -xp[5];
		case 1:	
			yp[2] = -xp[3];
			yp[3] = xp[2];
		case 0:	
			yp[0] = xp[0];
			yp[1] = xp[1];
		}	
		yp += 2*(n+5-k);
		xp += 2*(n+5-k);
	}
}

void Rz_minus_pi2(int p, _FLOAT_ *x, _FLOAT_ *y)
{
	_FLOAT_ *yp = y;
	_FLOAT_ *xp = x;
	int k,n;
	
	for (n=0; n<=p; n++) {
		for (k=4; k<=n; k+=4) {
			yp[0] = xp[0];
			yp[1] = xp[1];
			yp[2] = xp[3];
			yp[3] = -xp[2];
			yp[4] = -xp[4];
			yp[5] = -xp[5];
			yp[6] = -xp[7];
			yp[7] = xp[6];
			yp += 8;
			xp += 8;
		}	
		switch (n+4-k) {/* Note: no break statements! */
		case 3:	
			yp[6] = -xp[7];
			yp[7] = xp[6];
		case 2:	
			yp[4] = -xp[4];
			yp[5] = -xp[5];
		case 1:	
			yp[2] = xp[3];
			yp[3] = -xp[2];
		case 0:	
			yp[0] = xp[0];
			yp[1] = xp[1];
		}	
		yp += 2*(n+5-k);
		xp += 2*(n+5-k);
	}
}
	

