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
#include <math.h> /* ceil, cos, sin */
#include <string.h> /* memset */

#define PI  3.1415926535897932385

void init_FFT_M2X_2( void );
void init_FFT_X2L_2( void );
void init_FFT_M2X_4( void );
void init_FFT_X2L_4( void );
void init_FFT_M2X_8( void );
void init_FFT_X2L_8( void );
void init_FFT_M2X_12( void );
void init_FFT_X2L_12( void );
void init_FFT_M2X_16( void );
void init_FFT_X2L_16( void );
void init_FFT_M2X_20( void );
void init_FFT_X2L_20( void );
void init_FFT_M2X_24( void );
void init_FFT_X2L_24( void );
void init_FFT_M2X_28( void );
void init_FFT_X2L_28( void );
void init_FFT_M2X_32( void );
void init_FFT_X2L_32( void );
void init_FFT_M2X_36( void );
void init_FFT_X2L_36( void );
void init_FFT_M2X_40( void );
void init_FFT_X2L_40( void );
void init_FFT_M2X_44( void );
void init_FFT_X2L_44( void );
void init_FFT_M2X_48( void );
void init_FFT_X2L_48( void );
void init_FFT_M2X_52( void );
void init_FFT_X2L_52( void );
void init_FFT_M2X_56( void );
void init_FFT_X2L_56( void );
void init_FFT_M2X_60( void );
void init_FFT_X2L_60( void );
void init_FFT_M2X_64( void );
void init_FFT_X2L_64( void );




static void (*init_FFT_M2X[])(void) = {
	0, /* 0 */
	init_FFT_M2X_2,
	init_FFT_M2X_4,
	0, /* 6 */
	init_FFT_M2X_8,
	0, /* 10 */
	init_FFT_M2X_12,
	0, /* 14 */
	init_FFT_M2X_16,
	0, /* 18 */
	init_FFT_M2X_20,
	0, /* 22 */
	init_FFT_M2X_24,
	0, /* 26 */
	init_FFT_M2X_28,
	0, /* 30 */
	init_FFT_M2X_32,
	0, /* 34 */
	init_FFT_M2X_36,
	0, /* 38 */
	init_FFT_M2X_40,
	0, /* 42 */
	init_FFT_M2X_44,
	0, /* 46 */
	init_FFT_M2X_48,
	0, /* 50 */
	init_FFT_M2X_52,
	0, /* 54 */
	init_FFT_M2X_56,
	0, /* 58 */
	init_FFT_M2X_60,
	0, /* 62 */
	init_FFT_M2X_64,
	0, /* 66 */
};

static void (*init_FFT_X2L[])(void) = {
	0, /* 0 */
	init_FFT_X2L_2,
	init_FFT_X2L_4,
	0, /* 6 */
	init_FFT_X2L_8,
	0, /* 10 */
	init_FFT_X2L_12,
	0, /* 14 */
	init_FFT_X2L_16,
	0, /* 18 */
	init_FFT_X2L_20,
	0, /* 22 */
	init_FFT_X2L_24,
	0, /* 26 */
	init_FFT_X2L_28,
	0, /* 30 */
	init_FFT_X2L_32,
	0, /* 34 */
	init_FFT_X2L_36,
	0, /* 38 */
	init_FFT_X2L_40,
	0, /* 42 */
	init_FFT_X2L_44,
	0, /* 46 */
	init_FFT_X2L_48,
	0, /* 50 */
	init_FFT_X2L_52,
	0, /* 54 */
	init_FFT_X2L_56,
	0, /* 58 */
	init_FFT_X2L_60,
	0, /* 62 */
	init_FFT_X2L_64,
	0, /* 66 */
};



void init_F(FmmvHandle *FMMV)
{
	int *M = FMMV->M;
	int s_eps = FMMV->s_eps;
	int pM = FMMV->pM;
	int *FFT_rep;
	int i, j;


	int alreadyInitialized[128/4]; 
	
	memset(alreadyInitialized, 0, (128/4)*sizeof(int));

	FFT_rep = FMMV_MALLOC(FMMV, s_eps*sizeof(int));
	FMMV->FFT_rep = FFT_rep;
	


	FMMV->len_F = 0;
	for (i=0; i<s_eps; i++) {
		j = M[i]/4;
		if (!alreadyInitialized[j]) {
			alreadyInitialized[j] = 1;
			if (init_FFT_M2X[j]!=0) {
				init_FFT_M2X[j]();
			}	
			if (init_FFT_X2L[j]!=0) {
				init_FFT_X2L[j]();
			}
		}
		FFT_rep[i] = ceil(((double) (2*pM+1))/((double) M[i]));
		FMMV->len_F += FFT_rep[i]*(M[i]/2+1);
	}
}	
	


void finish_F(FmmvHandle *FMMV)
{
	int s_eps = FMMV->s_eps;

	FMMV_FREE(FMMV, FMMV->FFT_rep, s_eps*sizeof(int));
}
