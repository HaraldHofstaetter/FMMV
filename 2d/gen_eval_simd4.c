/*
 * FMMV - the Fastest Multipole Method of Vienna
 * Copyright (c) 2006-2010 Harald Hofstaetter
 * Institute of Mathematics 
 * University of Vienna, Austria
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

#include<math.h>
#include"_fmmv.h"

#include "simd.h"
#if (FMM_PRECISION==0)
   #include"simd4s.h"
#else
   #include"simd4d.h"
#endif


void gen_M_base_simd4(FmmvHandle *FMMV, Box *box, V4_TYPE x, V4_TYPE y, V4_TYPE q, V4_TYPE mx, V4_TYPE my, int dipole)
{
	int p = FMMV->pM;
	V4_BASETYPE *M = box->M;
       	V4_BASETYPE*recip = FMMV->A[0];
        V4_YPE m_re, m_im, h_re, h, scale;
	
	scale = V4_SET1(ldexp(1.0, box->level));

	z_re = V4_MUL(scale, z_re);
	z_im = V4_MUL(scale, z_im);
		
	M[0] -= q; // changed sign of M[0] ( => M2M simpler ) 
	/* M[1] += 0; */   // for imaginary part of log(z) charge	
	m_re = z_re;
	m_im = z_im;
	M[2] -= q*m_re;
	M[3] -= q*m_im;
	for (j=2; j<=p; j++) {
		h_re = m_re*z_re - m_im*z_im;
		m_im = m_re*z_im + m_im*z_re;
		m_re = h_re;	
		h = q*recip[j];
		M[2*j  ] -= h*m_re;
		M[2*j+1] -= h*m_im;
	}

        if (dipole) {
	        mx = scale*mx;  // 1/z charge corresponds to dipole charge => additional factor 2^level.
	        my = scale*my;  // 1/z charge corresponds to dipole charge => additional factor 2^level.

	        M[2] += mx;
	        M[3] += my;

	        m_re = z_re;
	        m_im = z_im;
	        M[4] += mx*m_re - my*m_im;
	        M[5] += mx*m_im + my*m_re;

	        for (j=3; j<=p; j++) {
		        h_re = m_re*z_re - m_im*z_im;
		        m_im = m_re*z_im + m_im*z_re;
		        m_re = h_re;	
		        M[2*j  ] += mx*m_re - my*m_im;
		        M[2*j+1] += mx*m_im + my*m_re;
	        }
        }

}

void gen_L_base_simd4(FmmvHandle *FMMV, Box *target, V4_TYPE x, V4_TYPE y, V4_TYPE q, V4_TYPE mx, V4_TYPE my, int dipole)
{
}




void eval_L_base_simd4(FmmvHandle *FMMV, Box *box, V4_TYPE x, V4_TYPE y, V4_TYPE *pot, V4_TYPE *dx, V4_TYPE *dy, int grad)
{
}


void eval_M_base_simd4(FmmvHandle *FMMV, Box *source, V4_TYPE x, V4_TYPE y, V4_TYPE *pot, V4_TYPE *dx, V4_TYPE *dy, int grad)
{
}



void gen_L_eval_M_base_simd4(FmmvHandle *FMMV, Box *list3, 
                        V4_TYPE x, V4_TYPE y, V4_TYPE q, V4_TYPE mx, V4_TYPE my,
                        V4_TYPE *pot, V4_TYPE *dx, V4_TYPE *dy, 
                        int dipole, int grad)

{

}
