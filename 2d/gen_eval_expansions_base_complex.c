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
#include<math.h>


void gen_M_base_complex(FmmvHandle *FMMV, Box *box, _FLOAT_ z_re, _FLOAT_ z_im, _FLOAT_ qx, _FLOAT_ qy, _FLOAT_ mx, _FLOAT_ my, int charge, int dipole)
{
	int p = FMMV->pM;
	_FLOAT_ *M = box->M;
       	_FLOAT_*recip = FMMV->recip;
        _FLOAT_ m_re, m_im, h_re, h, scale;
        int j;

	scale = ldexp(1.0, box->level);
        z_re *= scale;
        z_im *= scale;
	if (charge&&dipole) {
	    mx *= scale;  // 1/z charge corresponds to dipole charge => additional factor 2^level.
	    my *= scale;  // 1/z charge corresponds to dipole charge => additional factor 2^level.
	    M[0] -= qx; // changed sign of M[0] ( => M2M simpler ) 
	    M[1] -= qy;    // for imaginary part of log(z) charge	
	    M[2] += mx;
	    M[3] += my;
	    m_re = z_re;
	    m_im = z_im;
	    M[2] -= qx*m_re - qy*m_im;
	    M[3] -= qx*m_im + qy*m_re;
	    M[4] += mx*m_re - my*m_im;
	    M[5] += mx*m_im + my*m_re;
	    for (j=2; j<=p-1; j++) {
		h_re = m_re*z_re - m_im*z_im;
		m_im = m_re*z_im + m_im*z_re;
		m_re = h_re;	
		h = recip[j];
		M[2*j  ] -= h*(qx*m_re - qy*m_im);
		M[2*j+1] -= h*(qx*m_im + qy*m_re);
	        M[2*j+2] += mx*m_re - my*m_im;
	        M[2*j+3] += mx*m_im + my*m_re;
            }
	    h_re = m_re*z_re - m_im*z_im;
	    m_im = m_re*z_im + m_im*z_re;
	    m_re = h_re;	
	    h = recip[p];
	    M[2*p  ] -= h*(qx*m_re - qy*m_im);
	    M[2*p+1] -= h*(qx*m_im + qy*m_re);
        }
        else if (charge) {
	    M[0] -= qx; // changed sign of M[0] ( => M2M simpler ) 
	    M[1] -= qy;    // for imaginary part of log(z) charge	
	    m_re = z_re;
	    m_im = z_im;
	    M[2] -= qx*m_re - qy*m_im;
	    M[3] -= qx*m_im + qy*m_re;
	    for (j=2; j<=p; j++) {
		h_re = m_re*z_re - m_im*z_im;
		m_im = m_re*z_im + m_im*z_re;
		m_re = h_re;	
		h = recip[j];
		M[2*j  ] -= h*(qx*m_re - qy*m_im);
		M[2*j+1] -= h*(qx*m_im + qy*m_re);
	    }
        }    
        else if (dipole) {
	    mx *= scale;  // 1/z charge corresponds to dipole charge => additional factor 2^level.
	    my *= scale;  // 1/z charge corresponds to dipole charge => additional factor 2^level.

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

void gen_L_base_complex(FmmvHandle *FMMV, Box *target, _FLOAT_ z_re, _FLOAT_ z_im, _FLOAT_ qx, _FLOAT_ qy, _FLOAT_ mx, _FLOAT_ my, int charge, int dipole)
{
	int p = FMMV->pL;
	_FLOAT_ *L = target->L;
       	_FLOAT_* recip = FMMV->recip;
       	_FLOAT_ m_re, m_im, h_re, h, inv_scale;
        int j;
    	_FLOAT_ r2, one_over_r2, one_over_scale_r2;

	inv_scale = ldexp(1.0, -target->level);
	r2 = z_re*z_re + z_im*z_im;
	one_over_r2 = RECIP(r2);
	one_over_scale_r2 = inv_scale*one_over_r2;

        if (charge&&dipole) {
            m_re = 0.5*LOG(r2);
            m_im = ATAN2(z_im, z_re);
	    L[0] += qx*m_re - qy*m_im;
	    L[1] += qx*m_im + qy*m_re; 
	    m_re =  z_re*one_over_scale_r2; /* m = 1/z = conj(z)/abs(z)^2 */
            m_im = -z_im*one_over_scale_r2;
	    L[0] -= mx*m_re - my*m_im;
	    L[1] -= mx*m_im + my*m_re; 
	    L[2] -= qx*m_re - qy*m_im;
	    L[3] -= qx*m_im + qy*m_re; 
	    for (j=2; j<=p; j++) {  /* m = m/z = m*conj(z)/abs(z)^2 */
		    h_re = (m_re*z_re + m_im*z_im)*one_over_scale_r2;
		    m_im = (m_im*z_re - m_re*z_im)*one_over_scale_r2;
		    m_re = h_re;	
		    h = recip[j];
		    L[2*j-2] -= mx*m_re - my*m_im;
		    L[2*j-1] -= mx*m_im + my*m_re;
	            L[2*j  ] -= h*(qx*m_re - qy*m_im);
	            L[2*j+1] -= h*(qx*m_im + qy*m_re); 
	    }
	    h_re = (m_re*z_re + m_im*z_im)*one_over_scale_r2;
            m_im = (m_im*z_re - m_re*z_im)*one_over_scale_r2;
	    m_re = h_re;	
	    L[2*p  ] -= mx*m_re - my*m_im;
	    L[2*p+1] -= mx*m_im + my*m_re;
        }
        else if (charge) {
            m_re = 0.5*LOG(r2);
            m_im = ATAN2(z_im, z_re);
	    L[0] += qx*m_re - qy*m_im;
	    L[1] += qx*m_im + qy*m_re; 
	    m_re =  z_re*one_over_scale_r2; /* m = 1/z = conj(z)/abs(z)^2 */
            m_im = -z_im*one_over_scale_r2;
	    L[2] -= qx*m_re - qy*m_im;
	    L[3] -= qx*m_im + qy*m_re; 
	    for (j=2; j<=p; j++) {  /* m = m/z = m*conj(z)/abs(z)^2 */
		    h_re = (m_re*z_re + m_im*z_im)*one_over_scale_r2;
		    m_im = (m_im*z_re - m_re*z_im)*one_over_scale_r2;
		    m_re = h_re;	
		    h = recip[j];
	            L[2*j  ] -= h*(qx*m_re - qy*m_im);
	            L[2*j+1] -= h*(qx*m_im + qy*m_re); 
	    }
        }   
        else if (dipole) {
	    m_re = z_re*one_over_r2;
            m_im = -z_im*one_over_r2;
	    L[0] -= mx*m_re - my*m_im;
	    L[1] -= mx*m_im + my*m_re; 
	    for (j=1; j<=p; j++) {
		    h_re = (m_re*z_re + m_im*z_im)*one_over_scale_r2;
		    m_im = (m_im*z_re - m_re*z_im)*one_over_scale_r2;
		    m_re = h_re;	
		    L[2*j  ] -= mx*m_re - my*m_im;
		    L[2*j+1] -= mx*m_im + my*m_re;
	    }
       }    
}


void eval_L_base_complex(FmmvHandle *FMMV, Box *box, _FLOAT_ z_re, _FLOAT_ z_im, _FLOAT_ *potx, _FLOAT_ *poty, _FLOAT_ *dx, _FLOAT_ *dy, int grad)
{
	int p = FMMV->pL;
	_FLOAT_ *L = box->L;	
        _FLOAT_ f_re, f_im, h_re, scale;
        int j;
	
        scale = ldexp(1.0, box->level); // 2^level
	
	z_re = scale*z_re;
	z_im = scale*z_im;
	/* Horner scheme */
	f_re = L[2*p];
	f_im = L[2*p+1];
	for (j=p-1; j>=0; j--) {
	    h_re = L[2*j  ] + f_re*z_re - f_im*z_im;
	    f_im = L[2*j+1] + f_re*z_im + f_im*z_re;
	    f_re = h_re;
	}
	*potx = f_re;
	*poty = f_im;

        if (grad) {  
   	    f_re = p*L[2*p];
	    f_im = p*L[2*p+1];
	    for (j=p-1; j>=1; j--) {
		h_re = j*L[2*j  ] + f_re*z_re - f_im*z_im;
		f_im = j*L[2*j+1] + f_re*z_im + f_im*z_re;
		f_re = h_re;
	    }
	    *dx = scale*f_re;
	    *dy = -scale*f_im;
	}
}

void eval_M_base_complex(FmmvHandle *FMMV, Box *source, _FLOAT_ z_re, _FLOAT_ z_im, _FLOAT_ *potx, _FLOAT_ *poty, _FLOAT_ *dx, _FLOAT_ *dy, int grad)
{

	int p = FMMV->pM;
	_FLOAT_ *M = source->M;
       	_FLOAT_ f_re, f_im, h_re, r2, one_over_r2, one_over_scale_r2, inv_scale;
        int j;

	inv_scale = ldexp(1.0, -source->level); /* inv_scale = 1/scale */
	r2 = z_re*z_re + z_im*z_im;
	f_re = 0.5*LOG(r2);
	f_im = ATAN2(z_im, z_re);
	*potx = -M[0]*f_re - M[1]*f_im;   /* changed sign of M[0] ...  */
        *poty = -M[0]*f_im + M[1]*f_re;   /* TODO: CHECK SIGNS !!!!! */
	one_over_r2 = RECIP(r2);
	one_over_scale_r2 = inv_scale*one_over_r2;

	/* Horner scheme */ 
	f_re = M[2*p];
	f_im = M[2*p+1];
	for (j=p-1; j>=1; j--) {
	    h_re = (M[2*j  ] + (f_re*z_re + f_im*z_im)*one_over_scale_r2);
	    f_im = (M[2*j+1] + (f_im*z_re - f_re*z_im)*one_over_scale_r2);
	    f_re = h_re;
	}
	*potx += (f_re*z_re + f_im*z_im)*one_over_scale_r2;
	*poty += (f_im*z_re - f_re*z_im)*one_over_scale_r2;

        if (grad) {      
 	    f_re = p*M[2*p  ];
	    f_im = p*M[2*p+1];
	    for (j=p-1; j>=1; j--) {
		h_re = j*M[2*j  ] + (f_re*z_re + f_im*z_im)*one_over_scale_r2;
		f_im = j*M[2*j+1] + (f_im*z_re - f_re*z_im)*one_over_scale_r2;
		f_re = h_re;
	    }
	    h_re = M[0] + (f_re*z_re + f_im*z_im)*one_over_scale_r2;
	    f_im = M[1] + (f_im*z_re - f_re*z_im)*one_over_scale_r2; /* TODO: CHECK SIGNS !!!!!! */
	    f_re = h_re;
	    *dx = -(f_re*z_re + f_im*z_im)*one_over_r2;
	    *dy = -(f_re*z_im - f_im*z_re)*one_over_r2;  
        } 
}


void gen_L_eval_M_base_complex(FmmvHandle *FMMV, Box *list3, 
                        _FLOAT_ z_re, _FLOAT_ z_im, 
                        _FLOAT_ qx, _FLOAT_ qy, _FLOAT_ mx, _FLOAT_ my,
                        _FLOAT_ *potx, _FLOAT_ *poty, _FLOAT_ *dx, _FLOAT_ *dy,
                        int charge, int dipole, int grad)
{
	int pL = FMMV->pL;
	int pM = FMMV->pM;
	int len2M = (pM+1)*(pM+2)/2;
	int p = (pL>pM?pL:pM);
	_FLOAT_ *L = list3->L;
	_FLOAT_ *M = list3->M;
       	_FLOAT_*recip = FMMV->recip;
       	_FLOAT_ f_re, f_im, m_re, m_im, h_re, h, r2, one_over_r2, one_over_scale_r2, inv_scale;
        int j;
	
	inv_scale = ldexp(1.0, -list3->level);

        r2 = z_re*z_re + z_im*z_im;
	m_re = 0.5*LOG(r2);
	m_im = ATAN2(z_im, z_re);
	*potx = -M[0]*m_re - M[1]*m_im;   /* changed sign of M[0] ...  */
        *poty = -M[0]*m_im + M[1]*m_re;   /* TODO: CHECK SIGNS !!!!! */
	one_over_r2 = RECIP(r2);
	one_over_scale_r2 = inv_scale*one_over_r2;

        if (charge&&dipole) {
	    L[0] += qx*m_re - qy*m_im;
	    L[1] += qx*m_im + qy*m_re; 
	    m_re =  z_re*one_over_scale_r2; /* m = 1/z = conj(z)/abs(z)^2 */
            m_im = -z_im*one_over_scale_r2;
	    L[0] -= mx*m_re - my*m_im;
	    L[1] -= mx*m_im + my*m_re; 
	    L[2] -= qx*m_re - qy*m_im;
	    L[3] -= qx*m_im + qy*m_re; 
	    for (j=2; j<=p; j++) {  /* m = m/z = m*conj(z)/abs(z)^2 */
		    h_re = (m_re*z_re + m_im*z_im)*one_over_scale_r2;
		    m_im = (m_im*z_re - m_re*z_im)*one_over_scale_r2;
		    m_re = h_re;	
		    h = recip[j];
		    L[2*j-2] -= mx*m_re - my*m_im;
		    L[2*j-1] -= mx*m_im + my*m_re;
	            L[2*j  ] -= h*(qx*m_re - qy*m_im);
	            L[2*j+1] -= h*(qx*m_im + qy*m_re); 
	    }
	    h_re = (m_re*z_re + m_im*z_im)*one_over_scale_r2;
            m_im = (m_im*z_re - m_re*z_im)*one_over_scale_r2;
	    m_re = h_re;	
	    L[2*p  ] -= mx*m_re - my*m_im;
	    L[2*p+1] -= mx*m_im + my*m_re;
        }
        else if (charge) {
	    L[0] += qx*m_re - qy*m_im;
	    L[1] += qx*m_im + qy*m_re; 
	    m_re =  z_re*one_over_scale_r2; /* m = 1/z = conj(z)/abs(z)^2 */
            m_im = -z_im*one_over_scale_r2;
	    L[2] -= qx*m_re - qy*m_im;
	    L[3] -= qx*m_im + qy*m_re; 
	    for (j=2; j<=p; j++) {  /* m = m/z = m*conj(z)/abs(z)^2 */
		    h_re = (m_re*z_re + m_im*z_im)*one_over_scale_r2;
		    m_im = (m_im*z_re - m_re*z_im)*one_over_scale_r2;
		    m_re = h_re;	
		    h = recip[j];
	            L[2*j  ] -= h*(qx*m_re - qy*m_im);
	            L[2*j+1] -= h*(qx*m_im + qy*m_re); 
	    }
        }   
        else if (dipole) {
	    m_re = z_re*one_over_r2;
            m_im = -z_im*one_over_r2;
	    L[0] -= mx*m_re - my*m_im;
	    L[1] -= mx*m_im + my*m_re; 
	    for (j=1; j<=p; j++) {
		    h_re = (m_re*z_re + m_im*z_im)*one_over_scale_r2;
		    m_im = (m_im*z_re - m_re*z_im)*one_over_scale_r2;
		    m_re = h_re;	
		    L[2*j  ] -= mx*m_re - my*m_im;
		    L[2*j+1] -= mx*m_im + my*m_re;
	    }
        }    

	f_re = M[2*p];  /* TODO: merge this with gen_L ... */
	f_im = M[2*p+1];
	for (j=p-1; j>=1; j--) {
		h_re = (M[2*j  ] + (f_re*z_re + f_im*z_im)*one_over_scale_r2);
		f_im = (M[2*j+1] + (f_im*z_re - f_re*z_im)*one_over_scale_r2);
		f_re = h_re;
	}
	*potx += (f_re*z_re + f_im*z_im)*one_over_scale_r2;
	*poty += (f_im*z_re - f_re*z_im)*one_over_scale_r2;
  
        if (grad) {      
	    f_re = p*M[2*p  ];
	    f_im = p*M[2*p+1];
	    for (j=p-1; j>=1; j--) {
		h_re = j*M[2*j  ] + (f_re*z_re + f_im*z_im)*one_over_scale_r2;
		f_im = j*M[2*j+1] + (f_im*z_re - f_re*z_im)*one_over_scale_r2;
		f_re = h_re;
	    }
	    h_re = M[0] + (f_re*z_re + f_im*z_im)*one_over_scale_r2;
	    f_im = /*M[1] + */ (f_im*z_re - f_re*z_im)*one_over_scale_r2;
	    h_re = M[0] + (f_re*z_re + f_im*z_im)*one_over_scale_r2;
	    f_im = M[1] + (f_im*z_re - f_re*z_im)*one_over_scale_r2; /* TODO: CHECK SIGNS !!!!!! */
  	    f_re = h_re;
	    *dx = (f_re*z_re + f_im*z_im)*one_over_r2;
	    *dy = (f_re*z_im - f_im*z_re)*one_over_r2;  
       } 
}

