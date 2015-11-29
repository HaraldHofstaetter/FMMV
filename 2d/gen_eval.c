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


void gen_M_base(FmmvHandle *FMMV, Box *box, _FLOAT_ z_re, _FLOAT_ z_im, _FLOAT_ q, _FLOAT_ mx, _FLOAT_ my, int dipole)
{
	int p = FMMV->pM;
        _FLOAT_ beta = FMMV->beta;
	_FLOAT_ *M = box->M;
       	_FLOAT_*recip = FMMV->recip;
        _FLOAT_ m_re, m_im, h_re, h, scale;
        int j;

	scale = ldexp(1.0, box->level);
        
        if (beta==0.0) {

	    z_re = scale*z_re;
	    z_im = scale*z_im;
		
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
        else {
           if (dipole) {
               _FLOAT_ r;
               _FLOAT_ f0,fm,fp;
               _FLOAT_ h0,hm,hp;
                _FLOAT_ RR[2*FMM_P_MAX+4];
                double *FF = FMMV->FF;

                r = SQRT(z_re*z_re + z_im*z_im);
                z_re *= scale;   
                z_im *= scale;
                bessel_i_scaled(r*beta, RR, p+2, 2*p+4);

                f0 = 2.0*q*r*beta;
                fm = beta*beta*beta*r/scale;
                fp = beta*r*scale;

	        M[0] -= 0.5*f0*RR[0]; 

                hp = fp*FF[0]*RR[0];
                M[2] += hp*mx;
                M[3] += hp*my;
                
	        m_re = z_re;
	        m_im = z_im;                
                h0 = f0*RR[1];
                hm = fm*RR[1];
                hp = fp*FF[1]*RR[1];
                M[0] +=hm*(mx*z_re + my*z_im); 
	        M[2] -= h0*m_re;
	        M[3] -= h0*m_im;
                M[4] += hp*(mx*m_re - my*m_im);
                M[5] += hp*(mx*m_im + my*m_re);
                
	        for (j=2; j<p; j++) {
		    h_re = m_re*z_re - m_im*z_im;
		    m_im = m_re*z_im + m_im*z_re;
		    m_re = h_re;	
                    h0 = f0*FF[j-1]*RR[j];
                    hm = fm*FF[j-2]*RR[j];
                    hp = fp*FF[j]*RR[j];
                    M[2*j-2] += hm*(mx*m_re + my*m_im);
                    M[2*j-1] += hm*(mx*m_im - my*m_re);
		    M[2*j  ] -= h0*m_re;
		    M[2*j+1] -= h0*m_im;
                    M[2*j+2] += hp*(mx*m_re - my*m_im);
                    M[2*j+3] += hp*(mx*m_im + my*m_re);
	        }
                if (p>=2) {
                    j = p;
		    h_re = m_re*z_re - m_im*z_im;
		    m_im = m_re*z_im + m_im*z_re;
		    m_re = h_re;	
                    h0 = f0*FF[j-1]*RR[j];
                    hm = fm*FF[j-2]*RR[j];
                    M[2*j-2] += hm*(mx*m_re + my*m_im);
                    M[2*j-1] += hm*(mx*m_im - my*m_re);
		    M[2*j  ] -= h0*m_re;
		    M[2*j+1] -= h0*m_im;
	        }
                
           }
           else {
                _FLOAT_ r, f0;
                _FLOAT_ RR[2*FMM_P_MAX+2];
                double *FF = FMMV->FF;

                r = SQRT(z_re*z_re + z_im*z_im);
                z_re *= scale;   
                z_im *= scale;
                bessel_i_scaled(r*beta, RR, p+1, 2*p+2);
                f0 = 2.0*q*r*beta;

	        M[0] -= 0.5*f0*RR[0]; 
	        /* M[1] += 0; */   
	        m_re = z_re;
	        m_im = z_im;
                h = f0*RR[1];
	        M[2] -= h*m_re;
	        M[3] -= h*m_im;
	        for (j=2; j<=p; j++) {
		    h_re = m_re*z_re - m_im*z_im;
		    m_im = m_re*z_im + m_im*z_re;
		    m_re = h_re;	
                    h = f0*FF[j-1]*RR[j];
		    M[2*j  ] -= h*m_re;
		    M[2*j+1] -= h*m_im;
                }
           }
     }
}	

void gen_L_base(FmmvHandle *FMMV, Box *target, _FLOAT_ z_re, _FLOAT_ z_im, _FLOAT_ q, _FLOAT_ mx, _FLOAT_ my, int dipole)
{
	int p = FMMV->pL;
        _FLOAT_ beta = FMMV->beta;
	_FLOAT_ *L = target->L;
       	_FLOAT_* recip = FMMV->recip;
       	_FLOAT_ m_re, m_im, h_re, h, inv_scale;
        int j;

	inv_scale = ldexp(1.0, -target->level);

        if (beta==0.0) {
       	    _FLOAT_ r2, one_over_r2, one_over_scale_r2;

	    r2 = z_re*z_re + z_im*z_im;
	    L[0] += q*0.5*LOG(r2);
	    /* L[1] += q*atan2(z_im, z_re); */
	    one_over_r2 = RECIP(r2);
	    one_over_scale_r2 = inv_scale*one_over_r2;
        
	    m_re =  z_re*one_over_scale_r2; /* m = 1/z = conj(z)/abs(z)^2 */
            m_im = -z_im*one_over_scale_r2;
	    L[2*1  ] -= q*m_re;
	    L[2*1+1] -= q*m_im;
	    for (j=2; j<=p; j++) {  /* m = m/z = m*conj(z)/abs(z)^2 */
		    h_re = (m_re*z_re + m_im*z_im)*one_over_scale_r2;
		    m_im = (m_im*z_re - m_re*z_im)*one_over_scale_r2;
		    m_re = h_re;	
		    h = q*recip[j];
		    L[2*j  ] -= h*m_re;
		    L[2*j+1] -= h*m_im;
	    }
		
            if (dipole) {
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
       else {
             if (dipole) {
                _FLOAT_ r, one_over_scale_r2;
                _FLOAT_ f0,fm,fp;
                _FLOAT_ h0,hm,hp;
                _FLOAT_ RR[FMM_P_MAX+1];
                double *FF_INV = FMMV->FF_INV;
    
                r = SQRT(z_re*z_re + z_im*z_im);
	        one_over_scale_r2 = inv_scale*RECIP(r*r); 

                f0 = 2.0*q/(r*beta);
                fm = 1.0/(inv_scale*beta*r);
                fp = inv_scale*beta/r;

                bessel_k_scaled(r*beta, RR, p+1);

	        L[0] -= 0.5*f0*RR[0] + q*FMMV->k0_correction;
	        /* L[1] += 0 */

                hp = fp*RR[0]*FF_INV[1];
                L[2] -= hp*mx;
                L[3] += hp*my;
            
	        m_re =  z_re*one_over_scale_r2; /* m = 1/z = conj(z)/abs(z)^2 */
                m_im = -z_im*one_over_scale_r2;
	        
                h0 = f0*RR[1]*FF_INV[1];                
                hm = fm*RR[1]*FF_INV[0];
                hp = fp*RR[1]*FF_INV[2];
                L[2*1-2] -= hm*(mx*m_re - my*m_im);
                L[2*1-1] -= hm*(mx*m_im + my*m_re);
	        L[2*1  ] -= h0*m_re;
	        L[2*1+1] -= h0*m_im;
                L[2*1+2] -= hp*(mx*m_re + my*m_im);
                L[2*1+3] -= hp*(mx*m_im - my*m_re);
                
	        for (j=2; j<p; j++) {  /* m = m/z = m*conj(z)/abs(z)^2 */
		    h_re = (m_re*z_re + m_im*z_im)*one_over_scale_r2;
		    m_im = (m_im*z_re - m_re*z_im)*one_over_scale_r2;
		    m_re = h_re;	
		    h0 = f0*RR[j]*FF_INV[j];
                    hm = fm*RR[j]*FF_INV[j-1];
                    hp = fp*RR[j]*FF_INV[j+1];
                    L[2*j-2] -= hm*(mx*m_re - my*m_im);
                    L[2*j-1] -= hm*(mx*m_im + my*m_re);
		    L[2*j  ] -= h0*m_re;
		    L[2*j+1] -= h0*m_im;
                    L[2*j+2] -= hp*(mx*m_re + my*m_im);
                    L[2*j+3] -= hp*(mx*m_im - my*m_re);
	        }
                if (p>=2) {
                    j = p;
		    h_re = (m_re*z_re + m_im*z_im)*one_over_scale_r2;
		    m_im = (m_im*z_re - m_re*z_im)*one_over_scale_r2;
		    m_re = h_re;	
		    h0 = f0*RR[j]*FF_INV[j];
                    hm = fm*RR[j]*FF_INV[j-1];
                    L[2*j-2] -= hm*(mx*m_re - my*m_im);
                    L[2*j-1] -= hm*(mx*m_im + my*m_re);
		    L[2*j  ] -= h0*m_re;
		    L[2*j+1] -= h0*m_im;
	        }
             }
             else {
                _FLOAT_ r, one_over_scale_r2, f0, h0;
                _FLOAT_ RR[FMM_P_MAX+1];
                double *FF_INV = FMMV->FF_INV;
    
                r = SQRT(z_re*z_re + z_im*z_im);
	        one_over_scale_r2 = inv_scale*RECIP(r*r); 
                f0 = 2.0*q/(r*beta);
    
                bessel_k_scaled(r*beta, RR, p+1);
    
	        L[0] -= 0.5*f0*RR[0] + q*FMMV->k0_correction;
	        /* L[1] += 0 */
        
	        m_re =  z_re*one_over_scale_r2; /* m = 1/z = conj(z)/abs(z)^2 */
                m_im = -z_im*one_over_scale_r2;
	        h0 = f0*RR[1]*FF_INV[1];
	        L[2*1  ] -= h0*m_re;
	        L[2*1+1] -= h0*m_im;
	        for (j=2; j<=p; j++) {  /* m = m/z = m*conj(z)/abs(z)^2 */
		    h_re = (m_re*z_re + m_im*z_im)*one_over_scale_r2;
		    m_im = (m_im*z_re - m_re*z_im)*one_over_scale_r2;
		    m_re = h_re;	
		    h0 = f0*RR[j]*FF_INV[j];
		    L[2*j  ] -= h0*m_re;
		    L[2*j+1] -= h0*m_im;
	        }
             }
             
      }
}


void eval_L_base(FmmvHandle *FMMV, Box *box, _FLOAT_ z_re, _FLOAT_ z_im, _FLOAT_ *pot, _FLOAT_ *dx, _FLOAT_ *dy, int grad)
{
	int p = FMMV->pL;
        _FLOAT_ beta = FMMV->beta;
	_FLOAT_ *L = box->L;	
        _FLOAT_ f_re, f_im, h_re, scale;
        int j;
	
        scale = ldexp(1.0, box->level); // 2^level
	
        if (beta==0.0) {
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
	    *pot = f_re;

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
        else {
          if (grad) {
               _FLOAT_ r;
               _FLOAT_ m_re, m_im;
               _FLOAT_ pot_re;
               _FLOAT_ dm_re, dm_im;
               _FLOAT_ dp_re, dp_im;
               _FLOAT_ f0,fm,fp;
               _FLOAT_ h0,hm,hp;
                _FLOAT_ RR[2*FMM_P_MAX+4];
                double *FF = FMMV->FF;


                r = SQRT(z_re*z_re + z_im*z_im);
                f0 = r*beta;
                fm = 0.5*beta*beta*beta*r/scale;
                fp = 0.5*beta*r*scale;

                z_re *= scale;   
                z_im *= scale;
                bessel_i_scaled(r*beta, RR, p+2, 2*p+4);

                hp = FF[1]*RR[0];
                hm = RR[0];
	        pot_re = -RR[0]*L[0];
                dp_re = hp*L[2];
                dp_im = hp*L[3];
                
	        m_re = z_re;
	        m_im = z_im;               

                h0 = FF[1]*RR[1];
                hm = 2.0*FF[0]*RR[1]; // Note factor 2 !!!
                hp = FF[2]*RR[1];
                pot_re -= h0*(m_re*L[2] - m_im*L[3]);
                
                dm_re = hm*m_re*L[0];
                dm_im = hm*m_im*L[0]; 
                dp_re += hp*(m_re*L[4] - m_im*L[5]);  
                dp_im += hp*(m_re*L[5] + m_im*L[4]); 
                
	        for (j=2; j<p; j++) {
		    h_re = m_re*z_re - m_im*z_im;
		    m_im = m_re*z_im + m_im*z_re;
		    m_re = h_re;	
                    h0 = FF[j]*RR[j];
                    hm = FF[j-1]*RR[j];
                    hp = FF[j+1]*RR[j];
                    pot_re -= h0*(m_re*L[2*j] - m_im*L[2*j+1]); 
                    
                    dm_re += hm*(m_re*L[2*j-2] - m_im*L[2*j-1]);
                    dm_im += hm*(m_re*L[2*j-1] + m_im*L[2*j-2]);
                    dp_re += hp*(m_re*L[2*j+2] - m_im*L[2*j+3]);
                    dp_im += hp*(m_re*L[2*j+3] + m_im*L[2*j+2]);
	        }
                if (p>=2) {
                    j = p;
		    h_re = m_re*z_re - m_im*z_im;
		    m_im = m_re*z_im + m_im*z_re;
		    m_re = h_re;	
                    h0 = FF[j]*RR[j];
                    hm = FF[j-1]*RR[j];
                    pot_re -= h0*(m_re*L[2*j] - m_im*L[2*j+1]); 
                    dm_re += hm*(m_re*L[2*j-2] - m_im*L[2*j-1]);
                    dm_im += hm*(m_re*L[2*j-1] + m_im*L[2*j-2]);
                }
                *pot = f0*pot_re;
                *dx = -(fm*dm_re + fp*dp_re);
                *dy = -fm*dm_im + fp*dp_im;
          }
          else {
            _FLOAT_ r, r_beta, h;
            _FLOAT_ RR[2*FMM_P_MAX+2];
            double *FF = FMMV->FF;

            r = SQRT(z_re*z_re + z_im*z_im);
	    z_re *= scale;
	    z_im *= scale;
            bessel_i_scaled(r*beta, RR, p+1, 2*p+2);
            r_beta = r*beta;
	    /* Horner scheme */
            
            h = r_beta*FF[p]*RR[p]; /* *recip[j]; */
	    f_re = h*L[2*p];
	    f_im = h*L[2*p+1];
	    for (j=p-1; j>=0; j--) {
                h = r_beta*FF[j]*RR[j]; /* *recip[j]; */
		h_re = h*L[2*j  ] + f_re*z_re - f_im*z_im;
		f_im = h*L[2*j+1] + f_re*z_im + f_im*z_re;
		f_re = h_re;
	    }
	    *pot = -f_re;  

        
          }
        }        

}

void eval_M_base(FmmvHandle *FMMV, Box *source, _FLOAT_ z_re, _FLOAT_ z_im, _FLOAT_ *pot, _FLOAT_ *dx, _FLOAT_ *dy, int grad)
{

	int p = FMMV->pM;
        _FLOAT_ beta = FMMV->beta;
	_FLOAT_ *M = source->M;
       	_FLOAT_ f_re, f_im, h_re, r2, one_over_r2, one_over_scale_r2, inv_scale;
        int j;

	inv_scale = ldexp(1.0, -source->level); /* inv_scale = 1/scale */
        if (beta==0.0) {
	    r2 = z_re*z_re + z_im*z_im;
	    f_re = 0.5*LOG(r2);
	    // f_im = atan2(z_im, z_re);
	    *pot = -M[0]*f_re; /* + M[1]*f_im; */  // changed sign of M[0] ...
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
	    *pot += (f_re*z_re + f_im*z_im)*one_over_scale_r2;

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
		f_re = h_re;
		*dx = -(f_re*z_re + f_im*z_im)*one_over_r2;
		*dy = -(f_re*z_im - f_im*z_re)*one_over_r2;  
           } 
       }
       else {
          if (grad) {
                _FLOAT_ r, one_over_scale_r2;
                _FLOAT_ f0,fm,fp;
                _FLOAT_ h0,hm,hp;
               _FLOAT_ m_re, m_im;
               _FLOAT_ pot_re;
                _FLOAT_ dm_re, dm_im;
                _FLOAT_ dp_re, dp_im;
                _FLOAT_ RR[FMM_P_MAX+1];
                double *FF_INV = FMMV->FF_INV;

                r = SQRT(z_re*z_re + z_im*z_im);

	        one_over_scale_r2 = inv_scale*RECIP(r*r); 

                f0 = 2.0/(r*beta);
                fm = 1.0/(inv_scale*beta*r);
                fp = inv_scale*beta/r;

                bessel_k_scaled(r*beta, RR, p+1);
            
                h0 = 0.5*RR[0]; //+ FMMV->k0_correction;
        	pot_re = h0*(-M[0]);  // changed sign of M[0] ...

                hp = RR[0]*FF_INV[1];
                dp_re = -hp*M[2];
                dp_im = +hp*M[3];
            
	        m_re =  z_re*one_over_scale_r2; /* m = 1/z = conj(z)/abs(z)^2 */
                m_im = -z_im*one_over_scale_r2;
                
                h0 = RR[1]*FF_INV[1];                
                
                hm = RR[1]*FF_INV[0];
                hp = 2*RR[1]*FF_INV[2];
                pot_re -= h0*(m_re*M[2*1] - m_im*M[2*1+1]);   
                
                dm_re = -hm*(m_re*M[2*1-2] /* - m_im*M[2*1-1] */);
                dm_im = -hm*(/* m_re*M[2*1-1] */ + m_im*M[2*1-2] );

                dp_re -= hp*(m_re*M[2*1+2] - m_im*M[2*1+3]);
                dp_im += hp*(m_re*M[2*1+3] + m_im*M[2*1+2]);

	        for (j=2; j<p; j++) {  /* m = m/z = m*conj(z)/abs(z)^2 */
		    h_re = (m_re*z_re + m_im*z_im)*one_over_scale_r2;
		    m_im = (m_im*z_re - m_re*z_im)*one_over_scale_r2;
		    m_re = h_re;	
		    h0 = j*RR[j]*FF_INV[j];
                    hm = (j-1)*RR[j]*FF_INV[j-1];
                    hp = (j+1)*RR[j]*FF_INV[j+1];
                    pot_re -= h0*(m_re*M[2*j] - m_im*M[2*j+1]);
                    
                    dm_re -= hm*(m_re*M[2*j-2] - m_im*M[2*j-1]);
                    dm_im -= hm*(m_re*M[2*j-1] + m_im*M[2*j-2]);

                    dp_re -= hp*(m_re*M[2*j+2] - m_im*M[2*j+3]);
                    dp_im += hp*(m_re*M[2*j+3] + m_im*M[2*j+2]);

	        }
                if (p>=2) {
                    j = p;
		    h_re = (m_re*z_re + m_im*z_im)*one_over_scale_r2;
		    m_im = (m_im*z_re - m_re*z_im)*one_over_scale_r2;
		    m_re = h_re;	
		    h0 = j*RR[j]*FF_INV[j];
                    hm = (j-1)*RR[j]*FF_INV[j-1];
                    pot_re -= h0*(m_re*M[2*j] - m_im*M[2*j+1]);
                    dm_re -= hm*(m_re*M[2*j-2] - m_im*M[2*j-1]);
                    dm_im -= hm*(m_re*M[2*j-1] + m_im*M[2*j-2]);
                }
                *pot = f0*pot_re - M[0]*FMMV->k0_correction;
                *dx = -(fm*dm_re + fp*dp_re);
                *dy = +fm*dm_im + fp*dp_im;
          }
          else {
            _FLOAT_ r, h, one_over_scale_r2, one_over_r_beta;
            _FLOAT_ RR[FMM_P_MAX+1];
            double *FF_INV = FMMV->FF_INV;

            r = SQRT(z_re*z_re + z_im*z_im);
	    one_over_scale_r2 = inv_scale*RECIP(r*r); 
            one_over_r_beta = 1.0/(r*beta);

            bessel_k_scaled(r*beta, RR, p+1);

            h = one_over_r_beta*RR[0] + FMMV->k0_correction;
	    // f_im = atan2(z_im, z_re);
	    *pot = h*(-M[0]); // - M[1]);  // changed sign of M[0] ...
	    one_over_r2 = RECIP(r*r);
	    one_over_scale_r2 = inv_scale*one_over_r2;

	    /* Horner scheme */ 
	    h = one_over_r_beta*RR[p]*FF_INV[p-1];
	    f_re = h*M[2*p];
	    f_im = h*M[2*p+1];
	    for (j=p-1; j>=1; j--) {
		h = one_over_r_beta*RR[j]*FF_INV[j-1];
		h_re = (h*M[2*j  ] + (f_re*z_re + f_im*z_im)*one_over_scale_r2);
		f_im = (h*M[2*j+1] + (f_im*z_re - f_re*z_im)*one_over_scale_r2);
		f_re = h_re;
	    }
	    *pot -= (f_re*z_re + f_im*z_im)*one_over_scale_r2;
          }
       } 
}


void gen_L_eval_M_base(FmmvHandle *FMMV, Box *list3, 
                        _FLOAT_ z_re, _FLOAT_ z_im, 
                        _FLOAT_ q, _FLOAT_ mx, _FLOAT_ my,
                        _FLOAT_ *pot, _FLOAT_ *dx, _FLOAT_ *dy,
                        int dipole, int grad)
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
	f_re = 0.5*LOG(r2);
	// f_im = atan2(z_im, z_re);
	L[0] += q*f_re;
	*pot = -M[0]*f_re; /* + M[1]*f_im; */  // changed sign of M[0] ...
	one_over_r2 = RECIP(r2);
	one_over_scale_r2 = inv_scale*one_over_r2;

	m_re =  z_re*one_over_scale_r2; /* m = 1/z = conj(z)/abs(z)^2 */
        m_im = -z_im*one_over_scale_r2;
	L[2*1  ] -= q*m_re;
	L[2*1+1] -= q*m_im;
	for (j=2; j<=p; j++) {  /* m = m/z = m*conj(z)/abs(z)^2 */
		h_re = (m_re*z_re + m_im*z_im)*one_over_scale_r2;
		m_im = (m_im*z_re - m_re*z_im)*one_over_scale_r2;
		m_re = h_re;	
		h = q*recip[j];
		L[2*j  ] -= h*m_re;
		L[2*j+1] -= h*m_im;
	}

	f_re = M[2*p];  /* TODO: merge this with gen_L ... */
	f_im = M[2*p+1];
	for (j=p-1; j>=1; j--) {
		h_re = (M[2*j  ] + (f_re*z_re + f_im*z_im)*one_over_scale_r2);
		f_im = (M[2*j+1] + (f_im*z_re - f_re*z_im)*one_over_scale_r2);
		f_re = h_re;
	}
	*pot += (f_re*z_re + f_im*z_im)*one_over_scale_r2;
  
        if (dipole) {
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
		f_re = h_re;
		*dx = (f_re*z_re + f_im*z_im)*one_over_r2;
		*dy = (f_re*z_im - f_im*z_re)*one_over_r2;  
       } 



}

