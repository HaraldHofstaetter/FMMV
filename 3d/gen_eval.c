/*
 * FMMV - the Fastest Multipole Method of Vienna
 * Copyright (c) 2006-2010 Harald Hofstaetter
 * Institute of Mathematics
 * University of Vienna
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

extern void (*compute_spherical_harmonics[])(int p, _FLOAT_ *Y, _FLOAT_ sin_phi, _FLOAT_ cos_phi, _FLOAT_ cos_theta);
extern void bessel_i_scaled(_FLOAT_ x, _FLOAT_ i[], int n, int n1);
extern void bessel_k_scaled(_FLOAT_ x, _FLOAT_ k[], int n);


#if (FMM_PRECISION==0)
   #define FF FF_single
   #define FF_INV FF_INV_single
   #define R R_single
#endif

extern _FLOAT_ FF[];
extern _FLOAT_ FF_INV[];
extern _FLOAT_ R[];   

void gen_M_base(FmmvHandle *FMMV, Box *box, _FLOAT_ x, _FLOAT_ y, _FLOAT_ z, _FLOAT_ q, _FLOAT_ mx, _FLOAT_ my, _FLOAT_ mz, int dipole)
{
	int p = FMMV->pM;
        _FLOAT_ beta = FMMV->beta;

	_FLOAT_ rho, sin_phi, cos_phi, one_over_rho_sin_theta, cos_theta, scale, rho_scale;
	_FLOAT_ Y[(FMM_P_MAX+1)*(FMM_P_MAX+2)];
        _FLOAT_ RR[2*FMM_P_MAX+2];
	_FLOAT_ *M;
	int k;
	_FLOAT_ mx_m, my_m, mz_m;
	_FLOAT_ mx_p, my_p, mz_p;
	
	M = box->M;

	scale = ldexp(1.0, box->level);

	rho = SQRT(x*x + y*y + z*z);
	one_over_rho_sin_theta = 1.0/SQRT(x*x + y*y);
	sin_phi = y*one_over_rho_sin_theta;
	cos_phi = x*one_over_rho_sin_theta;
	cos_theta = z/rho;

        if (beta==0.0) {
		    compute_spherical_harmonics[p](p, Y, -sin_phi, cos_phi, cos_theta);
                    rho_scale = rho*scale;
                    RR[0] = 0;
                    RR[1] = 1.0;
                    RR[2] = rho_scale;
                    for (k=2; k<=p+1; k++) {
                        RR[k+1] = RR[k]*rho_scale;
                    }  

		    if (q != 0.0) {
			    core_gen_M_L(p, RR+1, q, Y, M);
		    }
		    if (dipole) {
			    mx *= scale;
			    my *= scale;
			    mz *= scale;
			    core_gen_M_L_dipole_minus(p, RR, mx, my, mz, Y, M);
		    }
        }
        else {
                    _FLOAT_ h;
                    _FLOAT_ scale_0[FMM_P_MAX+1];
                    _FLOAT_ rho_beta = rho*beta;
		    
		    if (dipole) {
                            _FLOAT_ scale_m[FMM_P_MAX+1];
                            _FLOAT_ scale_p[FMM_P_MAX+1];
                            _FLOAT_ rho_beta2 = rho_beta*rho_beta;

		            //compute_spherical_harmonics[p+1](p+1, Y, -sin_phi, cos_phi, cos_theta);
                            //bessel_i_scaled(rho_beta, RR, p+2, 2*p+4);
		            compute_spherical_harmonics[p](p, Y, -sin_phi, cos_phi, cos_theta);
                            bessel_i_scaled(rho_beta, RR, p+1, 2*p+2);

                            rho_scale = rho*scale;
                            h = rho_scale;
    
                            scale_0[0] = RR[0]*FF[0];
                            scale_0[1] = RR[1]*h*FF[1];
		            scale_m[0] = 0;
                            scale_p[0] = RR[1]*FF[0]; 
		            scale_m[1] = RR[0]*h*FF[1]*R[3]; 
                            scale_p[1] = RR[2]*h*FF[1]*R[3];
    
		            for (k=2;k<=p; k++){	
                                h *= rho_scale;
                                scale_0[k] = RR[k]*h*FF[k];
		                scale_m[k] = RR[k-1]*h*FF[k]*R[2*k+1];
		                scale_p[k] = RR[k+1]*h*FF[k]*R[2*k+1];
	                    }	

			    mx_m = mx*beta;
			    my_m = my*beta;
			    mz_m = mz*beta;
                            mx_p = -mx_m*rho_beta2;
                            my_p = -my_m*rho_beta2;
                            mz_p = -mz_m*rho_beta2;
                            
			    core_gen_M_L_dipole_minus(p, scale_m, mx_m, my_m, mz_m, Y, M);
			    //core_gen_M_L_dipole_plus(p, scale_p, mx_p, my_p, mz_p, Y, M);
			    core_gen_M_L_dipole_plus(p-1, scale_p, mx_p, my_p, mz_p, Y, M);

		            if (q != 0.0) {
			         core_gen_M_L(p, scale_0, rho_beta*q, Y, M);
		            }
		    }
                    else if (q!=0.0) {
		            compute_spherical_harmonics[p](p, Y, -sin_phi, cos_phi, cos_theta);
                            bessel_i_scaled(rho_beta, RR, p+1, 2*p+2);

                            rho_scale = rho*scale;
                            h = rho_scale;

                            scale_0[0] = RR[0]*FF[0];
                            scale_0[1] = RR[1]*h*FF[1];

		            for (k=2;k<=p; k++){	
                                h *= rho_scale;
                                scale_0[k] = RR[k]*h*FF[k];
	                    }	

		            core_gen_M_L(p, scale_0, rho_beta*q, Y, M);
                    }
        }    
}	

void gen_L_base(FmmvHandle *FMMV, Box *target, _FLOAT_ x, _FLOAT_ y, _FLOAT_ z, _FLOAT_ q, _FLOAT_ mx, _FLOAT_ my, _FLOAT_ mz, int dipole)
{
	int p = FMMV->pL;
        _FLOAT_ beta = FMMV->beta;

	_FLOAT_ one_over_rho_scale;
	_FLOAT_ rho, one_over_rho, sin_phi, cos_phi, one_over_rho_sin_theta, cos_theta, h, scale, inv_scale;
	_FLOAT_ Y[(FMM_P_MAX+2)*(FMM_P_MAX+3)];
        _FLOAT_ RR[FMM_P_MAX+2];
	_FLOAT_ *L;
	int k;
	
	scale = ldexp(1.0, target->level);
	inv_scale = ldexp(1.0, -target->level);

	L = target->L;

	rho = SQRT(x*x + y*y + z*z);
        one_over_rho = 1.0/rho;
	one_over_rho_sin_theta = 1.0/SQRT(x*x + y*y);
	sin_phi = y*one_over_rho_sin_theta;
	cos_phi = x*one_over_rho_sin_theta;
	cos_theta = z*one_over_rho;

        if(dipole) {
		compute_spherical_harmonics[p+1](p+1, Y, -sin_phi, cos_phi, cos_theta);
        }
        else{
		compute_spherical_harmonics[p](p, Y, -sin_phi, cos_phi, cos_theta);
        }

        if (beta==0.0) {
                    one_over_rho_scale = one_over_rho*inv_scale;
                    RR[0] = one_over_rho_scale;
                    for (k=1; k<=p+1; k++) {
                        RR[k] = RR[k-1]*one_over_rho_scale;
                    }    

		    if (q!=0) {
			    core_gen_M_L(p, RR, q*scale, Y, L);
		    }	
		    if (dipole) {
	                    _FLOAT_ scale2 = scale*scale;
			    mx *= scale2; 
			    my *= scale2;
			    mz *= scale2;
			    core_gen_M_L_dipole_plus(p, RR+1, mx, my, mz, Y, L);
		    }	
        }
        else {
              _FLOAT_ scale_0[FMM_P_MAX+1];

              if (dipole) {
                    _FLOAT_ c_m, c_p;
                    _FLOAT_ scale_m[FMM_P_MAX+1];
                    _FLOAT_ scale_p[FMM_P_MAX+1];

                    bessel_k_scaled(rho*beta, RR, p+2);

                    one_over_rho_scale = one_over_rho*inv_scale;
                    h = one_over_rho_scale*one_over_rho_scale; 

                    scale_0[0] = RR[0]*one_over_rho_scale*FF_INV[0];
                    scale_0[1] = RR[1]*h*FF_INV[1]*3;
		    scale_m[0] = 0;
                    scale_p[0] = RR[1]*one_over_rho_scale*FF_INV[0]; 
		    scale_m[1] = RR[0]*h*FF_INV[1];
                    scale_p[1] = RR[2]*h*FF_INV[1];

		    for (k=2; k<=p; k++){	
                        h *= one_over_rho_scale;
                        scale_0[k] = RR[k]*h*FF_INV[k]*((double) (2*k+1));
		        scale_m[k] = RR[k-1]*h*FF_INV[k];
		        scale_p[k] = RR[k+1]*h*FF_INV[k];
	            }

                    c_m = -beta*beta*rho*scale;
                    c_p = one_over_rho*scale;
		    core_gen_M_L_dipole_minus(p, scale_m, c_m*mx, c_m*my, c_m*mz, Y, L);
		    core_gen_M_L_dipole_plus(p, scale_p, c_p*mx, c_p*my, c_p*mz, Y, L);
		    if (q!=0) {
			    core_gen_M_L(p, scale_0, q*scale, Y, L);
		    }	

              }
              else if (q!=0) {
                    bessel_k_scaled(rho*beta, RR, p+2);

                    one_over_rho_scale = one_over_rho*inv_scale;
                    h = one_over_rho_scale*one_over_rho_scale; 

                    scale_0[0] = RR[0]*one_over_rho_scale*FF_INV[0];
                    scale_0[1] = RR[1]*h*FF_INV[1]*3;

		    for (k=2; k<=p; k++){	
                        h *= one_over_rho_scale;
                        scale_0[k] = RR[k]*h*FF_INV[k]*((double) (2*k+1));
	            }
		    core_gen_M_L(p, scale_0, q*scale, Y, L);
              }
	}	
}


void eval_L_base(FmmvHandle *FMMV, Box *box,  _FLOAT_ x, _FLOAT_ y, _FLOAT_ z, _FLOAT_ *pot, _FLOAT_ *dx, _FLOAT_ *dy, _FLOAT_ *dz, int grad)
{
	int p = FMMV->pL;
       _FLOAT_ beta = FMMV->beta;


	_FLOAT_ Y[(FMM_P_MAX+1)*(FMM_P_MAX+2)];
        _FLOAT_ RR[2*FMM_P_MAX+4];
	_FLOAT_ *L;
	_FLOAT_ rho, sin_phi, cos_phi, one_over_rho_sin_theta, cos_theta, h, scale;
        _FLOAT_ rho_scale;
	int k;
	
	scale = ldexp(1.0, box->level);

	L = box->L;
	rho = SQRT(x*x + y*y + z*z);
	one_over_rho_sin_theta = 1.0/SQRT(x*x + y*y);
	sin_phi = y*one_over_rho_sin_theta;
	cos_phi = x*one_over_rho_sin_theta;
	cos_theta = z/rho;

        if (beta==0.0) {
		    compute_spherical_harmonics[p](p, Y, sin_phi, cos_phi, cos_theta);
                    rho_scale = rho*scale;
                    RR[0] = 0.0;
                    RR[1] = 1.0;
                    RR[2] = rho_scale;
                    for (k=2; k<=p+1; k++) {
                        RR[k+1] = RR[k]*rho_scale;
                    }    
          	    *pot = core_eval_L_M(p, RR+1, L, Y);
		
		    if (grad) {  
		        core_eval_L_M_grad_minus(p, RR, L, Y, dx, dy, dz);
		        *dx *= scale;
		        *dy *= scale;
		        *dz *= scale;
                    }
       }
        else {
                    _FLOAT_ rho_beta = rho*beta;
                    _FLOAT_ scale_0[FMM_P_MAX+1];

                    if (grad) {
	                _FLOAT_ dx_m, dy_m, dz_m;
                	_FLOAT_ dx_p, dy_p, dz_p;
                        _FLOAT_ rho_beta2 = rho_beta*rho_beta;
                        _FLOAT_ scale_m[FMM_P_MAX+1];
                        _FLOAT_ scale_p[FMM_P_MAX+1];
                    
		        //compute_spherical_harmonics[p+1](p+1, Y, sin_phi, cos_phi, cos_theta);
                        //bessel_i_scaled(rho_beta, RR, p+2, 2*p+4);
		        compute_spherical_harmonics[p](p, Y, sin_phi, cos_phi, cos_theta);
                        bessel_i_scaled(rho_beta, RR, p+1, 2*p+2);
    
                        rho_scale = rho*scale;
                        h = rho_scale;
    
                        scale_0[0] = RR[0]*FF[0];
                        scale_0[1] = RR[1]*h*FF[1];
		        scale_m[0] = 0;
                        scale_p[0] = RR[1]*FF[0]; 
		        scale_m[1] = RR[0]*h*FF[1]*R[3]; 
                        scale_p[1] = RR[2]*h*FF[1]*R[3];
    
		        for (k=2;k<=p; k++){	
                            h *= rho_scale;
                            scale_0[k] = RR[k]*h*FF[k];
		            scale_m[k] = RR[k-1]*h*FF[k]*R[2*k+1];
		            scale_p[k] = RR[k+1]*h*FF[k]*R[2*k+1];
	                }		
    
		        h = core_eval_L_M(p, scale_0, L, Y);
		        *pot = rho_beta*h;
    
		        core_eval_L_M_grad_minus(p, scale_m, L, Y, &dx_m, &dy_m, &dz_m);
		        //core_eval_L_M_grad_plus(p, scale_p, L, Y, &dx_p, &dy_p, &dz_p);
		        core_eval_L_M_grad_plus(p-1, scale_p, L, Y, &dx_p, &dy_p, &dz_p);

		        *dx = beta*(dx_m + rho_beta2*dx_p);
		        *dy = beta*(dy_m + rho_beta2*dy_p);
		        *dz = beta*(dz_m + rho_beta2*dz_p);
                    }
                    else  {
		        compute_spherical_harmonics[p](p, Y, sin_phi, cos_phi, cos_theta);
                        bessel_i_scaled(rho_beta, RR, p+1, 2*p+2);

                        rho_scale = rho*scale;
                        h = rho_scale;

                        scale_0[0] = RR[0]*FF[0];
                        scale_0[1] = RR[1]*h*FF[1];

		        for (k=2;k<=p; k++){	
                            h *= rho_scale;
                            scale_0[k] = RR[k]*h*FF[k];
	                }		
    
		        h = core_eval_L_M(p, scale_0, L, Y);
		        *pot = rho_beta*h;
                    }
	}	
}

void eval_M_base(FmmvHandle *FMMV, Box *source, _FLOAT_ x, _FLOAT_ y, _FLOAT_ z, _FLOAT_ *pot, _FLOAT_ *dx, _FLOAT_ *dy, _FLOAT_ *dz, int grad)
{
	int p = FMMV->pM;
        _FLOAT_ beta = FMMV->beta;

	_FLOAT_ rho, one_over_rho, sin_phi, cos_phi, one_over_rho_sin_theta, cos_theta, h, inv_scale;
        _FLOAT_ scale, one_over_rho_scale;
	int k;
	_FLOAT_ Y[(FMM_P_MAX+2)*(FMM_P_MAX+3)];
        _FLOAT_ RR[FMM_P_MAX+2];
	_FLOAT_ *M;

	inv_scale = ldexp(1.0, -source->level); /* inv_scale = 1/scale */
	scale = ldexp(1.0, source->level); 

	M = source->M;

	rho = SQRT(x*x + y*y + z*z);
	one_over_rho = 1.0/rho;
	one_over_rho_sin_theta = 1.0/SQRT(x*x + y*y);
	sin_phi = y*one_over_rho_sin_theta;
	cos_phi = x*one_over_rho_sin_theta;
	cos_theta = z*one_over_rho;

        if (grad) {   
	    compute_spherical_harmonics[p+1](p+1, Y, sin_phi, cos_phi, cos_theta);
        }
        else {
	   compute_spherical_harmonics[p](p, Y, sin_phi, cos_phi, cos_theta);
        }

        if (beta==0.0) {
            one_over_rho_scale = one_over_rho*inv_scale;
            RR[0] = one_over_rho_scale;
            for (k=1; k<=p+1; k++) {
                 RR[k] = RR[k-1]*one_over_rho_scale;
            }    
	    h = core_eval_L_M(p, RR, M, Y);
	    *pot = scale*h;

            if (grad) {
	        _FLOAT_ scale2 = -ldexp(1.0, 2*source->level); 
	        core_eval_L_M_grad_plus(p, RR+1, M, Y, dx, dy, dz);
	        *dx *= scale2;
	        *dy *= scale2;
	        *dz *= scale2;
            }    
  
        }
        else {
            _FLOAT_ scale_0[FMM_P_MAX+1];

            if (grad) {
                    _FLOAT_ c_m, c_p;
	            _FLOAT_ dx_m, dy_m, dz_m;
         	    _FLOAT_ dx_p, dy_p, dz_p;
                    _FLOAT_ scale_m[FMM_P_MAX+1];
                    _FLOAT_ scale_p[FMM_P_MAX+1];

                    bessel_k_scaled(rho*beta, RR, p+2);
                    one_over_rho_scale = one_over_rho*inv_scale;
                    h = one_over_rho_scale*one_over_rho_scale; 

                    scale_0[0] = RR[0]*one_over_rho_scale*FF_INV[0];
                    scale_0[1] = RR[1]*h*FF_INV[1]*3;
		    scale_m[0] = 0;
                    scale_p[0] = RR[1]*one_over_rho_scale*FF_INV[0]; 
		    scale_m[1] = RR[0]*h*FF_INV[1];
                    scale_p[1] = RR[2]*h*FF_INV[1];

		    for (k=2; k<=p; k++){	
                        h *= one_over_rho_scale;
                        scale_0[k] = RR[k]*h*FF_INV[k]*((double)  (2*k+1));
		        scale_m[k] = RR[k-1]*h*FF_INV[k];
		        scale_p[k] = RR[k+1]*h*FF_INV[k];
	            }		

		    h = core_eval_L_M(p, scale_0, M, Y);
                    *pot = scale*h; 

		    core_eval_L_M_grad_minus(p, scale_m, M, Y, &dx_m, &dy_m, &dz_m);
		    core_eval_L_M_grad_plus(p, scale_p, M, Y, &dx_p, &dy_p, &dz_p);
                    c_m = -beta*beta*rho*scale;
                    c_p = -one_over_rho*scale;
		    *dx = c_m*dx_m + c_p*dx_p;
		    *dy = c_m*dy_m + c_p*dy_p;
		    *dz = c_m*dz_m + c_p*dz_p;
            }
            else {
                    bessel_k_scaled(rho*beta, RR, p+2);
                    one_over_rho_scale = one_over_rho*inv_scale;
                    h = one_over_rho_scale*one_over_rho_scale; 

                    scale_0[0] = RR[0]*one_over_rho_scale*FF_INV[0];
                    scale_0[1] = RR[1]*h*FF_INV[1]*3;

		    for (k=2; k<=p; k++){	
                        h *= one_over_rho_scale;
                        scale_0[k] = RR[k]*h*FF_INV[k]*((double)  (2*k+1));
	            }		

		    h = core_eval_L_M(p, scale_0, M, Y);
                    *pot = scale*h; 

            }
      }
}


void gen_L_eval_M_base(FmmvHandle *FMMV, Box *list3, 
                        _FLOAT_ x, _FLOAT_ y, _FLOAT_ z, 
                        _FLOAT_ q, _FLOAT_ mx, _FLOAT_ my, _FLOAT_ mz, 
                        _FLOAT_ *pot, _FLOAT_ *dx, _FLOAT_ *dy, _FLOAT_ *dz,
                        int dipole, int grad)
{
	int pL = FMMV->pL;
	int pM = FMMV->pM;
	int len2M = (pM+1)*(pM+2)/2;
	int p = (pL>pM?pL:pM);
         _FLOAT_ beta = FMMV->beta;


	_FLOAT_ rho, sin_phi, cos_phi, one_over_rho_sin_theta, cos_theta, one_over_rho, h, scale, inv_scale, one_over_rho_scale;
	_FLOAT_ Y[(FMM_P_MAX+2)*(FMM_P_MAX+3)];
        _FLOAT_ RR[FMM_P_MAX+1];
	_FLOAT_ *L;
	_FLOAT_ *M;
	int k;
	
	L = list3->L;
	M = list3->M;	
	
	scale = ldexp(1.0, list3->level);
	inv_scale = ldexp(1.0, -list3->level);

	rho = SQRT(x*x + y*y + z*z);
	one_over_rho = 1.0/rho;
	one_over_rho_sin_theta = 1.0/SQRT(x*x + y*y);
	sin_phi = y*one_over_rho_sin_theta;
	cos_phi = x*one_over_rho_sin_theta;
	cos_theta = z*one_over_rho;	

	if (dipole||grad) {
	        compute_spherical_harmonics[p+1](p+1, Y, -sin_phi, cos_phi, cos_theta);
        }
        else {
	        compute_spherical_harmonics[p](p, Y, -sin_phi, cos_phi, cos_theta);
        } 

        if (beta==0.0) {
                    one_over_rho_scale = one_over_rho*inv_scale;
                    RR[0] = one_over_rho_scale;
                    for (k=1; k<=p+1; k++) {
                        RR[k] = RR[k-1]*one_over_rho_scale;
                    }    

		    if (q!=0) {
			    core_gen_M_L(pL, RR, q*scale, Y, L);
		    }	
                    if (dipole) {
	                    _FLOAT_ scale2 = ldexp(1.0, 2*list3->level);
			    mx *= scale2; 
			    my *= scale2;
			    mz *= scale2;
			    core_gen_M_L_dipole_plus(pL, RR+1, mx, my, mz, Y, L);
		    }	

  		    VEC_conj(len2M, Y);

		    h = core_eval_L_M(pM, RR, M, Y);
		    *pot = scale*h;

                    if (grad) {
	                     _FLOAT_ scale2 = -ldexp(1.0, 2*list3->level);
		             core_eval_L_M_grad_plus(pM, RR+1, M, Y, dx, dy, dz);
		             *dx *= scale2;
		             *dy *= scale2;
		             *dz *= scale2;
                    }
        }
        else {
                _FLOAT_ scale_0[FMM_P_MAX+1];

                if (dipole||grad) {
                    _FLOAT_ c_m, c_p;
                    _FLOAT_ scale_m[FMM_P_MAX+1];
                    _FLOAT_ scale_p[FMM_P_MAX+1];

                    bessel_k_scaled(rho*beta, RR, p+2);
                    one_over_rho_scale = one_over_rho*inv_scale;
                    h = one_over_rho_scale*one_over_rho_scale; 

                    scale_0[0] = RR[0]*one_over_rho_scale*FF_INV[0];
                    scale_0[1] = RR[1]*h*FF_INV[1]*3;
		    scale_m[0] = 0;
                    scale_p[0] = RR[1]*one_over_rho_scale*FF_INV[0]; 
		    scale_m[1] = RR[0]*h*FF_INV[1];
                    scale_p[1] = RR[2]*h*FF_INV[1];

		    for (k=2; k<=p; k++){	
                        h *= one_over_rho_scale;
                        scale_0[k] = RR[k]*h*FF_INV[k]*((double)  (2*k+1));
		        scale_m[k] = RR[k-1]*h*FF_INV[k];
		        scale_p[k] = RR[k+1]*h*FF_INV[k];
	            }

		    if (q!=0) {
			    core_gen_M_L(pL, scale_0, q*scale, Y, L);
		    }	
                    if (dipole) {
                            c_m = -beta*beta*rho*scale;
                            c_p = one_over_rho*scale;
			    core_gen_M_L_dipole_minus(pL, scale_m, c_m*mx, c_m*my, c_m*mz, Y, L);
			    core_gen_M_L_dipole_plus(pL, scale_p, c_p*mx, c_p*my, c_p*mz, Y, L);
		    }	

  		    VEC_conj(len2M, Y);

		    h = core_eval_L_M(pM, scale_0, M, Y);
                    *pot = scale*h; 

                    if (grad) {
	                _FLOAT_ dx_m, dy_m, dz_m;
            	        _FLOAT_ dx_p, dy_p, dz_p;

		        core_eval_L_M_grad_minus(pM, scale_m, M, Y, &dx_m, &dy_m, &dz_m);
		        core_eval_L_M_grad_plus(pM, scale_p, M, Y, &dx_p, &dy_p, &dz_p);
                        c_m = -beta*beta*rho*scale;
                        c_p = -one_over_rho*scale;
		        *dx = c_m*dx_m + c_p*dx_p;
		        *dy = c_m*dy_m + c_p*dy_p;
		        *dz = c_m*dz_m + c_p*dz_p;
                    } 
                }
                else {
                    bessel_k_scaled(rho*beta, RR, p+2);
                    one_over_rho_scale = one_over_rho*inv_scale;
                    h = one_over_rho_scale*one_over_rho_scale; 

                    scale_0[0] = RR[0]*one_over_rho_scale*FF_INV[0];
                    scale_0[1] = RR[1]*h*FF_INV[1]*3;
		    for (k=2; k<=p; k++){	
                        h *= one_over_rho_scale;
                        scale_0[k] = RR[k]*h*FF_INV[k]*((double)  (2*k+1));
	            }

		    core_gen_M_L(pL, scale_0, q*scale, Y, L);

  		    VEC_conj(len2M, Y);

		    h = core_eval_L_M(pM, scale_0, M, Y);
                    *pot = scale*h; 
                }
        }        
}

