#include<math.h>
#include<stdio.h>
#define _FLOAT_ double
#define FMM_P_MAX 20

_FLOAT_ R[FMM_P_MAX+2];
_FLOAT_ F[2*FMM_P_MAX+8];
_FLOAT_ B[(FMM_P_MAX+2)*(FMM_P_MAX+4)+1];

#define J(n,m) ((n)*((n)+1)/2 + (m))

void init_coeffs(void)
{
	int n,m;

	/* init factorials */
	
	F[0] = 1.0;
	for (n=1; n<=2*FMM_P_MAX+4; n++) {
		F[n]=n*F[n-1];
	}
	
	for (n=0; n<=FMM_P_MAX+2; n++) {
		for (m=0; m<=n; m++) {
			B[J(n,m)] = sqrt(F[n-abs(m)]/F[n+abs(m)]);
		}
	}	

	/* init R */
	R[0] = 0.0;
	for (n=1; n<=FMM_P_MAX+1; n++) {
		R[n] = ((double) 1.0)/n;
	}
}

void compute_spherical_harmonics(int p, _FLOAT_ *Y, _FLOAT_ sin_phi, _FLOAT_ cos_phi, _FLOAT_ cos_theta)
{

	_FLOAT_	pmm;
	_FLOAT_	pm1;
	_FLOAT_	pm2;
	_FLOAT_	pml;
	_FLOAT_	c;
	_FLOAT_	s;
	_FLOAT_	h;
	_FLOAT_	alpha;
	_FLOAT_	beta;
	_FLOAT_	sqrt_1_minus_cos_theta_2 = sqrt(1.0 - cos_theta*cos_theta);

	int m, l, k, kk;

	pmm = 1.0;
	
	/* m==0: ***************/
	Y[0] = B[0] * pmm;
	Y[1] = 0.0;
	pm2 = pmm;
	pml = pmm*cos_theta;
	Y[2] = pml;
	Y[3] = 0.0 ;
	k = 1;
	for (l=2; l<=p; l++) {
		pm1 = pml;
		pml = R[l]*((2*l-1)*cos_theta*pm1 - (l-1)*pm2);
		pm2 = pm1;
		k += l;
		Y[2*k] = pml;
		Y[2*k+1] = 0.0;
	}

	/* m==1: ***************/
	pmm *= -sqrt_1_minus_cos_theta_2;
	s = sin_phi;
	c = cos_phi;
	alpha = 1-c;
	beta = s;
	h = B[2] * pmm;
	Y[4] = c*h;
	Y[5] = s*h;
	pm2 = pmm;
	pml = 3*pmm*cos_theta;
	h = B[4] * pml;
	Y[8] = c*h;
	Y[9] = s*h;
	k=4;
	for (l=3; l<=p; l++) {
		pm1 = pml;
		pml = R[l-1]*((2*l-1)*cos_theta*pm1 - l*pm2);
		pm2 = pm1;
		k += l;
		h = B[k] * pml;
		Y[2*k] = c*h;
		Y[2*k+1] = s*h;
	}

	/* 2 <=m <= p-1: ***************/
	kk = 1;
	for (m=2; m<p; m++) {
		pmm *= (1-2*m)*sqrt_1_minus_cos_theta_2;
		h = (alpha*c + beta*s);
		s = s - alpha*s + beta*c;
		c -= h;
		kk += m;
		k = kk + m;
		h = B[k] * pmm;
		Y[2*k] = c*h;
		Y[2*k+1] = s*h;
		pm2 = pmm;
		pml = (2*m+1)*pmm*cos_theta;
		k += m+1;
		h = B[k] * pml;
		Y[2*k] = c*h;
		Y[2*k+1] = s*h;
		for (l=m+2; l<=p; l++) {
			pm1 = pml;
			pml = R[l-m] * ((2*l-1)*cos_theta*pm1 - (l+m-1)*pm2);
			pm2 = pm1;
			k += l;
			h = B[k] * pml;
			Y[2*k] = c*h;
			Y[2*k+1] = s*h;
		}
	}

	/* m==p: ***************/
        pmm *= (1-2*p)*sqrt_1_minus_cos_theta_2;
        h = (alpha*c + beta*s);
        s = s - alpha*s + beta*c;
        c -= h;
	kk += p;
	k = kk + p;
        h = B[k] * pmm;
        Y[2*k] = c*h;
        Y[2*k+1] = s*h;
}


void compute_spherical_harmonics_xyz(int p, _FLOAT_ *Y, _FLOAT_ x, _FLOAT_ y, _FLOAT_ z)
{
    double sin_phi, cos_phi, cos_theta, rho, one_over_rho_sin_theta;

    rho = sqrt(x*x + y*y + z*z);
    one_over_rho_sin_theta = 1.0/sqrt(x*x + y*y);
    sin_phi = y*one_over_rho_sin_theta;
    cos_phi = x*one_over_rho_sin_theta;
    cos_theta = z/rho;


    compute_spherical_harmonics(p, Y, sin_phi, cos_phi, cos_theta);
}

void bessel_i(_FLOAT_ x, _FLOAT_ i[], int n, int n1)
{
     int l;
     _FLOAT_ one_over_x, s;

     one_over_x = 1.0/x;
     i[n1-1] = 0.0;
     i[n1-2] = 1.0;

     for (l=n1-3; l>=0; l--) {
         i[l] = i[l+2] + (2*l+3)*i[l+1]*one_over_x;
     }

     s = sinh(x)*one_over_x/i[0];
     for (l=0;l<n; l++) {
         i[l] *= s;
     }    
}

void bessel_k(_FLOAT_ x, _FLOAT_ k[], int n)
{
     int l;
     _FLOAT_ one_over_x;

     one_over_x = 1.0/x;
     k[0] = exp(-x)*one_over_x; // *1.5707963267948966192; /* pi/2 */
     k[1] = k[0]*one_over_x*(x+1.0);
     for (l=2; l<n; l++) {
         k[l] = (2*l-1)*k[l-1]*one_over_x + k[l-2];
     }
}

void compute_Q(int pp, _FLOAT_ *Q, _FLOAT_ kappa, _FLOAT_ x, _FLOAT_ y, _FLOAT_ z)
{
    double rho;
    int n,m, ii;
    double II[2*FMM_P_MAX+4];
    int p = abs(pp);
    
    rho = sqrt(x*x + y*y + z*z);
    compute_spherical_harmonics_xyz(p, Q, x, y, z);
    if (kappa==0) {
       if (pp<0) {
          for (n=0;n<=p;n++) {
              II[n] = pow(rho, -(n+1));
          }
       }
       else {
          for (n=0;n<=p;n++) {
              II[n] = pow(rho, n);
          }
       }
    }
    else {
       if (pp<0) {
           bessel_k(kappa*rho, II, p+1);
       }
       else {
           bessel_i(kappa*rho, II, p+1, 2*p+2);
       }
    }   
    
    ii = 0;
    for (n=0;n<=p;n++) {
       for (m=0;m<=n;m++) {
           Q[2*ii] *= II[n];
           Q[2*ii+1] *= II[n];
           ii++;
       }
    }
}    

#define Re(n,m) (2*((n)*((n)+1)/2 + (m)))
#define Im(n,m) (2*((n)*((n)+1)/2 + (m))+1)

void compute_dQ(int pp, _FLOAT_ *dxQ, _FLOAT_ *dyQ, _FLOAT_ *dzQ, 
                _FLOAT_ kappa, _FLOAT_ x, _FLOAT_ y, _FLOAT_ z)
{
    _FLOAT_ Y[(FMM_P_MAX+2)*(FMM_P_MAX+4)+1];

    int n,m;
    int p = abs(pp);
    

    if (kappa==0) {
    if (pp>0) {
    compute_Q(p, Y, kappa, x, y, z);
    for (n=0;n<=p;n++) {
      for (m=0;m<=n;m++) {
        if (n==0) {
          dxQ[Re(n,m)] = 0;
          dxQ[Im(n,m)] = 0;
          dyQ[Re(n,m)] = 0;
          dyQ[Im(n,m)] = 0.0;
          dzQ[Re(n,m)] = 0.0; 
          dzQ[Im(n,m)] = 0.0; 
        }
        else {
             if (m==0) {        
		dxQ[Re(n,m)] = sqrt(n*(n-1))*Y[Re(n-1,1)];
                dxQ[Im(n,m)] = 0;
		dyQ[Re(n,m)] = sqrt(n*(n-1))*Y[Im(n-1,1)];
                dyQ[Im(n,m)] = 0;
                dzQ[Re(n,m)] = n*Y[Re(n-1,0)];
                dzQ[Im(n,m)] = 0.0; 
             }
             else {
		dxQ[Re(n,m)] = 0.5*(+sqrt((n-m)*(n-m-1))*Y[Re(n-1,m+1)] - sqrt((n+m)*(n+m-1))*Y[Re(n-1,m-1)]);
		dxQ[Im(n,m)] = 0.5*(+sqrt((n-m)*(n-m-1))*Y[Im(n-1,m+1)] - sqrt((n+m)*(n+m-1))*Y[Im(n-1,m-1)]);
		dyQ[Re(n,m)] = 0.5*(+sqrt((n-m)*(n-m-1))*Y[Im(n-1,m+1)] + sqrt((n+m)*(n+m-1))*Y[Im(n-1,m-1)]);
		dyQ[Im(n,m)] = 0.5*(-sqrt((n-m)*(n-m-1))*Y[Re(n-1,m+1)] - sqrt((n+m)*(n+m-1))*Y[Re(n-1,m-1)]);
                dzQ[Re(n,m)] = sqrt((n-m)*(n+m))*Y[Re(n-1,m)];
                dzQ[Im(n,m)] = sqrt((n-m)*(n+m))*Y[Im(n-1,m)];
             }
        }
      }
    }
    }
    else { /************/
    compute_Q(-(p+1), Y, kappa, x, y, z);
    for (n=0;n<=p;n++) {
      for (m=0;m<=n;m++) {
        if (m==0) {
	   dxQ[Re(n,m)] = sqrt((n+2)*(n+1))*Y[Re(n+1,1)];
           dxQ[Im(n,m)] = 0;
	   dyQ[Re(n,m)] = sqrt((n+2)*(n+1))*Y[Im(n+1,1)];
           dyQ[Im(n,m)] = 0;
	   dzQ[Re(n,m)] = (-(n+1))*Y[Re(n+1,0)];
           dzQ[Im(n,m)] = 0;
        }
        else {
  	  dxQ[Re(n,m)] = 0.5*(+sqrt((n+m+2)*(n+m+1))*Y[Re(n+1,m+1)] - sqrt((n-m+2)*(n-m+1))*Y[Re(n+1,m-1)]);
	  dxQ[Im(n,m)] = 0.5*(+sqrt((n+m+2)*(n+m+1))*Y[Im(n+1,m+1)] - sqrt((n-m+2)*(n-m+1))*Y[Im(n+1,m-1)]);
	  dyQ[Re(n,m)] = 0.5*(+sqrt((n+m+2)*(n+m+1))*Y[Im(n+1,m+1)] + sqrt((n-m+2)*(n-m+1))*Y[Im(n+1,m-1)]);
	  dyQ[Im(n,m)] = 0.5*(-sqrt((n+m+2)*(n+m+1))*Y[Re(n+1,m+1)] - sqrt((n-m+2)*(n-m+1))*Y[Re(n+1,m-1)]);
	  dzQ[Re(n,m)] = -sqrt((n+m+1)*(n-m+1))*Y[Re(n+1,m)];
	  dzQ[Im(n,m)] = -sqrt((n+m+1)*(n-m+1))*Y[Im(n+1,m)];
        }
      }
    }  
    }      /************/
    }
    else {
    if (pp<0) {
       compute_Q(-(p+1), Y, kappa, x, y, z);
       kappa = -kappa;
    }
    else{
       compute_Q(p+1, Y, kappa, x, y, z);
    }

    for (n=0;n<=p;n++) {
      for (m=0;m<=n;m++) {
        if (n==0) {
          dxQ[Re(n,m)] = -kappa*sqrt(2)*Y[Re(n+1,1)];
          dxQ[Im(n,m)] = 0;
          dyQ[Re(n,m)] = -kappa*sqrt(2)*Y[Im(n+1,1)];
          dyQ[Im(n,m)] = 0.0;
          dzQ[Re(n,m)] = kappa*Y[Re(n+1,0)];
          dzQ[Im(n,m)] = 0.0; 
        }
        else {
             if (m==0) {        
		dxQ[Re(n,m)] = kappa/((double)(2*n+1))*(sqrt(n*(n-1))*Y[Re(n-1,1)] - sqrt((n+1)*(n+2))*Y[Re(n+1,1)]);
                dxQ[Im(n,m)] = 0;
		dyQ[Re(n,m)] = kappa/((double)(2*n+1))*(sqrt(n*(n-1))*Y[Im(n-1,1)] - sqrt((n+1)*(n+2))*Y[Im(n+1,1)]);
                dyQ[Im(n,m)] = 0;
                dzQ[Re(n,m)] = kappa/((double)(2*n+1))*(n*Y[Re(n-1,0)]+(n+1)*Y[Re(n+1,0)]);
                dzQ[Im(n,m)] = 0.0; 
             }
             else {
		dxQ[Re(n,m)] = 0.5*kappa/((double)(2*n+1))*(+sqrt((n-m)*(n-m-1))*Y[Re(n-1,m+1)] - sqrt((n+m+1)*(n+m+2))*Y[Re(n+1,m+1)] 
                                                            -sqrt((n+m)*(n+m-1))*Y[Re(n-1,m-1)] + sqrt((n-m+1)*(n-m+2))*Y[Re(n+1,m-1)]);
                                                            
		dxQ[Im(n,m)] = 0.5*kappa/((double)(2*n+1))*(+sqrt((n-m)*(n-m-1))*Y[Im(n-1,m+1)] - sqrt((n+m+1)*(n+m+2))*Y[Im(n+1,m+1)] 
                                                            -sqrt((n+m)*(n+m-1))*Y[Im(n-1,m-1)] + sqrt((n-m+1)*(n-m+2))*Y[Im(n+1,m-1)]); 
		dyQ[Re(n,m)] = .5*kappa/((double)(2*n+1))*(+sqrt((n-m)*(n-m-1))*Y[Im(n-1,m+1)]  - sqrt((n+m+1)*(n+m+2))*Y[Im(n+1,m+1)] 
                                                           +sqrt((n+m)*(n+m-1))*Y[Im(n-1,m-1)]  - sqrt((n-m+1)*(n-m+2))*Y[Im(n+1,m-1)]); 
		dyQ[Im(n,m)] = .5*kappa/((double)(2*n+1))*(-sqrt((n-m)*(n-m-1))*Y[Re(n-1,m+1)]  + sqrt((n+m+1)*(n+m+2))*Y[Re(n+1,m+1)] 
                                                           -sqrt((n+m)*(n+m-1))*Y[Re(n-1,m-1)]  + sqrt((n-m+1)*(n-m+2))*Y[Re(n+1,m-1)]); 
                dzQ[Re(n,m)] = kappa/((double)(2*n+1))*(sqrt((n-m)*(n+m))*Y[Re(n-1,m)]+sqrt((n+1-m)*(n+1+m))*Y[Re(n+1,m)]);
                dzQ[Im(n,m)] = kappa/((double)(2*n+1))*(sqrt((n-m)*(n+m))*Y[Im(n-1,m)]+sqrt((n+1-m)*(n+1+m))*Y[Im(n+1,m)]);
             }

             }
          }
      }


      }
}

void compute_dQ_nd(int pp, _FLOAT_ *dxQ, _FLOAT_ *dyQ, _FLOAT_ *dzQ, 
                _FLOAT_ dx, _FLOAT_ kappa, _FLOAT_ x, _FLOAT_ y, _FLOAT_ z)
{
    _FLOAT_ Q[(FMM_P_MAX+2)*(FMM_P_MAX+4)+1];
    _FLOAT_ L[(FMM_P_MAX+2)*(FMM_P_MAX+4)+1];
    int n,m;
    int p = abs(pp);
    
    compute_Q(pp, dxQ, kappa, x-dx, y, z);
    compute_Q(pp, Q, kappa, x+dx, y, z);
    for (n=0;n<=p;n++) {
      for (m=0;m<=n;m++) {
          dxQ[Re(n,m)] = (Q[Re(n,m)]-dxQ[Re(n,m)])/(2.0*dx);
          dxQ[Im(n,m)] = (Q[Im(n,m)]-dxQ[Im(n,m)])/(2.0*dx);
      }
    }  

    compute_Q(pp, dyQ, kappa, x, y-dx, z);
    compute_Q(pp, Q, kappa, x, y+dx, z);
    for (n=0;n<=p;n++) {
      for (m=0;m<=n;m++) {
          dyQ[Re(n,m)] = (Q[Re(n,m)]-dyQ[Re(n,m)])/(2.0*dx);
          dyQ[Im(n,m)] = (Q[Im(n,m)]-dyQ[Im(n,m)])/(2.0*dx);
      }
    }  

    compute_Q(pp, dzQ, kappa, x, y, z-dx);
    compute_Q(pp, Q, kappa, x, y, z+dx);
    for (n=0;n<=p;n++) {
      for (m=0;m<=n;m++) {
          dzQ[Re(n,m)] = (Q[Re(n,m)]-dzQ[Re(n,m)])/(2.0*dx);
          dzQ[Im(n,m)] = (Q[Im(n,m)]-dzQ[Im(n,m)])/(2.0*dx);
      }
    }  
}    
    



void main(int argc, char** argv) 
{
    double dxQ[(FMM_P_MAX+2)*(FMM_P_MAX+4)+1];
    double dyQ[(FMM_P_MAX+2)*(FMM_P_MAX+4)+1];
    double dzQ[(FMM_P_MAX+2)*(FMM_P_MAX+4)+1];
    double dxQ_nd[(FMM_P_MAX+2)*(FMM_P_MAX+4)+1];
    double dyQ_nd[(FMM_P_MAX+2)*(FMM_P_MAX+4)+1];
    double dzQ_nd[(FMM_P_MAX+2)*(FMM_P_MAX+4)+1];
    double x,y,z, kappa, dx;
    int pp,p,n,m;

    pp = -3;
    x = .23;
    y = -.17;
    z = .29;
    dx = 1e-4;
    kappa = 1;

    p = abs(pp);
    init_coeffs();

    
    compute_dQ_nd(pp, dxQ_nd, dyQ_nd, dzQ_nd, dx, kappa, x, y, z);
    compute_dQ(pp, dxQ, dyQ, dzQ, kappa, x, y, z);

    printf("*** dx: *********************************\n");
    for (n=0; n<=p; n++) {
    for (m=0; m<=n; m++) {
        printf("%3i %3i  %18.10e %18.10e    %18.10e %18.10e\n", n,m, dxQ[Re(n,m)], dxQ_nd[Re(n,m)], dxQ[Im(n,m)], dxQ_nd[Im(n,m)]);
    }
    }
    printf("*** dy: *********************************\n");
    for (n=0; n<=p; n++) {
    for (m=0; m<=n; m++) {
        printf("%3i %3i  %18.10e %18.10e    %18.10e %18.10e\n", n,m, dyQ[Re(n,m)], dyQ_nd[Re(n,m)], dyQ[Im(n,m)], dyQ_nd[Im(n,m)]);
    }
    }
    printf("*** dz: *********************************\n");
    for (n=0; n<=p; n++) {
    for (m=0; m<=n; m++) {
        printf("%3i %3i  %18.10e %18.10e    %18.10e %18.10e\n", n,m, dzQ[Re(n,m)], dzQ_nd[Re(n,m)], dzQ[Im(n,m)], dzQ_nd[Im(n,m)]);
    }
    }

}
