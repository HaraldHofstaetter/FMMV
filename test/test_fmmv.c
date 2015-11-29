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

#if (FMM_DIM==2)
#include"fmmv2d.h"
#else
#include"fmmv3d.h"
#endif

#include<assert.h>
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
//#include<fenv.h> 
//#include<sys/time.h>
//#include<sys/resource.h>
//#include<pthread.h>


#undef USE_PAPI

#ifdef USE_PAPI
  #include<papi.h>
#endif
#ifdef USE_SINGLE_PRECISION
typedef float _FLOAT_;
#else
typedef double _FLOAT_;
#endif
#include"test_utilities.h"

int main(int argc, char *argv[]) 
{
  struct FmmvOptions options = fmmvGetDefaultOptions();
  struct FmmvStatistics statistics;
  int NParticles=10000;
  int NTargets=0;  // if NTargets==0: targets == sources
  int with_dipoleSources=0;
  int with_gradients=0;
  int with_targets = 0;

  _FLOAT_ (*particles)[FMM_DIM];
  _FLOAT_ (*targets)[FMM_DIM] = NULL;
  _FLOAT_ *charges;
  _FLOAT_ (*dipoleMoments)[FMM_DIM] = NULL;
  _FLOAT_ x, y, z, q, mx, my, mz;
  _FLOAT_ *pot;
  _FLOAT_ (*grad)[FMM_DIM] = NULL;
  _FLOAT_ *exPot;
  _FLOAT_ (*exGrad)[FMM_DIM] = NULL;
  int i;
  char *err_message;
  _FLOAT_ pi = 4.0*atan(1.0);
  _FLOAT_ phi, theta;
  _FLOAT_ errL2, errL2grad;
  double exTime;
  int exAcc;
// #include<fenv.h> 

 // struct rusage usage;
  
  int NDirect = 100;
  int printResult = 0;
  int directEvalExtendedPrecision = 0;
  int noOfExample = 0;
  int positiveCharges = 0;

//  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );

  options.splitThreshold=190;
  options.levels=51;
  options.directEvalThreshold=-1;
  options.periodicBoundaryConditions = 0;
  options.extrinsicCorrection = 0;
  options.pM = 6;
  options.pL = 6;
  options.s = 8;
  options.ws = 1;
  options.reducedScheme = 0;
  options.useHilbertOrder = 0;
#ifdef USE_SINGLE_PRECISION
  options.directEvalAccuracy = 1;
  exAcc = 1;
#else
  options.directEvalAccuracy = 2;
  exAcc = 2;
#endif
  options.scale = 1.0;
  options.useFarfieldNearfieldThreads = 0;
  options.beta = 0.0;

  set_from_command_line_int(argc, argv, "example", &noOfExample);
  set_from_command_line_int(argc, argv, "positiveCharges", &positiveCharges);

  set_from_command_line_bool(argc, argv, "dipoleSources", &with_dipoleSources);
  set_from_command_line_bool(argc, argv, "gradients", &with_gradients);
  set_from_command_line_bool(argc, argv, "printResult", &printResult);
  set_from_command_line_bool(argc, argv, "directEvalExtendedPrecision", &directEvalExtendedPrecision);

  set_from_command_line_int(argc, argv, "NParticles", &NParticles);
  set_from_command_line_int(argc, argv, "NTargets", &NTargets);
  set_from_command_line_int(argc, argv, "NDirect", &NDirect);

  set_from_command_line_int(argc, argv, "pM", &options.pM);
  set_from_command_line_int(argc, argv, "pL", &options.pL);
  set_from_command_line_int(argc, argv, "s", &options.s);
    
  set_from_command_line_int(argc, argv, "ws", &options.ws);
  set_from_command_line_bool(argc, argv, "reducedScheme", &options.reducedScheme);
  
  set_from_command_line_double(argc, argv, "scale", &options.scale);
  set_from_command_line_double(argc, argv, "beta", &options.beta);
  set_from_command_line_int(argc, argv, "splitThreshold", &options.splitThreshold);
  set_from_command_line_int(argc, argv, "splitTargetThreshold", &options.splitTargetThreshold);
  set_from_command_line_int(argc, argv, "levels", &options.levels);
  set_from_command_line_int(argc, argv, "directEvalThreshold", &options.directEvalThreshold);
  set_from_command_line_bool(argc, argv, "periodicBoundaryConditions", &options.periodicBoundaryConditions);
  set_from_command_line_bool(argc, argv, "extrinsicCorrection", &options.extrinsicCorrection);
  set_from_command_line_bool(argc, argv, "useHilbertOrder", &options.useHilbertOrder);
  set_from_command_line_int(argc, argv, "directEvalAccuracy", &options.directEvalAccuracy);
  set_from_command_line_bool(argc, argv, "useFarfieldNearfieldThreads", &options.useFarfieldNearfieldThreads);

  
  with_targets = (NTargets>0);
  if (!with_targets) {
       NTargets = NParticles;
  }     

  if (NTargets < NDirect) {
    NDirect = NTargets;
  }	

  //options.PAPIeventSet = PAPI_NULL;
#ifdef USE_PAPI
  if (strstr(argv[1],"papi")) {
	int retval;	
	retval = PAPI_library_init(PAPI_VER_CURRENT);
	assert(retval == PAPI_VER_CURRENT);
	//retval = PAPI_thread_init(pthread_self);
	//assert(retval == PAPI_OK);
		
	options.PAPIeventSet = PAPI_NULL;	
	retval = PAPI_create_eventset(&options.PAPIeventSet);
		if (retval != PAPI_OK) 
			fprintf(stderr, "PAPI error %d: %s\n",retval,PAPI_strerror(retval));
		assert(retval == PAPI_OK);
	
//	retval = PAPI_add_event(options.PAPIeventSet, PAPI_TOT_INS);
//	retval = PAPI_add_event(options.PAPIeventSet, PAPI_TOT_CYC);
//	retval = PAPI_add_event(options.PAPIeventSet, PAPI_VEC_INS);
//	retval = PAPI_add_event(options.PAPIeventSet, PAPI_FP_INS);
	retval = PAPI_add_event(options.PAPIeventSet, PAPI_FP_OPS);
//	retval = PAPI_add_event(options.PAPIeventSet, PAPI_L1_DCM);
//	retval = PAPI_add_event(options.PAPIeventSet, PAPI_L2_DCM);
	/*
		if (retval != PAPI_OK) 
			fprintf(stderr, "PAPI error %d: %s\n",retval,PAPI_strerror(retval));
		assert(retval == PAPI_OK);
	*/	
	/* add more events here ... */	
  }	
#endif
  
  particles = (_FLOAT_ (*)[FMM_DIM]) calloc(NParticles, FMM_DIM*sizeof(_FLOAT_));  
  assert(particles!=NULL);  
  charges = (_FLOAT_ *) calloc(NParticles, sizeof(_FLOAT_));  
  assert(charges!=NULL);
  pot = (_FLOAT_ *) calloc(NTargets, sizeof(_FLOAT_));
  assert(pot!=NULL);
  if (with_gradients) {
	grad = (_FLOAT_ (*)[FMM_DIM]) calloc(NTargets, FMM_DIM*sizeof(_FLOAT_));  
	assert(grad!=NULL);
  }	  

  my_srand(1481765933);
  //my_srand(1111111111);

  for (i = 0; i < NParticles; i++) {
  	switch(noOfExample) {
	case 0:
	      x = 0.01 + 0.98*(my_rand() / (_FLOAT_) MY_RAND_MAX);
	      y = 0.01 + 0.98*(my_rand() / (_FLOAT_) MY_RAND_MAX);
#if (FMM_DIM==3)
	      z = 0.01 + 0.98*(my_rand() / (_FLOAT_) MY_RAND_MAX);
#endif
	      break;

	case 1:      
      	      phi = 2.0*pi*my_rand() / (_FLOAT_) MY_RAND_MAX;
              theta = pi*my_rand() / (_FLOAT_) MY_RAND_MAX;
              x = 0.5 + 0.49 * sin(theta) * cos(phi);
              y = 0.5 + 0.49 * sin(theta) * sin(phi);
#if (FMM_DIM==3)
              z = 0.5 + 0.49 * cos(theta);
#endif
	      break;
	case 2:
	      { _FLOAT_ r, r2;
	      do {	
	      x = 2.0*(my_rand() / (_FLOAT_) MY_RAND_MAX) - 1.0;
	      y = 2.0*(my_rand() / (_FLOAT_) MY_RAND_MAX) - 1.0;
#if (FMM_DIM==3)
	      z = 2.0*(my_rand() / (_FLOAT_) MY_RAND_MAX)  - 1.0;
#endif
	      r2 = x*x+y*y+z*z;
	      } while (r2>1.0);
	      r = 0.49*sqrt(r2);
	      x = 0.5 + r*x;
	      y = 0.5 + r*y;
#if (FMM_DIM==3)
	      z = 0.5 + r*z;
#endif
	      }
	}
	

	if (positiveCharges) {
              q = ((my_rand() / (_FLOAT_) MY_RAND_MAX));
	}
	else {
              q = 2.0*((my_rand() / (_FLOAT_) MY_RAND_MAX)) - 1.0 ;
	}

      particles[i][0] = x;
      particles[i][1] = y;
#if (FMM_DIM==3)
      particles[i][2] = z;
#endif
      //charges[i] = q;
      charges[i] = 1.0;
  }
  if (with_dipoleSources) {
  	dipoleMoments = (_FLOAT_ (*)[FMM_DIM]) calloc(NParticles, FMM_DIM*sizeof(_FLOAT_));  
  	assert(dipoleMoments!=NULL);

  	for (i = 0; i < NParticles; i++) {
	      mx = (my_rand() / (_FLOAT_) MY_RAND_MAX) - 0.5 ;
	      my = (my_rand() / (_FLOAT_) MY_RAND_MAX) - 0.5 ;
#if (FMM_DIM==3)
	      mz = (my_rand() / (_FLOAT_) MY_RAND_MAX) - 0.5 ;
#endif

	      dipoleMoments[i][0] = mx;
	      dipoleMoments[i][1] = my;
#if (FMM_DIM==3)
	      dipoleMoments[i][2] = mz;
#endif
	}
  }
  if (with_targets) {
     targets = (_FLOAT_ (*)[FMM_DIM]) calloc(NTargets, FMM_DIM*sizeof(_FLOAT_));  
     assert(targets!=NULL);  
     for (i = 0; i < NTargets; i++) {
  	switch(noOfExample) {
	case 0:
	      x = 0.01 + 0.98*(my_rand() / (_FLOAT_) MY_RAND_MAX);
	      y = 0.01 + 0.98*(my_rand() / (_FLOAT_) MY_RAND_MAX);
#if (FMM_DIM==3)
	      z = 0.01 + 0.98*(my_rand() / (_FLOAT_) MY_RAND_MAX);
#endif
	      break;

	case 1:      
      	      phi = 2.0*pi*my_rand() / (_FLOAT_) MY_RAND_MAX;
              theta = pi*my_rand() / (_FLOAT_) MY_RAND_MAX;
              x = 0.5 + 0.49 * sin(theta) * cos(phi);
              y = 0.5 + 0.49 * sin(theta) * sin(phi);
#if (FMM_DIM==3)
              z = 0.5 + 0.49 * cos(theta);
#endif
	      break;
	case 2:
	      { _FLOAT_ r, r2;
	      do {	
	      x = 2.0*(my_rand() / (_FLOAT_) MY_RAND_MAX) - 1.0;
	      y = 2.0*(my_rand() / (_FLOAT_) MY_RAND_MAX) - 1.0;
#if (FMM_DIM==3)
	      z = 2.0*(my_rand() / (_FLOAT_) MY_RAND_MAX)  - 1.0;
#endif
	      r2 = x*x+y*y+z*z;
	      } while (r2>1.0);
	      r = 0.49*sqrt(r2);
	      x = 0.5 + r*x;
	      y = 0.5 + r*y;
#if (FMM_DIM==3)
	      z = 0.5 + r*z;
#endif
	      }
	}
        targets[i][0] = x;
        targets[i][1] = y;
#if (FMM_DIM==3)
        targets[i][2] = z;
#endif
     }
  }   

/*  
  getrusage(RUSAGE_SELF, &usage);
  printf("Data from getrusage():\n");
//  printf("  ru_utime    (user time used) \n");
//  printf("  ru_stime    (system time used)\n");
  printf("  ru_maxrss   (maximum resident set size)    : %i\n", usage.ru_maxrss);
  printf("  ru_ixrss    (integral shared memory size)  : %i\n", usage.ru_ixrss);
  printf("  ru_idrss    (integral unshared data size)  : %i\n", usage.ru_idrss);
  printf("  ru_isrss    (integral unshared stack size) : %i\n", usage.ru_isrss);
  printf("  ru_minflt   (page reclaims)                : %i\n", usage.ru_minflt);
  printf("  ru_majflt   (page faults)                  : %i\n", usage.ru_majflt);
  printf("  ru_nswap    (swaps)                        : %i\n", usage.ru_nswap);
  printf("  ru_inblock  (block input operations)       : %i\n", usage.ru_inblock);
  printf("  ru_oublock  (block output operations)      : %i\n", usage.ru_oublock);
  printf("  ru_msgsnd   (messages sent)                : %i\n", usage.ru_msgsnd);
  printf("  ru_msgrcv   (messages received)            : %i\n", usage.ru_msgrcv);
  printf("  ru_nsignals (ignals received)              : %i\n", usage.ru_nsignals);
  printf("  ru_nvcsw    (voluntary context switches)   : %i\n", usage.ru_nvcsw);
  printf("  ru_nivcsw   (involuntary context switches) : %i\n", usage.ru_nivcsw);
  
*/

  printf("\nFMM STARTED.\n");
#if (FMM_DIM==2)
#ifdef USE_SINGLE_PRECISION
      fmmv2df
#else	  
      fmmv2d
#endif	  
#else
#ifdef USE_SINGLE_PRECISION
      fmmv3df
#else	  
      fmmv3d
#endif	  
#endif	  
  	(NParticles, particles, charges, dipoleMoments,
	NTargets, targets, 
	pot, grad, &options, &statistics, &err_message);
  if (err_message) {
      	printf("%s\n", err_message);
	exit(EXIT_FAILURE);
  }	
  printf("FMM FINISHED.\n\n");
  printFmmvStatistics(&statistics);

  exPot = (_FLOAT_ *) calloc(NDirect, sizeof(_FLOAT_));
  assert(exPot!=NULL);
  if (with_gradients) {
	exGrad = (_FLOAT_ (*)[FMM_DIM]) calloc(NDirect, FMM_DIM*sizeof(_FLOAT_));  
	assert(exGrad!=NULL);
  }	  

  printf("\nDIRECT METHOD STARTED.\n");
#if (FMM_DIM==2)
#ifdef USE_SINGLE_PRECISION
    fmmv2df_direct
#else
    fmmv2d_direct
#endif
#else
#ifdef USE_SINGLE_PRECISION
    fmmv3df_direct
#else
    fmmv3d_direct
#endif
#endif
    ( NParticles, particles, charges, dipoleMoments,
        NDirect, targets, 
        exPot, exGrad, exAcc, 
        options.beta, 
        &exTime, &err_message);
printf("DIRECT METHOD FINISHED.\n");

/*
  getrusage(RUSAGE_SELF, &usage);
  printf("Data from getrusage():\n");
//  printf("  ru_utime    (user time used) \n");
//  printf("  ru_stime    (system time used)\n");
  printf("  ru_maxrss   (maximum resident set size)    : %i\n", usage.ru_maxrss);
  printf("  ru_ixrss    (integral shared memory size)  : %i\n", usage.ru_ixrss);
  printf("  ru_idrss    (integral unshared data size)  : %i\n", usage.ru_idrss);
  printf("  ru_isrss    (integral unshared stack size) : %i\n", usage.ru_isrss);
  printf("  ru_minflt   (page reclaims)                : %i\n", usage.ru_minflt);
  printf("  ru_majflt   (page faults)                  : %i\n", usage.ru_majflt);
  printf("  ru_nswap    (swaps)                        : %i\n", usage.ru_nswap);
  printf("  ru_inblock  (block input operations)       : %i\n", usage.ru_inblock);
  printf("  ru_oublock  (block output operations)      : %i\n", usage.ru_oublock);
  printf("  ru_msgsnd   (messages sent)                : %i\n", usage.ru_msgsnd);
  printf("  ru_msgrcv   (messages received)            : %i\n", usage.ru_msgrcv);
  printf("  ru_nsignals (ignals received)              : %i\n", usage.ru_nsignals);
  printf("  ru_nvcsw    (voluntary context switches)   : %i\n", usage.ru_nvcsw);
  printf("  ru_nivcsw   (involuntary context switches) : %i\n", usage.ru_nivcsw);
*/
  

  if (with_gradients) {
	  relL2Err(0, 0, 0); /* init */
#if (FMM_DIM==2)
	  relL2Err2(0, NULL, NULL); /* init */
#else          
	  relL2Err3(0, NULL, NULL); /* init */
#endif          

	  if (printResult) {
  	  	printf("\n        Pot(FMM)   Pot(exact)     rel.err.  rel.err.(grad)\n");
	  	printf("============================================================\n");
	  }	
	  for (i=0; i<NDirect; i++) {
		relL2Err(1, pot[i], exPot[i]); /* accumulate */
#if (FMM_DIM==2)
		relL2Err2(1, grad[i], exGrad[i]); /* accumulate */
#else          
		relL2Err3(1, grad[i], exGrad[i]); /* accumulate */
#endif          
	  	if (printResult) {
			printf("%3i %12.4e %12.4e %12.4e %12.4e\n", i, 
                               (double) pot[i], (double) exPot[i], (double) relErr(pot[i], 
                               exPot[i]),
#if (FMM_DIM==2)
                               (_FLOAT_) relErr3(grad[i], exGrad[i])
#else          
                               (_FLOAT_) relErr3(grad[i], exGrad[i])
#endif          
                               );
		}	
  	}
	errL2 = relL2Err(2, 0, 0); /*finalize*/ 
#if (FMM_DIM==2)
	errL2grad = relL2Err2(2, NULL, NULL); /*finalize*/ 
#else          
	errL2grad = relL2Err3(2, NULL, NULL); /*finalize*/ 
#endif          
	printf ("\nerr_L2(pot) = %.4e  err_L2(grad) = %.4e\n", errL2, errL2grad);
  }	  
  else {	  
	  relL2Err(0, 0, 0); /* init */
	  if (printResult) {
	  	printf("\n        Pot(FMM)   Pot(exact)     rel.err.\n");
	  	printf("============================================\n");
	  }	
	  for (i=0; i<NDirect; i++) {
		relL2Err(1, pot[i], exPot[i]); /* accumulate */
	  	if (printResult) {
			printf("%3i %24.16e %24.16e %12.4e\n", i, 
                        (double) pot[i], (double) exPot[i], (double) relErr(pot[i], exPot[i]));
		}	
  	}
	errL2 = relL2Err(2, 0, 0); /*finalize*/ 
	printf ("\nerr_L2(pot) = %.4e\n", errL2);
  }

  printf("\n");

  
  free(particles);
  free(charges);
  free(dipoleMoments);
  free(targets);
  free(pot);
  free(exPot);
  free(grad);
  free(exGrad);
  return 0;
}	





