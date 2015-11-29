#ifdef USE_SINGLE_PRECISION 
#define FMM_PRECISION 0
#else
#define FMM_PRECISION 1
#endif
#include "_fmmv.h"
#include<stdio.h>
#include<math.h>
#include<time.h>
#include"test_utilities.h"

void init_M2M(FmmvHandle *FMMV, int level);
void init_L2L(FmmvHandle *FMMV, int level);
void finish_M2M(FmmvHandle *FMMV);
void finish_L2L(FmmvHandle *FMMV);


int main(int argc, char**argv)
{
    int err;
    int i;
    double t0, t1;
    FmmvHandle _FMMV;
    FmmvHandle *FMMV = &_FMMV;
    void (*GEN_M)(FmmvHandle *FMMV, Box *box);
    void (*EVAL_M)(FmmvHandle *FMMV, Box *target, Box *source);
    void (*EVAL_L)(FmmvHandle *FMMV, Box *box);
    void (*GEN_L)(FmmvHandle *FMMV, Box *target, Box *source);
    void (*GEN_L_EVAL_M)(FmmvHandle *FMMV, Box *target, Box *source);
    void (*EVAL_DIRECT)(FmmvHandle *FMMV, Box *target, Box *source) = FMMV->eval_direct;

    Box _box_A;
    Box _box_B;

    Box *box_A = &_box_A;
    Box *box_B = &_box_B;

    Box _child_A[8];
    Box _child_B[8];

    int NParticles = 100;
    int NTargets = 100;
     _FLOAT_ (*particles)[3];
     _FLOAT_ (*targets)[3] = NULL;
     _FLOAT_ *charges;
     _FLOAT_ (*dipoleMoments)[3] = NULL;
     _FLOAT_ x, y, z;
     _FLOAT_ *potM;
     _FLOAT_ *potL;
     _FLOAT_ (*gradM)[3] = NULL;
     _FLOAT_ (*gradL)[3] = NULL;
     _FLOAT_ *exPot;
     _FLOAT_ (*exGrad)[3] = NULL;
     _FLOAT_ errL2, errL2grad;
     double beta = 0.0;
     int pM = 8;
     int pL = 8;
     int s = 8;
     int *perm;
     int *permTargets;
     int with_dipoleSources=0;
     int with_gradients=0;
     int LM_combined=0;
     int with_M2M=0;
     int with_L2L=0;
     int printResult=0;
#ifdef USE_SINGLE_PRECISION 
     int exAcc=1;
#else     
     int exAcc=2;
#endif
     //int exAcc=0;

    _FLOAT_ x_A = .25;
    _FLOAT_ y_A = .25;
    _FLOAT_ z_A = .25;
    double size_A = .125;
    int level_A = 2;
    _FLOAT_ x_B = .75;
    _FLOAT_ y_B = .75;
    _FLOAT_ z_B = .75;
    double size_B = .125;
    int level_B = 2;

    set_from_command_line_double(argc, argv, "beta", &beta);
    set_from_command_line_bool(argc, argv, "dipoleSources", &with_dipoleSources);
    set_from_command_line_bool(argc, argv, "gradients", &with_gradients);
    set_from_command_line_bool(argc, argv, "printResult", &printResult);
    set_from_command_line_bool(argc, argv, "LM_combined", &LM_combined);
    set_from_command_line_bool(argc, argv, "M2M", &with_M2M);
    set_from_command_line_bool(argc, argv, "L2L", &with_L2L);
    set_from_command_line_int(argc, argv, "directEvalAccuracy", &exAcc);
 
        if (with_dipoleSources) {
            if (with_gradients) {
                FMMV->dataKind = FMM_ST_DIPOLE_GRAD;
                GEN_M = gen_M_ST_dipole_grad;
                EVAL_M = eval_M_ST_dipole_grad;
                GEN_L = gen_L_ST_dipole_grad;
                EVAL_L = eval_L_ST_dipole_grad;
                GEN_L_EVAL_M = gen_L_eval_M_dipole_grad;
                switch (exAcc) {
                    case 0: EVAL_DIRECT = eval_direct_ST_dipole_grad_acc0; break;
                    case 1: EVAL_DIRECT = eval_direct_ST_dipole_grad_acc1; break;
#if FMM_PRECISION==1
                    case 2: EVAL_DIRECT = eval_direct_ST_dipole_grad_acc2; break;
#endif                
                }    
            }
            else {
                FMMV->dataKind = FMM_ST_DIPOLE;
                GEN_M = gen_M_ST_dipole;
                EVAL_M = eval_M_ST_dipole;
                GEN_L = gen_L_ST_dipole;
                EVAL_L = eval_L_ST_dipole;
                GEN_L_EVAL_M = gen_L_eval_M_dipole;
                switch (exAcc) {
                    case 0: EVAL_DIRECT = eval_direct_ST_dipole_acc0; break;
                    case 1: EVAL_DIRECT = eval_direct_ST_dipole_acc1; break;
#if FMM_PRECISION==1
                    case 2: EVAL_DIRECT = eval_direct_ST_dipole_acc2; break;
#endif                
                }    
            }
        }
        else {
            if (with_gradients) {
                FMMV->dataKind = FMM_ST_GRAD;
                GEN_M = gen_M_ST_grad;
                EVAL_M = eval_M_ST_grad;
                GEN_L = gen_L_ST_grad;
                EVAL_L = eval_L_ST_grad;
                GEN_L_EVAL_M = gen_L_eval_M_grad;
                switch (exAcc) {
                    case 0: EVAL_DIRECT = eval_direct_ST_grad_acc0; break;
                    case 1: EVAL_DIRECT = eval_direct_ST_grad_acc1; break;
#if FMM_PRECISION==1
                    case 2: EVAL_DIRECT = eval_direct_ST_grad_acc2; break;
#endif                
                }    
            }
            else {
                FMMV->dataKind = FMM_ST_STANDARD;
                GEN_M = gen_M_ST_standard;
                EVAL_M = eval_M_ST_standard;
                GEN_L = gen_L_ST_standard;
                EVAL_L = eval_L_ST_standard;
                GEN_L_EVAL_M = gen_L_eval_M_standard;
                switch (exAcc) {
                    case 0: EVAL_DIRECT = eval_direct_ST_standard_acc0; break;
                    case 1: EVAL_DIRECT = eval_direct_ST_standard_acc1; break;
#if FMM_PRECISION==1
                    case 2: EVAL_DIRECT = eval_direct_ST_standard_acc2; break;
#endif                
                }    
            }
        }

    set_from_command_line_int(argc, argv, "NParticles", &NParticles);
    set_from_command_line_int(argc, argv, "NTargets", &NTargets);
    set_from_command_line_int(argc, argv, "pM", &pM);
    set_from_command_line_int(argc, argv, "pL", &pL);
    set_from_command_line_int(argc, argv, "s", &s);
    set_from_command_line_int(argc, argv, "level_A", &level_A);
    set_from_command_line_int(argc, argv, "level_B", &level_B);
//    set_from_command_line_double(argc, argv, "size_A", &size_A);
//    set_from_command_line_double(argc, argv, "size_B", &size_B);
    size_A = ldexp(1.0, -level_A);
    size_B = ldexp(1.0, -level_B);
    
    particles = (_FLOAT_ (*)[3]) calloc(NParticles, 3*sizeof(_FLOAT_));  
    charges = (_FLOAT_ *) calloc(NParticles, sizeof(_FLOAT_));  
    perm = (int *) calloc(NParticles, sizeof(int));  
    permTargets = (int *) calloc(NTargets, sizeof(int));  
    targets = (_FLOAT_ (*)[3]) calloc(NTargets, 3*sizeof(_FLOAT_));  
    potM = (_FLOAT_ *) calloc(NTargets, sizeof(_FLOAT_));
    potL = (_FLOAT_ *) calloc(NTargets, sizeof(_FLOAT_));
    exPot = (_FLOAT_ *) calloc(NTargets, sizeof(_FLOAT_));
    if (with_gradients) {
	gradM = (_FLOAT_ (*)[3]) calloc(NTargets, 3*sizeof(_FLOAT_));  
	gradL = (_FLOAT_ (*)[3]) calloc(NTargets, 3*sizeof(_FLOAT_));  
	exGrad = (_FLOAT_ (*)[3]) calloc(NTargets, 3*sizeof(_FLOAT_));  
    }	   
    if (with_dipoleSources) {
	dipoleMoments = (_FLOAT_ (*)[3]) calloc(NParticles, 3*sizeof(_FLOAT_));  
    }

    FMMV->beta = beta;
    FMMV->pM = pM;
    FMMV->pL = pL;
    FMMV->s_eps = s;
    FMMV->particles = particles;
    FMMV->charges = charges;
    FMMV->dipoleMoments = dipoleMoments;
    FMMV->perm = perm;
    FMMV->NParticles = NParticles;
    FMMV->targets = targets;
    FMMV->permTargets = permTargets;
    FMMV->NTargets = NTargets;
    FMMV->scale = 1.0;
    FMMV->lambda = 0;

    my_srand(1481765933);

    box_A->firstParticle = 0;
    box_A->noOfParticles = NParticles;
    box_A->x = x_A;
    box_A->y = y_A;
    box_A->z = z_A;
    box_A->level = level_A;
    box_A->M = 0; 
    /* if (with_M2M) */ {
        int m = NParticles/8;
        double d = 0.25*size_A;
        int c;
        for (c=0; c<8; c++) {
            box_A->child[c] = &_child_A[c];
            box_A->child[c]->level = level_A + 1;
            box_A->child[c]->M = 0;
            box_A->child[c]->firstParticle = c*m;
            box_A->child[c]->noOfParticles = m;
        }
        box_A->child[7]->noOfParticles = NParticles - box_A->child[7]->firstParticle;
         
        box_A->child[SWD]->x = x_A - d;
        box_A->child[SWD]->y = y_A - d;
        box_A->child[SWD]->z = z_A - d;

        box_A->child[NWD]->x = x_A - d;
        box_A->child[NWD]->y = y_A + d;
        box_A->child[NWD]->z = z_A - d;

        box_A->child[SED]->x = x_A + d;
        box_A->child[SED]->y = y_A - d;
        box_A->child[SED]->z = z_A - d;

        box_A->child[NED]->x = x_A + d;
        box_A->child[NED]->y = y_A + d;
        box_A->child[NED]->z = z_A - d;

        box_A->child[SWU]->x = x_A - d;
        box_A->child[SWU]->y = y_A - d;
        box_A->child[SWU]->z = z_A + d;

        box_A->child[NWU]->x = x_A - d;
        box_A->child[NWU]->y = y_A + d;
        box_A->child[NWU]->z = z_A + d;

        box_A->child[SEU]->x = x_A + d;
        box_A->child[SEU]->y = y_A - d;
        box_A->child[SEU]->z = z_A + d;

        box_A->child[NEU]->x = x_A + d;
        box_A->child[NEU]->y = y_A + d;
        box_A->child[NEU]->z = z_A + d;
        for (c=0; c<8; c++) {
            for (i = box_A->child[c]->firstParticle; i < box_A->child[c]->firstParticle+box_A->child[c]->noOfParticles; i++) {
                x = box_A->child[c]->x + 2.0*d*(my_rand() / (_FLOAT_) MY_RAND_MAX) - d;
                y = box_A->child[c]->y + 2.0*d*(my_rand() / (_FLOAT_) MY_RAND_MAX) - d;
                z = box_A->child[c]->z + 2.0*d*(my_rand() / (_FLOAT_) MY_RAND_MAX) - d;
                particles[i][0] = x;
                particles[i][1] = y;
                particles[i][2] = z;
            }
        }
    }
    /* else {
        double d = size_A/2.0;
        for (i=0;i<NParticles;i++) {
            x = x_A + size_A*(my_rand() / (_FLOAT_) MY_RAND_MAX) - d;
            y = y_A + size_A*(my_rand() / (_FLOAT_) MY_RAND_MAX) - d;
            z = y_A + size_A*(my_rand() / (_FLOAT_) MY_RAND_MAX) - d;
            particles[i][0] = x;
            particles[i][1] = y;
            particles[i][2] = z;
        }    
    }*/
    for (i=0;i<NParticles;i++) {
        perm[i] = i;
        charges[i] = 1.0;
    }
    if (with_dipoleSources) {
        for (i=0;i<NParticles;i++) {
            x = (my_rand() / (_FLOAT_) MY_RAND_MAX) - 0.5;
            y = (my_rand() / (_FLOAT_) MY_RAND_MAX) - 0.5;
            z = (my_rand() / (_FLOAT_) MY_RAND_MAX) - 0.5;
            dipoleMoments[i][0] = x;
            dipoleMoments[i][1] = y;
            dipoleMoments[i][2] = z;

        }
    }

    box_B->firstTarget = 0;
    box_B->noOfTargets = NTargets;
    box_B->x = x_B;
    box_B->y = y_B;
    box_B->z = z_B;
    box_B->level = level_B;
    box_B->L = 0; 
    /* if (with_L2L)  */ { 
        int m = NTargets/8;
        double d = 0.25*size_B;
        int c;
        for (c=0; c<8; c++) {
            box_B->child[c] = &_child_B[c];
            box_B->child[c]->level = level_B + 1;
            box_B->child[c]->L = 0;
            box_B->child[c]->firstTarget = c*m;
            box_B->child[c]->noOfTargets = m;
        }
        box_B->child[7]->noOfTargets = NTargets - box_B->child[7]->firstTarget;
         
        box_B->child[SWD]->x = x_B - d;
        box_B->child[SWD]->y = y_B - d;
        box_B->child[SWD]->z = z_B - d;

        box_B->child[NWD]->x = x_B - d;
        box_B->child[NWD]->y = y_B + d;
        box_B->child[NWD]->z = z_B - d;

        box_B->child[SED]->x = x_B + d;
        box_B->child[SED]->y = y_B - d;
        box_B->child[SED]->z = z_B - d;

        box_B->child[NED]->x = x_B + d;
        box_B->child[NED]->y = y_B + d;
        box_B->child[NED]->z = z_B - d;

        box_B->child[SWU]->x = x_B - d;
        box_B->child[SWU]->y = y_B - d;
        box_B->child[SWU]->z = z_B + d;

        box_B->child[NWU]->x = x_B - d;
        box_B->child[NWU]->y = y_B + d;
        box_B->child[NWU]->z = z_B + d;

        box_B->child[SEU]->x = x_B + d;
        box_B->child[SEU]->y = y_B - d;
        box_B->child[SEU]->z = z_B + d;

        box_B->child[NEU]->x = x_B + d;
        box_B->child[NEU]->y = y_B + d;
        box_B->child[NEU]->z = z_B + d;
        for (c=0; c<8; c++) {
            for (i = box_B->child[c]->firstTarget; i < box_B->child[c]->firstTarget+box_B->child[c]->noOfTargets; i++) {
                x = box_B->child[c]->x + 2.0*d*(my_rand() / (_FLOAT_) MY_RAND_MAX) - d;
                y = box_B->child[c]->y + 2.0*d*(my_rand() / (_FLOAT_) MY_RAND_MAX) - d;
                z = box_B->child[c]->z + 2.0*d*(my_rand() / (_FLOAT_) MY_RAND_MAX) - d;
                targets[i][0] = x;
                targets[i][1] = y;
                targets[i][2] = z;
            }
        }
    }
    /* else {
        double d = size_B/2.0;
        for (i=0;i<NTargets;i++) {
            x = x_B + size_B*(my_rand() / (_FLOAT_) MY_RAND_MAX) - d;
            y = y_B + size_B*(my_rand() / (_FLOAT_) MY_RAND_MAX) - d;
            z = y_B + size_B*(my_rand() / (_FLOAT_) MY_RAND_MAX) - d;
            targets[i][0] = x;
            targets[i][1] = y;
            targets[i][2] = z;
        }    
    } */
    for (i=0;i<NTargets;i++) {
        permTargets[i] = i;
    }
    


    init_all(FMMV);

    err = ida_allocate(FMMV);
    copy_particles(FMMV);
    copy_charges(FMMV);
    
    zero_pot(FMMV);
    t0 = (double) clock()/CLOCKS_PER_SEC;
    EVAL_DIRECT(FMMV, box_B, box_A);
    t1 = (double) clock()/CLOCKS_PER_SEC - t0;
    printf("time (eval_direct)  : %7.4f\n", t1);
    FMMV->potentials = exPot;
    FMMV->gradients = exGrad;
    backcopy_pot(FMMV);

    if (with_M2M) {
        int c;
        t0 = (double) clock()/CLOCKS_PER_SEC;
        for (c=0; c<8; c++) {
            GEN_M(FMMV, box_A->child[c]);    
        }
        t1 = (double) clock()/CLOCKS_PER_SEC - t0;
        printf("time (GEN_M)        : %7.4f\n", t1);
        init_M2M(FMMV, -1);
        init_M2M(FMMV, level_A);
        t0 = (double) clock()/CLOCKS_PER_SEC;
        M2M(FMMV, box_A);
        t1 = (double) clock()/CLOCKS_PER_SEC - t0;
        printf("time (M2M)          : %7.4f\n", t1);
    }
    else {
        t0 = (double) clock()/CLOCKS_PER_SEC;
        GEN_M(FMMV, box_A);
        t1 = (double) clock()/CLOCKS_PER_SEC - t0;
        printf("time (GEN_M)        : %7.4f\n", t1);
    }    
    		
#if 0
                { int m,n, ii;
                  double h;
                  double xxx[] = { 1.0, 3.0, 15.0, 105.0, 945.0, 10395.0, 135135.0, 2027025.0,
	                           34459425.0, 654729075.0, 1.374931057500000e+10};

                ii = 0;
                for (n=0;n<=pM;n++) {
                    for (m=0;m<=n;m++) {
                        if(FMMV->beta==0) {
                            h = 1.0;
                        }
                        else 
                        {
                            h = xxx[n]*pow(FMMV->beta,-(n+1));
                        }
                        printf(" %i %i %18.10e %18.10e\n", n, m, h*box_A->M[ii], h*box_A->M[ii+1]);
                        ii+=2;
                    }
                }    
                
                }
#endif
/*    
    if(LM_combined) {
        zero_pot(FMMV);
        t0 = (double) clock()/CLOCKS_PER_SEC;
        GEN_L_EVAL_M(FMMV, box_B, box_A);
        t1 = (double) clock()/CLOCKS_PER_SEC - t0;
        printf("time (GEN_L_EVAL_M): %7.4f\n", t1);
        FMMV->potentials = potM;
        backcopy_pot(FMMV);
    }
    else {
*/    
        zero_pot(FMMV);
        t0 = (double) clock()/CLOCKS_PER_SEC;
        EVAL_M(FMMV, box_B, box_A);
        t1 = (double) clock()/CLOCKS_PER_SEC - t0;
        printf("time (EVAL_M)       : %7.4f\n", t1);
        FMMV->potentials = potM;
        FMMV->gradients = gradM;
        backcopy_pot(FMMV);

        t0 = (double) clock()/CLOCKS_PER_SEC;
        GEN_L(FMMV, box_B, box_A);
        t1 = (double) clock()/CLOCKS_PER_SEC - t0;
        printf("time (GEN_L)        : %7.4f\n", t1);
/*    } */
    zero_pot(FMMV);
    if (with_L2L) {
        int c;
        init_L2L(FMMV, -1);
        init_L2L(FMMV, level_B);
        t0 = (double) clock()/CLOCKS_PER_SEC;
        L2L(FMMV, box_B);
        t1 = (double) clock()/CLOCKS_PER_SEC - t0;
        printf("time (L2L)          : %7.4f\n", t1);
        t0 = (double) clock()/CLOCKS_PER_SEC;
        for (c=0; c<8; c++) {
            EVAL_L(FMMV, box_B->child[c]);    
        }
        t1 = (double) clock()/CLOCKS_PER_SEC - t0;
        printf("time (EVAL_L)       : %7.4f\n", t1);
    }
    else {
        t0 = (double) clock()/CLOCKS_PER_SEC;
        EVAL_L(FMMV, box_B);
        t1 = (double) clock()/CLOCKS_PER_SEC - t0;
        printf("time (EVAL_L)       : %7.4f\n", t1);
    }
    FMMV->potentials = potL;
    FMMV->gradients = gradL;
    backcopy_pot(FMMV);

    ida_free(FMMV);

    if (with_gradients) {
	relL2Err(0, 0, 0); /* init */
	relL2Err3(0, NULL, NULL); /* init */

	if (printResult) {
  	    printf("\n          Pot(M)   Pot(exact)     rel.err.  rel.err.(grad)\n");
	    printf("============================================================\n");
	}	
	for (i=0; i<NTargets; i++) {
	    relL2Err(1, potM[i], exPot[i]); /* accumulate */
	    relL2Err3(1, gradM[i], exGrad[i]); /* accumulate */
	    if (printResult) {
	        printf("%3i %12.4e %12.4e %12.4e %12.4e\n", i, (double) potM[i], (double) exPot[i], (double) relErr(potM[i], exPot[i]),(_FLOAT_) relErr3(gradM[i], exGrad[i]));
	    }	
  	}
	errL2 = relL2Err(2, 0, 0); /*finalize*/ 
	errL2grad = relL2Err3(2, NULL, NULL); /*finalize*/ 
	printf ("\nerr_L2(potM) = %.4e  err_L2(gradM) = %.4e\n", errL2, errL2grad);
        
	relL2Err(0, 0, 0); /* init */
	relL2Err3(0, NULL, NULL); /* init */

	if (printResult) {
  	    printf("\n          Pot(L)   Pot(exact)     rel.err.  rel.err.(grad)\n");
	    printf("============================================================\n");
	}	
	for (i=0; i<NTargets; i++) {
	    relL2Err(1, potL[i], exPot[i]); /* accumulate */
	    relL2Err3(1, gradL[i], exGrad[i]); /* accumulate */
	    if (printResult) {
	        // printf("%3i %12.4e %12.4e %12.4e %12.4e\n", i, (double) potL[i], (double) exPot[i], (double) relErr(potL[i], exPot[i]),(_FLOAT_) relErr3(gradL[i], exGrad[i]));
	        printf("%3i %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e\n",
                       i, (double) potL[i], (double) exPot[i], (double) relErr(potL[i], exPot[i]),
                        (double) gradL[i][0]/  exGrad[i][0],
                        (double) gradL[i][1]/  exGrad[i][1],
                        (double) gradL[i][2]/  exGrad[i][2],
                       (_FLOAT_) relErr3(gradL[i], exGrad[i]));
	    }	
  	}
	errL2 = relL2Err(2, 0, 0); /*finalize*/ 
	errL2grad = relL2Err3(2, NULL, NULL); /*finalize*/ 
	printf ("\nerr_L2(potL) = %.4e  err_L2(gradL) = %.4e\n", errL2, errL2grad);
    }
    else {
        relL2Err(0, 0, 0);
        if (printResult) {
            printf("\n          Pot(M)   Pot(exact)     rel.err.\n");
            printf("============================================\n");
        }
        for (i=0; i<NTargets; i++) {
            relL2Err(1, potM[i], exPot[i]); /* accumulate */
            if (printResult) {
	       printf("%3i %12.4e %12.4e %12.4e\n", i, (double) potM[i], (double) exPot[i], (double) relErr(potM[i], exPot[i]));
            }   
        }
        errL2 = relL2Err(2, 0, 0); /*finalize*/ 
        printf ("\nerr_L2(potM) = %.4e\n", errL2);
    
        relL2Err(0, 0, 0);
        if (printResult) {
            printf("\n          Pot(L)   Pot(exact)     rel.err.\n");
            printf("============================================\n");
        }
        for (i=0; i<NTargets; i++) {
            relL2Err(1, potL[i], exPot[i]); /* accumulate */
            if (printResult) {
	       printf("%3i %12.4e %12.4e %12.4e\n", i, (double) potL[i], (double) exPot[i], (double) relErr(potL[i], exPot[i]));
            }   
        }
        errL2 = relL2Err(2, 0, 0); /*finalize*/ 
        printf ("\nerr_L2(potL) = %.4e\n", errL2);
    }
    
    // TODO: deallocate...
    finish_all(FMMV);

    return 0;
}
