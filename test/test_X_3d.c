#ifdef USE_SINGLE_PRECISION 
#define FMM_PRECISION 0
#else
#define FMM_PRECISION 1
#endif

#define X_RRII
#include "_fmmv.h"

#include<stdio.h>
#include<math.h>
#include<time.h>
#include<fenv.h> 
#include"test_utilities.h"

void M2X(FmmvHandle *FMMV, int dir, Box *box, _FLOAT_ *X1, _FLOAT_ *X2);
void X2L(FmmvHandle *FMMV, int dir, Box *box, int reduced);
void init_M2L(FmmvHandle *FMMV, int level);
void finish_M2L(FmmvHandle *FMMV);

#ifdef DOUBLE_FOR_INIT
  #define _FLOAT1_ double
#else
  #define _FLOAT1_ _FLOAT_
#endif

void gen_diag_X2X(_FLOAT1_ beta, _FLOAT1_ *lambda, int *M, int s_eps, int s_exp, _FLOAT1_ x, _FLOAT1_ y, _FLOAT1_ z, _FLOAT1_ *d);



int main(int argc, char**argv)
{
    int err;
    int i;
    FmmvHandle _FMMV;
    FmmvHandle *FMMV = &_FMMV;
    void (*GEN_M)(FmmvHandle *FMMV, Box *box);
    void (*GEN_L)(FmmvHandle *FMMV, Box *target, Box *source);
    void (*EVAL_L)(FmmvHandle *FMMV, Box *box);
    void (*EVAL_DIRECT)(FmmvHandle *FMMV, Box *target, Box *source) = FMMV->eval_direct;

    _FLOAT_ X1[8*FMM_S_EXP_MAX];
    _FLOAT_ X2[8*FMM_S_EXP_MAX];
    _FLOAT_ X[FMM_S_EXP_MAX];
    _FLOAT_ D[FMM_S_EXP_MAX];


    Box _box_A;
    Box _box_B;
    Box _parent_A;
    Box _parent_B;
    int which_child_A = SWD;
    int which_child_B = SWD;


    Box *box_A = &_box_A;
    Box *box_B = &_box_B;
    Box *parent_A = &_parent_A;
    Box *parent_B = &_parent_B;

    int NParticles = 100;
    int NTargets = 100;
     _FLOAT_ (*particles)[3];
     _FLOAT_ (*targets)[3] = NULL;
     _FLOAT_ *charges;
     _FLOAT_ x, y, z;
     _FLOAT_ x_B, y_B, z_B;
     _FLOAT_ x_A, y_A, z_A;
     _FLOAT_ *potX;
     _FLOAT_ *potL;
     _FLOAT_ *exPot;
     _FLOAT_ errL2;
     double beta = 0.0;
     int pM = 8;
     int pL = 8;
     int s = 8;
     int *perm;
     int *permTargets;
     int printResult=0;

    int level = 2;
    _FLOAT_ size;
    
    _FLOAT_ dx = 1;
    _FLOAT_ dy = 1;
    _FLOAT_ dz = 2;

  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );


    FMMV->dataKind = FMM_ST_STANDARD;
    GEN_M = gen_M_ST_standard;
    GEN_L = gen_L_ST_standard;
    EVAL_L = eval_L_ST_standard;
#if (FMM_PRECISION==0)
    EVAL_DIRECT = eval_direct_ST_standard_acc1; 
#else
    EVAL_DIRECT = eval_direct_ST_standard_acc2;
#endif

    set_from_command_line_double(argc, argv, "beta", &beta);
    set_from_command_line_bool(argc, argv, "printResult", &printResult);
    set_from_command_line_int(argc, argv, "NParticles", &NParticles);
    set_from_command_line_int(argc, argv, "NTargets", &NTargets);
    set_from_command_line_int(argc, argv, "pM", &pM);
    set_from_command_line_int(argc, argv, "pL", &pL);
    set_from_command_line_int(argc, argv, "s", &s);
    set_from_command_line_int(argc, argv, "level", &level);
    size = ldexp(1.0, -level);
    
    particles = (_FLOAT_ (*)[3]) calloc(NParticles, 3*sizeof(_FLOAT_));  
    charges = (_FLOAT_ *) calloc(NParticles, sizeof(_FLOAT_));  
    perm = (int *) calloc(NParticles, sizeof(int));  
    permTargets = (int *) calloc(NTargets, sizeof(int));  
    targets = (_FLOAT_ (*)[3]) calloc(NTargets, 3*sizeof(_FLOAT_));  
    potX = (_FLOAT_ *) calloc(NTargets, sizeof(_FLOAT_));
    potL = (_FLOAT_ *) calloc(NTargets, sizeof(_FLOAT_));
    exPot = (_FLOAT_ *) calloc(NTargets, sizeof(_FLOAT_));

    FMMV->beta = beta;
    FMMV->pM = pM;
    FMMV->pL = pL;
    FMMV->s_eps = s;
    FMMV->particles = particles;
    FMMV->charges = charges;
    FMMV->perm = perm;
    FMMV->NParticles = NParticles;
    FMMV->targets = targets;
    FMMV->permTargets = permTargets;
    FMMV->NTargets = NTargets;
    FMMV->scale = 1.0;
    FMMV->lambda = 0;

    my_srand(1481765933);
    x_A = 0.5*size;
    y_A = 0.5*size;
    z_A = 0.5*size;

    box_A->firstParticle = 0;
    box_A->noOfParticles = NParticles;
    parent_A->firstParticle = 0;
    parent_A->noOfParticles = NParticles;
    box_A->x = x_A;
    box_A->y = y_A;
    box_A->z = z_A;
    box_A->level = level;
    parent_A->level = level-1;
    box_A->M = 0; 
    
    x_B = x_A + size*dx;
    y_B = y_A + size*dy;
    z_B = z_A + size*dz;

    box_B->firstTarget = 0;
    box_B->noOfTargets = NTargets;
    parent_B->firstTarget = 0;
    parent_B->noOfTargets = NTargets;
    box_B->x = x_B;
    box_B->y = y_B;
    box_B->z = z_B;
    box_B->level = level;
    parent_B->level = level-1;
    box_B->L = 0; 
    
    for (i=0; i<8; i++) {
        parent_A->child[i] = 0;
        parent_B->child[i] = 0;
    }
    
    parent_A->child[which_child_A] = box_A;
    parent_B->child[which_child_B] = box_B;
    
    for (i=0;i<NParticles;i++) {
        x = x_A + size*(my_rand() / (_FLOAT_) MY_RAND_MAX - 0.5);
        y = y_A + size*(my_rand() / (_FLOAT_) MY_RAND_MAX - 0.5);
        z = z_A + size*(my_rand() / (_FLOAT_) MY_RAND_MAX - 0.5);
        particles[i][0] = x;
        particles[i][1] = y;
        particles[i][2] = z;
    }    
    for (i=0;i<NParticles;i++) {
        perm[i] = i;
        charges[i] = 1.0;
    }

    for (i=0;i<NTargets;i++) {
        x = x_B + size*(my_rand() / (_FLOAT_) MY_RAND_MAX - 0.5);
        y = y_B + size*(my_rand() / (_FLOAT_) MY_RAND_MAX - 0.5);
        z = z_B + size*(my_rand() / (_FLOAT_) MY_RAND_MAX - 0.5);
        targets[i][0] = x;
        targets[i][1] = y;
        targets[i][2] = z;
    }    

    for (i=0;i<NTargets;i++) {
        permTargets[i] = i;
    }
    
    init_all(FMMV);

    err = ida_allocate(FMMV);
    copy_particles(FMMV);
    copy_charges(FMMV);
    
    zero_pot(FMMV);
    EVAL_DIRECT(FMMV, box_B, box_A);
    FMMV->potentials = exPot;
    backcopy_pot(FMMV);

    GEN_M(FMMV, box_A);

    GEN_L(FMMV, box_B, box_A);

    zero_pot(FMMV);
    EVAL_L(FMMV, box_B);
    FMMV->potentials = potL;
    backcopy_pot(FMMV);

    /**********************************************/
    for (i=0;i<6;i++) {
       box_B->X[i] = 0;
    }
    box_B->X[XU] = X;
    for (i=0;i<pL*(pL+1);i++) {
       box_B->L[i] = 0.0;
    }   

    init_M2L(FMMV, -1);
    init_M2L(FMMV, level);

    M2X(FMMV, 0, parent_A, X1, X2);
    gen_diag_X2X(ldexp(beta,-level), FMMV->lambda, FMMV->M, FMMV->s_eps, FMMV->s_exp, dx, dy, dz, D); 
    VEC_MUL_C(FMMV->s_exp/2, D, X1+FMMV->s_exp*which_child_A, box_B->X[XU]);
    X2L(FMMV, 0, parent_B, 0);

    finish_M2L(FMMV);

    zero_pot(FMMV);
    EVAL_L(FMMV, box_B);
    FMMV->potentials = potX;
    backcopy_pot(FMMV);
    /**********************************************/


    ida_free(FMMV);

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

    relL2Err(0, 0, 0);
    if (printResult) {
        printf("\n          Pot(X)   Pot(exact)     rel.err.\n");
        printf("============================================\n");
    }
    for (i=0; i<NTargets; i++) {
        relL2Err(1, potX[i], exPot[i]); /* accumulate */
        if (printResult) {
            printf("%3i %12.4e %12.4e %12.4e\n", i, (double) potX[i], (double) exPot[i], (double) relErr(potX[i], exPot[i]));
        }   
    }
    errL2 = relL2Err(2, 0, 0); /*finalize*/ 
    printf ("\nerr_L2(potX) = %.4e\n", errL2);
  
    // TODO: deallocate...
    finish_all(FMMV);

    return 0;
}
