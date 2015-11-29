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

#ifndef _FMMV_COMMON_H_
#define _FMMV_COMMON_H_

/* The following should be (un-)defined in Makefile.inc :
   #define USE_PTHREADS
*/   
/* The following should be (un-)defined in statistics.c :
   #define USE_PAPI
*/

#if (FMM_DIM == 2)
	#define FMM_CHILDS_PER_BOX 4
	enum Childs { SW, NW, SE, NE};

	#define FMM_NEIGHBORS_PER_BOX 8 
	enum Neighbors { SW_, MW_, NW_, SM_,  NM_, SE_, ME_, NE_};

	#define FMM_X_DIRECTIONS 4
	enum Boundarys { B_S=1, B_N=2, B_W=4, B_E=8 };
	enum Xdirections { XE, XN, XW, XS };

	/* Note: sizes = # of real numbers = 2 * (# of complex numbers ) */
	#define FMM_SIZE_OF_M(P) (2*((P)+1))
	#define FMM_SIZE_OF_L(P) (2*((P)+1))

#elif (FMM_DIM == 3)
	#define FMM_CHILDS_PER_BOX 8
	enum Childs { SWD, NWD, SED, NED, SWU, NWU, SEU, NEU};

	#define FMM_NEIGHBORS_PER_BOX 26
	enum Neighbors {  SWD_, MWD_, NWD_, SMD_, MMD_, NMD_, SED_, MED_, NED_,
		          SWM_, MWM_, NWM_, SMM_,       NMM_, SEM_, MEM_, NEM_,
	        	  SWU_, MWU_, NWU_, SMU_, MMU_, NMU_, SEU_, MEU_, NEU_};

	#define FMM_X_DIRECTIONS 6
	enum Boundarys { B_S=1, B_N=2, B_W=4, B_E=8, B_D=16, B_U=32};
	enum Xdirections { XE, XN, XW, XS, XU, XD };

	/* Note: sizes = # of real numbers = 2 * (# of complex numbers ) */
	#define FMM_SIZE_OF_M(P) (((P)+2)*((P)+1))
	#define FMM_SIZE_OF_L(P) (((P)+2)*((P)+1))
#endif

typedef struct box_ *BoxPtr;

#ifdef USE_PTHREADS
	#define _REENTRANT
	#include<pthread.h>
#endif 

typedef struct box_ {
    BoxPtr nextSourceBox;    /* next Box in list */
    BoxPtr nextTargetBox;
    BoxPtr parent;
    BoxPtr child[FMM_CHILDS_PER_BOX];
    BoxPtr neighbor[FMM_NEIGHBORS_PER_BOX];

#if (FMM_DIM == 1) /* coordinates of center */
    _FLOAT_ x; 
#elif (FMM_DIM == 2)
    _FLOAT_ x, y; 
#elif (FMM_DIM == 3)
    _FLOAT_ x, y, z; 
#endif
    
    short int level;      /* the level to which this box belongs */
    unsigned char whichChild;

    unsigned char atBoundary;

    unsigned int firstParticle;  
    unsigned int noOfParticles;
    
    unsigned int firstTarget;  
    unsigned int noOfTargets;
    
#ifdef USE_PTHREADS
    pthread_mutex_t mutex;
#endif    

    _FLOAT_ *M;  /* multipole expansion */
    _FLOAT_ *L;  /* local expansion */
    
    int noOfXin;  
    _FLOAT_ *X[FMM_X_DIRECTIONS]; /* exponential expansion */

    int noOfXin2; 
    _FLOAT_ *X2[FMM_X_DIRECTIONS];  /* exponential expansion for ws=2 */

    FMM_ADDITIONAL_DATA_IN_BOX

} Box;


#define isSource(box) ((box) && (box)->noOfParticles)
#define isTarget(box) ((box) && (box)->noOfTargets)

#if (FMM_DIM == 2)

#define isChildless(box) \
(!((box)->child[0] || (box)->child[1] || (box)->child[2] || (box)->child[3])) 

#define hasSourceChilds(box) \
  (isSource((box)->child[0]) || isSource((box)->child[1]) \
|| isSource((box)->child[2]) || isSource((box)->child[3]))

#define hasTargetChilds(box) \
  (isTarget((box)->child[0]) || isTarget((box)->child[1]) \
|| isTarget((box)->child[2]) || isTarget((box)->child[3])) 

#elif (FMM_DIM == 3)

#define isChildless(box) \
(!((box)->child[0] || (box)->child[1] || (box)->child[2] || (box)->child[3] || \
   (box)->child[4] || (box)->child[5] || (box)->child[6] || (box)->child[7]))  

#define hasSourceChilds(box) \
  (isSource((box)->child[0]) || isSource((box)->child[1]) \
|| isSource((box)->child[2]) || isSource((box)->child[3]) \
|| isSource((box)->child[4]) || isSource((box)->child[5]) \
|| isSource((box)->child[6]) || isSource((box)->child[7]))

#define hasTargetChilds(box) \
  (isTarget((box)->child[0]) || isTarget((box)->child[1]) \
|| isTarget((box)->child[2]) || isTarget((box)->child[3]) \
|| isTarget((box)->child[4]) || isTarget((box)->child[5]) \
|| isTarget((box)->child[6]) || isTarget((box)->child[7]))

#endif

#define isChildlessSource(box) \
(isSource(box) && !hasSourceChilds(box))

#define isChildlessTarget(box) \
(isTarget(box) && !hasTargetChilds(box))

typedef struct fmmvHandle_ *FmmvHandlePtr;

typedef struct  fmmvHandle_ {
	int magicNumber;
	int NParticles;
	int NTargets;
	int typeSources;
	int computeGradient;
	int splitThreshold;
	int splitTargetThreshold;
	int maxLevel;
	int maxTargetLevel;
	int directEvalThreshold;
	int periodicBoundaryConditions;
	int extrinsicCorrection;
	int useHilbertOrder;
	int directEvalAccuracy;
	int useFarfieldNearfieldThreads;
        int reducedScheme;
        _FLOAT_ screeningParameter;  /* parameter in exponent of Yukawa potential */
	_FLOAT_ scale;

        int maxNoOfStoredXin; 
        int noOfStoredXin; 
	unsigned int maxAllocatedMemory;	
	unsigned int allocatedMemory;
	unsigned long long int noOfDirectInteractions;
	struct FmmvStatistics statistics;
	
	int pM;      /* length of multipole expansions */
	int pL;      /* length of local expansions */
        int s_eps;   /* parameter for length of exponential expansions */
        int s_exp;   /* length of exponential expansions */
	
	int ws; /* 1 or 2 */

	Box *firstSourceBoxOfLevel[52];
	Box *firstTargetBoxOfLevel[52];

	_FLOAT_ (*particles)[FMM_DIM];
	_FLOAT_ (*targets)[FMM_DIM];
	_FLOAT_ (*particles0)[FMM_DIM];
	_FLOAT_ (*targets0)[FMM_DIM];

	int *perm;
	int *permTargets;

	int dataKind; 
	
        void (*gen_M)(FmmvHandlePtr, Box *box);
        void (*eval_L)(FmmvHandlePtr, Box *box);
        void (*gen_L)(FmmvHandlePtr, Box *target, Box *source);
        void (*eval_M)(FmmvHandlePtr, Box *target, Box *source);
        void (*gen_L_eval_M)(FmmvHandlePtr, Box *list3, Box *list4);
        void (*eval_direct)(FmmvHandlePtr, Box *target, Box *source);
        void (*eval_direct_periodic)(FmmvHandlePtr, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy
		#if (FMM_DIM>=3)
		, _FLOAT_ dz
		#endif
		);
        void (*extrinsic_correction)(FmmvHandlePtr);
	
	_FLOAT_ *D_X2X; /* for storing the diagonal X2X translation operators */

	FMM_ADDITIONAL_GLOBAL_DATA  /* from defs.h */

	void *DATA;    
} FmmvHandle;


char* fmmv_initialize(FmmvHandle **FMMV, struct FmmvOptions *options, struct FmmvStatistics *statistics);
char* fmmv_evaluate(FmmvHandle *FMMV, struct FmmvStatistics *statistics);
char* fmmv_finalize(FmmvHandle *FMMV, struct FmmvStatistics *statistics);

int init_all(FmmvHandle *fh);
int finish_all(FmmvHandle *fh);
void double2single(int n, double *source, float *target);

void genList1(Box *box, Box **list1, int *n);
void genList3(Box *box, Box **list3, int *n);
void genList4(Box *box, Box **list4, int *n);
void genList1_ST(Box *box, Box **list1, int *n);
void genList3_ST(Box *box, Box **list3, int *n);
void genList4_ST(Box *box, Box **list4, int *n);

void gen_neighbor2_list(Box *box, Box **list);
void genNeighborRelations(Box *box);
void genPeriodicNeighborRelations(Box *root);

int buildTree(FmmvHandle *fh);
int buildTree_ST(FmmvHandle *fh);
void genTreeStatistics(FmmvHandle *fh);
void genTreeStatistics_ST(FmmvHandle *fh);
void freeTree(FmmvHandle *fh, Box *box);
void freeTree_ST(FmmvHandle *fh, Box *box);

enum KindOfBox {SOURCE=1, TARGET=2};

Box* createBox(FmmvHandle *fh, Box *parent, 
	  _FLOAT_ x, 
	  _FLOAT_ y,
	#if (FMM_DIM>=3)
	  _FLOAT_ z,
	#endif	  
	  unsigned int first, 
	  unsigned int noOf,
	  unsigned int level, 
	  unsigned int whichChild,
	  enum KindOfBox kindOFBox); 	  

typedef struct {
	int thread;
	FmmvHandle *fh;
} GenericFmmThreadArg;	

void* non_adaptive_fmm(GenericFmmThreadArg *arg);
void* non_adaptive_fmm_periodic(GenericFmmThreadArg *arg);
void* adaptive_fmm(GenericFmmThreadArg *arg);
void* adaptive_fmm_ST(GenericFmmThreadArg *arg);
void* non_adaptive_fmm_ws2(GenericFmmThreadArg *arg);
void* non_adaptive_fmm_periodic_ws2(GenericFmmThreadArg *arg);
void* adaptive_fmm_ST_ws2(GenericFmmThreadArg *arg);
void* adaptive_fmm_ST_dipole_grad_ws2(GenericFmmThreadArg *arg);

void stat_init(FmmvHandle *fh);
void stat_start(FmmvHandle *fh, enum StatStep step);	
void stat_stop(FmmvHandle *fh, enum StatStep step);

void init_M2M(FmmvHandle *fh, int level);
void init_M2L(FmmvHandle *fh, int level);
void init_L2L(FmmvHandle *fh, int level);
void finish_M2M(FmmvHandle *fh);
void finish_M2L(FmmvHandle *fh);
void finish_L2L(FmmvHandle *fh);

void M2M(FmmvHandle *fh, Box *box);
void L2L(FmmvHandle *fh, Box *box);
void M2L(FmmvHandle *fh, Box *box);
void M2L_ws2(FmmvHandle *fh, Box *box);
void M2L_ws2_reduced(FmmvHandle *fh, Box *box);
void periodic_lattice_M2L(FmmvHandle *fh, _FLOAT_ *M, _FLOAT_ *L);

int ida_allocate(FmmvHandle *fh); 
void ida_free(FmmvHandle *fh); 
void copy_particles(FmmvHandle *fh);
void copy_charges(FmmvHandle *fh);
void zero_pot(FmmvHandle *fh);
void backcopy_pot(FmmvHandle *fh);


#include<stdlib.h> /* malloc */

#define SIMD_ALIGN\

#define FMMV_MALLOC(FMMV, n) malloc(n); \
	FMMV->allocatedMemory+=(n); \
	if (FMMV->allocatedMemory>FMMV->maxAllocatedMemory) \
		FMMV->maxAllocatedMemory=FMMV->allocatedMemory;

#define FMMV_FREE(FMMV, x, n) if (x) FMMV->allocatedMemory-=(n); \
	free(x); \


void VEC_mul(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_scale(const int n, const _FLOAT_ d, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_mul_simd2(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_mul2_simd2(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_scale_simd2(const int n, const _FLOAT_ d, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_add(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_add_simd2(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y);
_FLOAT_* VEC_add2(FmmvHandle *fh, const int n, const _FLOAT_ *x, _FLOAT_ *y);
_FLOAT_* VEC_add2_simd2(FmmvHandle *fh, const int n, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_sub(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_sub_simd2(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_mul_c(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_mul_cj(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_mul_ccj(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y, _FLOAT_ *z);
void VEC_mul_c_rrii(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_mul_cj_rrii(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_mul_ccj_rrii(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y, _FLOAT_ *z);
void VEC_mul_c_rrii_simd2(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_mul_cj_rrii_simd2(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_mul_ccj_rrii_simd2(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y, _FLOAT_ *z);

void VEC_mul_simd4(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_mul4_simd4(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_scale_simd4(const int n, const _FLOAT_ d, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_add_simd4(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y);
_FLOAT_* VEC_add2_simd4(FmmvHandle *fh, const int n, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_add_simd4(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_sub_simd4(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_mul_c_rrii_simd4(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_mul_cj_rrii_simd4(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y);
void VEC_mul_ccj_rrii_simd4(const int n, const _FLOAT_ *d, const _FLOAT_ *x, _FLOAT_ *y, _FLOAT_ *z);

void neg(const int n, const int *vec, _FLOAT_ *x);
void neg_simd2(const int n, const int *vec, _FLOAT_ *x);
void neg_simd4(const int n, const int *vec, _FLOAT_ *x);

void VEC_conj(const int n, _FLOAT_ *y);
void VEC_conj2_simd2(const int n, _FLOAT_ *y);
void VEC_conj4_simd4(const int n, _FLOAT_ *y);

#if (FMM_SIMD_VECTOR_LENGTH==4)
  #define VEC_MUL(n, d, x, y) VEC_mul_simd4(n, d, x, y)
  #define VEC_SCALE(n, d, x, y) VEC_scale_simd4(n, d, x, y)
  #define VEC_ADD(n, d, x, y) VEC_add_simd4(n, d, x, y)
  #define VEC_ADD2(FMMV, n, x, y) VEC_add2_simd4(FMMV, n, x, y)
  #define VEC_SUB(n, d, x, y) VEC_sub_simd4(n, d, x, y)
#elif (FMM_SIMD_VECTOR_LENGTH==2)
  #define VEC_MUL(n, d, x, y) VEC_mul_simd2(n, d, x, y)
  #define VEC_SCALE(n, d, x, y) VEC_scale_simd2(n, d, x, y)
  #define VEC_ADD(n, d, x, y) VEC_add_simd2(n, d, x, y)
  #define VEC_ADD2(FMMV, n, x, y) VEC_add2_simd2(FMMV, n, x, y)
  #define VEC_SUB(n, d, x, y) VEC_sub_simd2(n, d, x, y)
#else
  #define VEC_MUL(n, d, x, y) VEC_mul(n, d, x, y)
  #define VEC_SCALE(n, d, x, y) VEC_scale(n, d, x, y)
  #define VEC_ADD(n, d, x, y) VEC_add(n, d, x, y)
  #define VEC_ADD2(FMMV, n, x, y) VEC_add2(FMMV, n, x, y)
  #define VEC_SUB(n, d, x, y) VEC_sub(n, d, x, y)
#endif

#ifdef X_RRII
#if (FMM_SIMD_VECTOR_LENGTH==4)
  #define VEC_MUL_C(n, d, x, y) VEC_mul_c_rrii_simd4(n, d, x, y)
  #define VEC_MUL_CJ(n, d, x, y) VEC_mul_cj_rrii_simd4(n, d, x, y)
  #define VEC_MUL_CCJ(n, d, x, y, z) VEC_mul_ccj_rrii_simd4(n, d, x, y, z)
#elif (FMM_SIMD_VECTOR_LENGTH==2)
  #define VEC_MUL_C(n, d, x, y) VEC_mul_c_rrii_simd2(n, d, x, y)
  #define VEC_MUL_CJ(n, d, x, y) VEC_mul_cj_rrii_simd2(n, d, x, y)
  #define VEC_MUL_CCJ(n, d, x, y, z) VEC_mul_ccj_rrii_simd2(n, d, x, y, z)
#else
  #define VEC_MUL_C(n, d, x, y) VEC_mul_c_rrii(n, d, x, y)
  #define VEC_MUL_CJ(n, d, x, y) VEC_mul_cj_rrii(n, d, x, y)
  #define VEC_MUL_CCJ(n, d, x, y, z) VEC_mul_ccj_rrii(n, d, x, y, z)
#endif
#else
  #define VEC_MUL_C(n, d, x, y) VEC_mul_c(n, d, x, y)
  #define VEC_MUL_CJ(n, d, x, y) VEC_mul_cj(n, d, x, y)
  #define VEC_MUL_CCJ(n, d, x, y, z) VEC_mul_ccj(n, d, x, y, z)
#endif

#include <string.h>
#define VEC_COPY(n, x , y) memcpy(y, x, (n)*sizeof(_FLOAT_))


#endif /* _FMMV_COMMON_H_ */
