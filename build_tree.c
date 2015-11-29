/*
 * FMMV - the Fastest Multipole Method of Vienna
 * Copyright (c) 2006 Harald Hofstaetter
 * Institute for Analysis and Scientific Computing
 * Vienna University of Technology
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


#include "_fmmv.h"

#include<assert.h>
#include<stdlib.h>
#include<string.h>  /* memset */
#include<math.h>    /* ldexp */  
#ifdef USE_PTHREADS
   #include<pthread.h>
#endif


/* with gcc for the getbit functions the option -fno-strict-aliases
 * should be used,
 */

/* TODO: comment why  9 resp. 12 instead of 8 resp. 11 ... */

#if (FMM_PRECISION == 0)
// #define GETBIT(x, k) getbitf(x, k)
   #define GETBIT(x, k) ((*((unsigned int*) &(x)) << (9 + (k))) & 0x80000000L )
#else
// #define GETBIT(x, k) getbit(x, k)
   #define GETBIT(x, k) ((*((unsigned long long int*) &(x)) << (12 + (k))) & 0x8000000000000000LL )
#endif

static int 
radixSplit(_FLOAT_ particles[][FMM_DIM], int perm[], int index, int left, int right, int level, int reverse)
{
     int i=left;
     int j=right-1;
     
     _FLOAT_ tx, ty;
     #if (FMM_DIM>=3) 
     _FLOAT_ tz;
     #endif
     int tp;
     
     if (reverse) { /* GETBIT -> !GETBIT, !GETBIT -> GETBIT */ 
         while (1){
            while (i<j && (GETBIT(particles[i][index],level))) i++;
            if (i>=j) break;
            while (j>i && (!GETBIT(particles[j][index],level))) j--;
	    if (i>=j) break;

            tx = particles[i][0]; 
            ty = particles[i][1];
	    #if (FMM_DIM>=3) 
            tz = particles[i][2]; 
	    #endif
	    
	    particles[i][0] = particles[j][0]; 
	    particles[i][1] = particles[j][1]; 
	    #if (FMM_DIM>=3) 
	    particles[i][2] = particles[j][2]; 
	    #endif

	    particles[j][0] = tx;
	    particles[j][1] = ty;
	    #if (FMM_DIM>=3) 
	    particles[j][2] = tz;
	    #endif

            tp = perm[i];
	    perm[i] = perm[j]; 
	    perm[j] = tp;
         }
	 for (j=left; (j<right) && GETBIT(particles[j][index],level) ; j++) ;
     } 
     else {
         while (1){
            while (i<j && (!GETBIT(particles[i][index],level))) i++;
            if (i>=j) break;
            while (j>i && (GETBIT(particles[j][index],level))) j--;
	    if (i>=j) break;
	
            tx = particles[i][0]; 
            ty = particles[i][1]; 
	    #if (FMM_DIM>=3) 
            tz = particles[i][2]; 
	    #endif
	    
	    particles[i][0] = particles[j][0]; 
	    particles[i][1] = particles[j][1]; 
	    #if (FMM_DIM>=3) 
	    particles[i][2] = particles[j][2]; 
	    #endif

	    particles[j][0] = tx;
	    particles[j][1] = ty;
	    #if (FMM_DIM>=3) 
	    particles[j][2] = tz;
	    #endif

            tp = perm[i];
	    perm[i] = perm[j]; 
	    perm[j] = tp;
         }
	 for (j=left; (j<right) && !GETBIT(particles[j][index],level) ; j++) ;
     }
     return j;
}	

/********************************************************************************/
#if (FMM_DIM==2)
static int
createChilds(FmmvHandle *FMMV, Box *parent, int level, enum KindOfBox kindOfBox)
{
   Box * child; 
   int err = 0;
   int n0, n2, n4, n6, n8;
   _FLOAT_ x,y;
   _FLOAT_ h;

   if (level==FMMV->maxLevel) return 0;

   x = parent->x;
   y = parent->y;
   h = ldexp(1.0, -(level+2)) / FMMV->scale;

   if (kindOfBox & SOURCE) {
       if (parent->noOfParticles < FMMV->splitThreshold) {
           return 0;
       }
       n0 = parent->firstParticle;
       n8 = n0 + parent->noOfParticles;

       n4 = radixSplit(FMMV->particles0, FMMV->perm, 0, n0, n8, level, 0);
       n2 = radixSplit(FMMV->particles0, FMMV->perm, 1, n0, n4, level, 0);
       n6 = radixSplit(FMMV->particles0, FMMV->perm, 1, n4, n8, level, 0);
   }
   else { /* kindOfBox == TARGET */
       if (parent->noOfTargets < FMMV->splitTargetThreshold) {
           return 0;
       }
       n0 = parent->firstTarget;
       n8 = n0 + parent->noOfTargets;

       n4 = radixSplit(FMMV->targets0, FMMV->permTargets, 0, n0, n8, level, 0);
       n2 = radixSplit(FMMV->targets0, FMMV->permTargets, 1, n0, n4, level, 0);
       n6 = radixSplit(FMMV->targets0, FMMV->permTargets, 1, n4, n8, level, 0);
   }	   

   if (n2-n0>0) {
       child = createBox(FMMV, parent, x-h, y-h, n0, n2-n0, level+1, SW, kindOfBox);
       if (!child) return -1;
       parent->child[SW] = child;
       err = createChilds(FMMV, child, level+1, kindOfBox);
       if (err) return -1;			  
   }

   if (n4-n2>0) {
       child = createBox(FMMV, parent, x-h, y+h, n2, n4-n2, level+1, NW, kindOfBox);
       if (!child) return -1;
       parent->child[NW] = child;
       err = createChilds(FMMV, child, level+1, kindOfBox);
       if (err) return -1;			  
   }   
   
   if (n6-n4>0) {
       child = createBox(FMMV, parent, x+h, y-h, n4, n6-n4, level+1, SE, kindOfBox);
       if (!child) return -1;
       parent->child[SE] = child;
       err = createChilds(FMMV, child, level+1, kindOfBox);
       if (err) return -1;			  
   }

   if (n8-n6>0) {
       child = createBox(FMMV, parent, x+h, y+h, n6, n8-n6, level+1, NE, kindOfBox);
       if (!child) return -1;
       parent->child[NE] = child;
       err = createChilds(FMMV, child, level+1, kindOfBox);
       if (err) return -1;			  
   }   
   
   return 0;
   
}

/* the following python code genreated the tables below... */
//TODO: Check this for 2D !!!!!
#if 0
/*
m = ['SW', 'NW', 'SE', 'NE']
		
s1 = [ 0, 2, 3, 1 ]
s2 = [ 0, 1, 3, 2 ]

block = {}
for s in [s1, s2]:
	a = s[:]
	for i in range(2):
		b = a[2:]+a[:2]
		block[(a[0],a[3])] = a
		block[(b[0],b[3])] = b
		a = a[1:2]+[a[0], a[3]]+a[2:3]
print 'enum Edges {',
n = 0
for b in block.keys():
	if n==7:
		print 'E_%s_%s' % (m[b[0]],m[b[1]]), 
	else:
		print 'E_%s_%s,' % (m[b[0]],m[b[1]]), 
	n += 1
print "};"	

def mm(a,b):
	M = {2:1, 1:2}
	x= a ^ b
	if x & a: 
		return -M[x]
	else:	
		return M[x]

n = 0
for b in block.keys():
	print '{ / * %i E_%s_%s: * /\n\t{' % (n, m[b[0]],m[b[1]]), 
	for i in range(4):
		if i==3:
			print '%s' %m[block[b][i]],
		else:
			print '%s,' %m[block[b][i]],
	print '},\n\t',
	print '{ E_%s_%s,' % (m[block[b][0]], m[block[b][1]]), 
	print 'E_%s_%s,' % (m[block[b][2]], m[block[b][3]]),
	print 'E_%s_%s},\n\t' % (m[block[b][2]], m[block[b][3]]),
	print '%+i, %+i },' % ( 
		mm(block[b][1], block[b][2]), mm(block[b][0], block[b][1]))
	n += 1


*/
#endif

enum Edges { E_SW_NW, E_NE_SE, 
	E_NW_NE, E_NE_NW, 
	E_SE_SW, E_SE_NE, 
	E_NW_SW, E_SW_SE };

typedef struct {
	int child_order[4];
	int sub_path[3]; 
	int b0, b1;
} HilbertPath;	

static HilbertPath hilbert_path[] =
{
	{ /* 0 E_SW_NW: */
        { SW, SE, NE, NW },
        { E_SW_SE, E_NE_NW, E_NE_NW},
        +2, +1 },
		
	{ /* 1 E_NE_SE: */
        { NE, NW, SW, SE },
        { E_NE_NW, E_SW_SE, E_SW_SE},
        -2, -1 },
	
	{ /* 2 E_NW_NE: */
        { NW, SW, SE, NE },
        { E_NW_SW, E_SE_NE, E_SE_NE},
        +1, -2 },
	
	{ /* 3 E_NE_NW: */
        { NE, SE, SW, NW },
        { E_NE_SE, E_SW_NW, E_SW_NW},
        -1, -2 },
	
	{ /* 4 E_SE_SW: */
        { SE, NE, NW, SW },
        { E_SE_NE, E_NW_SW, E_NW_SW},
        -1, +2 },
	
	{ /* 5 E_SE_NE: */
        { SE, SW, NW, NE },
        { E_SE_SW, E_NW_NE, E_NW_NE},
        +2, -1 },
	
	{ /* 6 E_NW_SW: */
        { NW, NE, SE, SW },
        { E_NW_NE, E_SE_SW, E_SE_SW},
        -2, +1 },
	
	{ /* 7 E_SW_SE: */
        { SW, NW, NE, SE },
        { E_SW_NW, E_NE_SE, E_NE_SE},
        +1, +2 },	
};	

typedef struct {
	_FLOAT_ dx, dy;
} ChildDirection;

static ChildDirection child_direction[] = {
   /* SW */   { -1, -1},
   /* NW */   { -1, +1},
   /* SE */   { +1, -1},
   /* NE */   { +1, +1},
};

static int
createChildsHilbert(FmmvHandle *FMMV, Box *parent, int level, enum KindOfBox kindOfBox, int path_no)
{
   Box * child; 
   int err = 0;
   int n0, n2, n4, n6, n8;
   int c;
   _FLOAT_ x,y;
   _FLOAT_ x1, y1;
   _FLOAT_ h;
   HilbertPath *path = hilbert_path + path_no;

   if (level==FMMV->maxLevel) return 0;

   x = parent->x;
   y = parent->y;
   h = ldexp(1.0, -(level+2)) / FMMV->scale;

   if (kindOfBox & SOURCE) {
       if (parent->noOfParticles < FMMV->splitThreshold) {
           return 0;
       }
       n0 = parent->firstParticle;
       n8 = n0 + parent->noOfParticles;

       n4 = radixSplit(FMMV->particles0, FMMV->perm, abs(path->b0)-1, n0, n8, level, path->b0 < 0);
    
       n2 = radixSplit(FMMV->particles0, FMMV->perm, abs(path->b1)-1, n0, n4, level, path->b1 < 0);

       n6 = radixSplit(FMMV->particles0, FMMV->perm, abs(path->b1)-1, n4, n8, level, path->b1 > 0);
   }
   else { /* kindOfBox == TARGET */
       if (parent->noOfTargets < FMMV->splitTargetThreshold) {
           return 0;
       }
       n0 = parent->firstTarget;
       n8 = n0 + parent->noOfTargets;

       n4 = radixSplit(FMMV->targets0, FMMV->permTargets, abs(path->b0)-1, n0, n8, level, path->b0 < 0);

       n2 = radixSplit(FMMV->targets0, FMMV->permTargets, abs(path->b1)-1, n0, n4, level, path->b1 < 0);
    
       n6 = radixSplit(FMMV->targets0, FMMV->permTargets, abs(path->b1)-1, n4, n8, level, path->b1 > 0);
   }	   

   if (n2-n0>0) {
       c = path->child_order[0];
       x1 = x + h*child_direction[c].dx;
       y1 = y + h*child_direction[c].dy;
       child = createBox(FMMV, parent, x1, y1, n0, n2-n0, level+1, c, kindOfBox);
       if (!child) return -1;
       parent->child[c] = child;
       err = createChildsHilbert(FMMV, child, level+1, kindOfBox, path->sub_path[0]);
       if (err) return -1;			  
   }

   if (n4-n2>0) {
       c = path->child_order[1];
       x1 = x + h*child_direction[c].dx;
       y1 = y + h*child_direction[c].dy;
       child = createBox(FMMV, parent, x1, y1, n2, n4-n2, level+1, c, kindOfBox);
       if (!child) return -1;
       parent->child[c] = child;
       err = createChildsHilbert(FMMV, child, level+1, kindOfBox, path->sub_path[1]);
       if (err) return -1;			  
   }   
   
   if (n6-n4>0) {
       c = path->child_order[2];
       x1 = x + h*child_direction[c].dx;
       y1 = y + h*child_direction[c].dy;
       child = createBox(FMMV, parent, x1, y1, n4, n6-n4, level+1, c, kindOfBox);
       if (!child) return -1;
       parent->child[c] = child;
       err = createChildsHilbert(FMMV, child, level+1, kindOfBox, path->sub_path[1]);
       if (err) return -1;			  
   }


   if (n8-n6>0) {
       c = path->child_order[3];
       x1 = x + h*child_direction[c].dx;
       y1 = y + h*child_direction[c].dy;
       child = createBox(FMMV, parent, x1, y1, n6, n8-n6, level+1, c, kindOfBox);
       if (!child) return -1;
       parent->child[c] = child;
       err = createChildsHilbert(FMMV, child, level+1, kindOfBox, path->sub_path[2]);
       if (err) return -1;			  
   }   
   
   return 0;
   
}

#endif /* FMM_DIM==2 */
/********************************************************************************/
#if (FMM_DIM==3)

static int
createChilds(FmmvHandle *FMMV, Box *parent, int level, enum KindOfBox kindOfBox)
{
   Box * child; 
   int err = 0;
   int n0, n1, n2, n3, n4, n5, n6, n7, n8;
   _FLOAT_ x,y,z;
   _FLOAT_ h;

   if (level==FMMV->maxLevel) return 0;

   x = parent->x;
   y = parent->y;
   z = parent->z;
   h = ldexp(1.0, -(level+2)) / FMMV->scale;

   if (kindOfBox & SOURCE) {
       if (parent->noOfParticles < FMMV->splitThreshold) {
           return 0;
       }
       n0 = parent->firstParticle;
       n8 = n0 + parent->noOfParticles;

       n4 = radixSplit(FMMV->particles0, FMMV->perm, 0, n0, n8, level, 0);
    
       n2 = radixSplit(FMMV->particles0, FMMV->perm, 1, n0, n4, level, 0);
       n1 = radixSplit(FMMV->particles0, FMMV->perm, 2, n0, n2, level, 0);
       n3 = radixSplit(FMMV->particles0, FMMV->perm, 2, n2, n4, level, 0);

       n6 = radixSplit(FMMV->particles0, FMMV->perm, 1, n4, n8, level, 0);
       n5 = radixSplit(FMMV->particles0, FMMV->perm, 2, n4, n6, level, 0);
       n7 = radixSplit(FMMV->particles0, FMMV->perm, 2, n6, n8, level, 0);   
   }
   else { /* kindOfBox == TARGET */
       if (parent->noOfTargets < FMMV->splitTargetThreshold) {
           return 0;
       }
       n0 = parent->firstTarget;
       n8 = n0 + parent->noOfTargets;

       n4 = radixSplit(FMMV->targets0, FMMV->permTargets, 0, n0, n8, level, 0);

       n2 = radixSplit(FMMV->targets0, FMMV->permTargets, 1, n0, n4, level, 0);
       n1 = radixSplit(FMMV->targets0, FMMV->permTargets, 2, n0, n2, level, 0);
       n3 = radixSplit(FMMV->targets0, FMMV->permTargets, 2, n2, n4, level, 0);
    
       n6 = radixSplit(FMMV->targets0, FMMV->permTargets, 1, n4, n8, level, 0);
       n5 = radixSplit(FMMV->targets0, FMMV->permTargets, 2, n4, n6, level, 0);
       n7 = radixSplit(FMMV->targets0, FMMV->permTargets, 2, n6, n8, level, 0);   
   }	   

   if (n1-n0>0) {
       child = createBox(FMMV, parent, x-h, y-h, z-h, n0, n1-n0, level+1, SWD, kindOfBox);
       if (!child) return -1;
       parent->child[SWD] = child;
       err = createChilds(FMMV, child, level+1, kindOfBox);
       if (err) return -1;			  
   }

   if (n2-n1>0) {
       child = createBox(FMMV, parent, x-h, y-h, z+h, n1, n2-n1, level+1, SWU, kindOfBox);
       if (!child) return -1;
       parent->child[SWU] = child;
       err = createChilds(FMMV, child, level+1, kindOfBox);
       if (err) return -1;			  
   }   

   if (n3-n2>0) {
       child = createBox(FMMV, parent, x-h, y+h, z-h, n2, n3-n2, level+1, NWD, kindOfBox);
       if (!child) return -1;
       parent->child[NWD] = child;
       err = createChilds(FMMV, child, level+1, kindOfBox);
       if (err) return -1;			  
   }   
   
   if (n4-n3>0) {
       child = createBox(FMMV, parent, x-h, y+h, z+h, n3, n4-n3, level+1, NWU, kindOfBox);
       if (!child) return -1;
       parent->child[NWU] = child;
       err = createChilds(FMMV, child, level+1, kindOfBox);
       if (err) return -1;			  
   }      

   if (n5-n4>0) {
       child = createBox(FMMV, parent, x+h, y-h, z-h, n4, n5-n4, level+1, SED, kindOfBox);
       if (!child) return -1;
       parent->child[SED] = child;
       err = createChilds(FMMV, child, level+1, kindOfBox);
       if (err) return -1;			  
   }

   if (n6-n5>0) {
       child = createBox(FMMV, parent, x+h, y-h, z+h, n5, n6-n5, level+1, SEU, kindOfBox);
       if (!child) return -1;
       parent->child[SEU] = child;
       err = createChilds(FMMV, child, level+1, kindOfBox);
       if (err) return -1;			  
   }   

   if (n7-n6>0) {
       child = createBox(FMMV, parent, x+h, y+h, z-h, n6, n7-n6, level+1, NED, kindOfBox);
       if (!child) return -1;
       parent->child[NED] = child;
       err = createChilds(FMMV, child, level+1, kindOfBox);
       if (err) return -1;			  
   }   
   
   if (n8-n7>0) {
       child = createBox(FMMV, parent, x+h, y+h, z+h, n7, n8-n7, level+1, NEU, kindOfBox);
       if (!child) return -1;
       parent->child[NEU] = child;
       err = createChilds(FMMV, child, level+1, kindOfBox);
       if (err) return -1;			  
   }         
   
   return 0;
   
}

/* the following python code genreated the tables below... */
#if 0
/*

m = ['SWD', 'SWU', 'NWD', 'NWU', 'SED', 'SEU', 'NED', 'NEU']
		
s1 = [ 0, 4, 6, 2,  3, 7, 5, 1]
s2 = [ 0, 1, 5, 4,  6, 7, 3, 2]
s3 = [ 0, 2, 3, 1,  5, 7, 6, 4]

block = {}
for s in [s1, s2, s3]:
	a = s[:]
	for i in range(4):
		b = a[4:]+a[:4]
		block[(a[0],a[7])] = a
		block[(b[0],b[7])] = b
		a = a[1:4]+[a[0], a[7]]+a[4:7]
print 'enum Edges {',
n = 0
for b in block.keys():
	if n==23:
		print 'E_%s_%s' % (m[b[0]],m[b[1]]), 
	else:
		print 'E_%s_%s,' % (m[b[0]],m[b[1]]), 
	n += 1
print "};"	

def mm(a,b):
	M = {4:1, 2:2, 1:3}
	x= a ^ b
	if x & a: 
		return -M[x]
	else:	
		return M[x]

n = 0
for b in block.keys():
	print '{ / * %i E_%s_%s: * /\n\t{' % (n, m[b[0]],m[b[1]]), 
	for i in range(8):
		if i==7:
			print '%s' %m[block[b][i]],
		else:
			print '%s,' %m[block[b][i]],
	print '},\n\t',
	print '{ E_%s_%s,' % (m[block[b][0]], m[block[b][1]]), 
	print 'E_%s_%s,' % (m[block[b][0]], m[block[b][3]]),
	print 'E_%s_%s,' % (m[block[b][2]], m[block[b][5]]),
	print 'E_%s_%s,' % (m[block[b][4]], m[block[b][7]]),
	print 'E_%s_%s},\n\t' % (m[block[b][6]], m[block[b][7]]),
	print '%+i, %+i, %+i },' % (mm(block[b][3], block[b][4]), 
		mm(block[b][1], block[b][2]), mm(block[b][0], block[b][1]))
	n += 1
*/
#endif

enum Edges { E_NEU_NWU, E_SWU_NWU, E_SEU_SED, E_NED_NWD, E_SEU_SWU, 
	E_NWU_NEU, E_SED_SWD, E_NWU_SWU, E_NED_NEU, E_NEU_NED, 
	E_SWU_SEU, E_SWD_SED, E_NED_SED, E_NWU_NWD, E_NWD_NED, 
	E_SED_SEU, E_NEU_SEU, E_NWD_NWU, E_SWU_SWD, E_SWD_SWU, 
	E_SED_NED, E_SEU_NEU, E_NWD_SWD, E_SWD_NWD };

typedef struct {
	int child_order[8];
	int sub_path[5];
	int b0, b1, b2;
} HilbertPath;	

static HilbertPath hilbert_path[] =
{
	{ /* 0 E_NEU_NWU: */
        { NEU, NED, SED, SEU, SWU, SWD, NWD, NWU },
        { E_NEU_NED, E_NEU_SEU, E_SED_SWD, E_SWU_NWU, E_NWD_NWU},
        -1, -2, -3 },
	
	{ /* 1 E_SWU_NWU: */
        { SWU, SEU, SED, SWD, NWD, NED, NEU, NWU },
        { E_SWU_SEU, E_SWU_SWD, E_SED_NED, E_NWD_NWU, E_NEU_NWU},
        +2, -3, +1 },
	
	{ /* 2 E_SEU_SED: */
        { SEU, SWU, NWU, NEU, NED, NWD, SWD, SED },
        { E_SEU_SWU, E_SEU_NEU, E_NWU_NWD, E_NED_SED, E_SWD_SED},
        -3, +2, -1 },
	
	{ /* 3 E_NED_NWD: */
        { NED, SED, SEU, NEU, NWU, SWU, SWD, NWD },
        { E_NED_SED, E_NED_NEU, E_SEU_SWU, E_NWU_NWD, E_SWD_NWD},
        -1, +3, -2 },
	
	{ /* 4 E_SEU_SWU: */
        { SEU, NEU, NED, SED, SWD, NWD, NWU, SWU },
        { E_SEU_NEU, E_SEU_SED, E_NED_NWD, E_SWD_SWU, E_NWU_SWU},
        -1, -3, +2 },
	
	{ /* 5 E_NWU_NEU: */
        { NWU, SWU, SWD, NWD, NED, SED, SEU, NEU },
        { E_NWU_SWU, E_NWU_NWD, E_SWD_SED, E_NED_NEU, E_SEU_NEU},
        +1, -3, -2 },
	
	{ /* 6 E_SED_SWD: */
        { SED, SEU, NEU, NED, NWD, NWU, SWU, SWD },
        { E_SED_SEU, E_SED_NED, E_NEU_NWU, E_NWD_SWD, E_SWU_SWD},
        -1, +2, +3 },
	
	{ /* 7 E_NWU_SWU: */
        { NWU, NWD, NED, NEU, SEU, SED, SWD, SWU },
        { E_NWU_NWD, E_NWU_NEU, E_NED_SED, E_SEU_SWU, E_SWD_SWU},
        -2, +1, -3 },
	
	{ /* 8 E_NED_NEU: */
        { NED, NWD, SWD, SED, SEU, SWU, NWU, NEU },
        { E_NED_NWD, E_NED_SED, E_SWD_SWU, E_SEU_NEU, E_NWU_NEU},
        +3, -2, -1 },
	
	{ /* 9 E_NEU_NED: */
        { NEU, SEU, SWU, NWU, NWD, SWD, SED, NED },
        { E_NEU_SEU, E_NEU_NWU, E_SWU_SWD, E_NWD_NED, E_SED_NED},
        -3, -1, -2 },
	
	{ /* 10 E_SWU_SEU: */
        { SWU, SWD, NWD, NWU, NEU, NED, SED, SEU },
        { E_SWU_SWD, E_SWU_NWU, E_NWD_NED, E_NEU_SEU, E_SED_SEU},
        +1, +2, -3 },
	
	{ /* 11 E_SWD_SED: */
        { SWD, NWD, NWU, SWU, SEU, NEU, NED, SED },
        { E_SWD_NWD, E_SWD_SWU, E_NWU_NEU, E_SEU_SED, E_NED_SED},
        +1, +3, +2 },
	
	{ /* 12 E_NED_SED: */
        { NED, NEU, NWU, NWD, SWD, SWU, SEU, SED },
        { E_NED_NEU, E_NED_NWD, E_NWU_SWU, E_SWD_SED, E_SEU_SED},
        -2, -1, +3 },
	
	{ /* 13 E_NWU_NWD: */
        { NWU, NEU, SEU, SWU, SWD, SED, NED, NWD },
        { E_NWU_NEU, E_NWU_SWU, E_SEU_SED, E_SWD_NWD, E_NED_NWD},
        -3, -2, +1 },
	
	{ /* 14 E_NWD_NED: */
        { NWD, NWU, SWU, SWD, SED, SEU, NEU, NED },
        { E_NWD_NWU, E_NWD_SWD, E_SWU_SEU, E_SED_NED, E_NEU_NED},
        +1, -2, +3 },
	
	{ /* 15 E_SED_SEU: */
        { SED, NED, NWD, SWD, SWU, NWU, NEU, SEU },
        { E_SED_NED, E_SED_SWD, E_NWD_NWU, E_SWU_SEU, E_NEU_SEU},
        +3, -1, +2 },
	
	{ /* 16 E_NEU_SEU: */
        { NEU, NWU, NWD, NED, SED, SWD, SWU, SEU },
        { E_NEU_NWU, E_NEU_NED, E_NWD_SWD, E_SED_SEU, E_SWU_SEU},
        -2, -3, -1 },
	
	{ /* 17 E_NWD_NWU: */
        { NWD, SWD, SED, NED, NEU, SEU, SWU, NWU },
        { E_NWD_SWD, E_NWD_NED, E_SED_SEU, E_NEU_NWU, E_SWU_NWU},
        +3, +1, -2 },
	
	{ /* 18 E_SWU_SWD: */
        { SWU, NWU, NEU, SEU, SED, NED, NWD, SWD },
        { E_SWU_NWU, E_SWU_SEU, E_NEU_NED, E_SED_SWD, E_NWD_SWD},
        -3, +1, +2 },
	
	{ /* 19 E_SWD_SWU: */
        { SWD, SED, NED, NWD, NWU, NEU, SEU, SWU },
        { E_SWD_SED, E_SWD_NWD, E_NED_NEU, E_NWU_SWU, E_SEU_SWU},
        +3, +2, +1 },
	
	{ /* 20 E_SED_NED: */
        { SED, SWD, SWU, SEU, NEU, NWU, NWD, NED },
        { E_SED_SWD, E_SED_SEU, E_SWU_NWU, E_NEU_NED, E_NWD_NED},
        +2, +3, -1 },
	
	{ /* 21 E_SEU_NEU: */
        { SEU, SED, SWD, SWU, NWU, NWD, NED, NEU },
        { E_SEU_SED, E_SEU_SWU, E_SWD_NWD, E_NWU_NEU, E_NED_NEU},
        +2, -1, -3 },
	
	{ /* 22 E_NWD_SWD: */
        { NWD, NED, NEU, NWU, SWU, SEU, SED, SWD },
        { E_NWD_NED, E_NWD_NWU, E_NEU_SEU, E_SWU_SWD, E_SED_SWD},
        -2, +3, +1 },
	
	{ /* 23 E_SWD_NWD: */
        { SWD, SWU, SEU, SED, NED, NEU, NWU, NWD },
        { E_SWD_SWU, E_SWD_SED, E_SEU_NEU, E_NED_NWD, E_NWU_NWD},
        +2, +1, +3 }
};	

typedef struct {
	_FLOAT_ dx, dy, dz;
} ChildDirection;

static ChildDirection child_direction[] = {
   /* SWD */   { -1, -1, -1},
   /* NWD */   { -1, +1, -1},
   /* SED */   { +1, -1, -1},
   /* NED */   { +1, +1, -1},
   /* SWU */   { -1, -1, +1},
   /* NWU */   { -1, +1, +1},
   /* SEU */   { +1, -1, +1},
   /* NEU */   { +1, +1, +1}
};

static int
createChildsHilbert(FmmvHandle *FMMV, Box *parent, int level, enum KindOfBox kindOfBox, int path_no)
{
   Box * child; 
   int err = 0;
   int n0, n1, n2, n3, n4, n5, n6, n7, n8;
   int c;
   _FLOAT_ x,y,z;
   _FLOAT_ x1, y1, z1;
   _FLOAT_ h;
   HilbertPath *path = hilbert_path + path_no;

   if (level==FMMV->maxLevel) return 0;

   x = parent->x;
   y = parent->y;
   z = parent->z;
   h = ldexp(1.0, -(level+2)) / FMMV->scale;

   if (kindOfBox & SOURCE) {
       if (parent->noOfParticles < FMMV->splitThreshold) {
           return 0;
       }
       n0 = parent->firstParticle;
       n8 = n0 + parent->noOfParticles;

       n4 = radixSplit(FMMV->particles0, FMMV->perm, abs(path->b0)-1, n0, n8, level, path->b0 < 0);
    
       n2 = radixSplit(FMMV->particles0, FMMV->perm, abs(path->b1)-1, n0, n4, level, path->b1 < 0);
       n1 = radixSplit(FMMV->particles0, FMMV->perm, abs(path->b2)-1, n0, n2, level, path->b2 < 0);
       n3 = radixSplit(FMMV->particles0, FMMV->perm, abs(path->b2)-1, n2, n4, level, path->b2 > 0);

       n6 = radixSplit(FMMV->particles0, FMMV->perm, abs(path->b1)-1, n4, n8, level, path->b1 > 0);
       n5 = radixSplit(FMMV->particles0, FMMV->perm, abs(path->b2)-1, n4, n6, level, path->b2 < 0);
       n7 = radixSplit(FMMV->particles0, FMMV->perm, abs(path->b2)-1, n6, n8, level, path->b2 > 0);   
   }
   else { /* kindOfBox == TARGET */
       if (parent->noOfTargets < FMMV->splitTargetThreshold) {
           return 0;
       }
       n0 = parent->firstTarget;
       n8 = n0 + parent->noOfTargets;

       n4 = radixSplit(FMMV->targets0, FMMV->permTargets, abs(path->b0)-1, n0, n8, level, path->b0 < 0);

       n2 = radixSplit(FMMV->targets0, FMMV->permTargets, abs(path->b1)-1, n0, n4, level, path->b1 < 0);
       n1 = radixSplit(FMMV->targets0, FMMV->permTargets, abs(path->b2)-1, n0, n2, level, path->b2 < 0);
       n3 = radixSplit(FMMV->targets0, FMMV->permTargets, abs(path->b2)-1, n2, n4, level, path->b2 > 0);
    
       n6 = radixSplit(FMMV->targets0, FMMV->permTargets, abs(path->b1)-1, n4, n8, level, path->b1 > 0);
       n5 = radixSplit(FMMV->targets0, FMMV->permTargets, abs(path->b2)-1, n4, n6, level, path->b2 < 0);
       n7 = radixSplit(FMMV->targets0, FMMV->permTargets, abs(path->b2)-1, n6, n8, level, path->b2 > 0);   
   }	   

   if (n1-n0>0) {
       c = path->child_order[0];
       x1 = x + h*child_direction[c].dx;
       y1 = y + h*child_direction[c].dy;
       z1 = z + h*child_direction[c].dz;
       child = createBox(FMMV, parent, x1, y1, z1, n0, n1-n0, level+1, c, kindOfBox);
       if (!child) return -1;
       parent->child[c] = child;
       err = createChildsHilbert(FMMV, child, level+1, kindOfBox, path->sub_path[0]);
       if (err) return -1;			  
   }

   if (n2-n1>0) {
       c = path->child_order[1];
       x1 = x + h*child_direction[c].dx;
       y1 = y + h*child_direction[c].dy;
       z1 = z + h*child_direction[c].dz;
       child = createBox(FMMV, parent, x1, y1, z1, n1, n2-n1, level+1, c, kindOfBox);
       if (!child) return -1;
       parent->child[c] = child;
       err = createChildsHilbert(FMMV, child, level+1, kindOfBox, path->sub_path[1]);
       if (err) return -1;			  
   }   

   if (n3-n2>0) {
       c = path->child_order[2];
       x1 = x + h*child_direction[c].dx;
       y1 = y + h*child_direction[c].dy;
       z1 = z + h*child_direction[c].dz;
       child = createBox(FMMV, parent, x1, y1, z1, n2, n3-n2, level+1, c, kindOfBox);
       if (!child) return -1;
       parent->child[c] = child;
       err = createChildsHilbert(FMMV, child, level+1, kindOfBox, path->sub_path[1]);
       if (err) return -1;			  
   }   
   
   if (n4-n3>0) {
       c = path->child_order[3];
       x1 = x + h*child_direction[c].dx;
       y1 = y + h*child_direction[c].dy;
       z1 = z + h*child_direction[c].dz;
       child = createBox(FMMV, parent, x1, y1, z1, n3, n4-n3, level+1, c, kindOfBox);
       if (!child) return -1;
       parent->child[c] = child;
       err = createChildsHilbert(FMMV, child, level+1, kindOfBox, path->sub_path[2]);
       if (err) return -1;			  
   }      

   if (n5-n4>0) {
       c = path->child_order[4];
       x1 = x + h*child_direction[c].dx;
       y1 = y + h*child_direction[c].dy;
       z1 = z + h*child_direction[c].dz;
       child = createBox(FMMV, parent, x1, y1, z1, n4, n5-n4, level+1, c, kindOfBox);
       if (!child) return -1;
       parent->child[c] = child;
       err = createChildsHilbert(FMMV, child, level+1, kindOfBox, path->sub_path[2]);
       if (err) return -1;			  
   }

   if (n6-n5>0) {
       c = path->child_order[5];
       x1 = x + h*child_direction[c].dx;
       y1 = y + h*child_direction[c].dy;
       z1 = z + h*child_direction[c].dz;
       child = createBox(FMMV, parent, x1, y1, z1, n5, n6-n5, level+1, c, kindOfBox);
       if (!child) return -1;
       parent->child[c] = child;
       err = createChildsHilbert(FMMV, child, level+1, kindOfBox, path->sub_path[3]);
       if (err) return -1;			  
   }   

   if (n7-n6>0) {
       c = path->child_order[6];
       x1 = x + h*child_direction[c].dx;
       y1 = y + h*child_direction[c].dy;
       z1 = z + h*child_direction[c].dz;
       child = createBox(FMMV, parent, x1, y1, z1, n6, n7-n6, level+1, c, kindOfBox);
       if (!child) return -1;
       parent->child[c] = child;
       err = createChildsHilbert(FMMV, child, level+1, kindOfBox, path->sub_path[3]);
       if (err) return -1;			  
   }   
   
   if (n8-n7>0) {
       c = path->child_order[7];
       x1 = x + h*child_direction[c].dx;
       y1 = y + h*child_direction[c].dy;
       z1 = z + h*child_direction[c].dz;
       child = createBox(FMMV, parent, x1, y1, z1, n7, n8-n7, level+1, c, kindOfBox);
       if (!child) return -1;
       parent->child[c] = child;
       err = createChildsHilbert(FMMV, child, level+1, kindOfBox, path->sub_path[4]);
       if (err) return -1;			  
   }         
   
   return 0;
   
}

#endif /* FMM_DIM==3 */
/********************************************************************************/



int buildTree(FmmvHandle *FMMV)
{
    int err;
    int i;
    Box *root;
    _FLOAT_ m;
    _FLOAT_ scale;

    for (i=0; i<=FMMV->maxLevel; i++) {
	FMMV->firstSourceBoxOfLevel[i] = 0;
	FMMV->firstTargetBoxOfLevel[i] = 0;
    }
    
    FMMV->particles0 = (_FLOAT_ (*)[FMM_DIM]) FMMV_MALLOC(FMMV, FMM_DIM*FMMV->NParticles*sizeof(_FLOAT_));
    if (FMMV->particles0==0) return -1;

    scale=0.5*FMMV->scale; /* TODO: comment factor 0.5 etc.*/
    for (i=0; i<FMMV->NParticles; i++) {
        FMMV->perm[i] = i;
	FMMV->particles0[i][0] = scale * FMMV->particles[i][0]+0.5;
	FMMV->particles0[i][1] = scale * FMMV->particles[i][1]+0.5;
	#if (FMM_DIM>=3)
	FMMV->particles0[i][2] = scale * FMMV->particles[i][2]+0.5;
	#endif
    }
    
    m = 0.5/FMMV->scale;	
    root = createBox(FMMV, 0, m, m, 
		    #if (FMM_DIM>=3)
		    m, 
	            #endif
		    0, FMMV->NParticles, 0, 0, SOURCE|TARGET);
    if (!root)  return -1;

    root->firstTarget = 0;
    root->noOfTargets = FMMV->NTargets;
    root->nextTargetBox = 0;
    FMMV->firstTargetBoxOfLevel[0] = root;

    if (FMMV->useHilbertOrder) {
	    err = createChildsHilbert(FMMV, root, 0, SOURCE|TARGET, 0);
    }	    
    else {
	    err = createChilds(FMMV, root, 0, SOURCE|TARGET);
    }	    
    if (err) return -1;

    genNeighborRelations(root); 
    if (FMMV->periodicBoundaryConditions) {
        genPeriodicNeighborRelations(root);
    }	
    
    FMMV_FREE(FMMV, FMMV->particles0, FMM_DIM*FMMV->NParticles*sizeof(_FLOAT_));
    return 0;
}


int buildTree_ST(FmmvHandle *FMMV)
{
    int err;
    int i;
    Box *root;
    _FLOAT_ m;
    _FLOAT_ scale;

    for (i=0; i<=FMMV->maxLevel; i++) {
	FMMV->firstSourceBoxOfLevel[i] = 0;
	FMMV->firstTargetBoxOfLevel[i] = 0;
    }

    FMMV->particles0 = (_FLOAT_ (*)[FMM_DIM]) FMMV_MALLOC(FMMV, FMM_DIM*FMMV->NParticles*sizeof(_FLOAT_));
    if (FMMV->particles0==0) return -1;

    FMMV->targets0 = (_FLOAT_ (*)[FMM_DIM]) FMMV_MALLOC(FMMV, FMM_DIM*FMMV->NTargets*sizeof(_FLOAT_));
    if (FMMV->targets0==0) return -1;

    scale = 0.5*FMMV->scale;
    for (i=0; i<FMMV->NParticles; i++) {
        FMMV->perm[i] = i;
	FMMV->particles0[i][0] = scale * FMMV->particles[i][0]+0.5;
	FMMV->particles0[i][1] = scale * FMMV->particles[i][1]+0.5;
	#if (FMM_DIM>=3)
	FMMV->particles0[i][2] = scale * FMMV->particles[i][2]+0.5;
	#endif
    }

    for (i=0; i<FMMV->NTargets; i++) {
        FMMV->permTargets[i] = i;
	FMMV->targets0[i][0] = scale * FMMV->targets[i][0]+0.5;
	FMMV->targets0[i][1] = scale * FMMV->targets[i][1]+0.5;
	#if (FMM_DIM>=3)
	FMMV->targets0[i][2] = scale * FMMV->targets[i][2]+0.5;
	#endif
    }
    
    m = 0.5/FMMV->scale;	
    root = createBox(FMMV, 0, m, m, 
		    #if (FMM_DIM>=3)
		    m, 
	            #endif
		    0, FMMV->NParticles, 0, 0, SOURCE|TARGET);
    if (!root) return -1;

    root->firstTarget = 0;
    root->noOfTargets = FMMV->NTargets;
    root->nextTargetBox = 0;
    FMMV->firstTargetBoxOfLevel[0] = root;

    if (FMMV->useHilbertOrder) {
	    err = createChildsHilbert(FMMV, root, 0, SOURCE, 0);
    }	    
    else {
    	    err = createChilds(FMMV, root, 0, SOURCE);
    }	    
    if (err) return -1;

    if (FMMV->useHilbertOrder) {
	    err = createChildsHilbert(FMMV, root, 0, TARGET, 0);
    }
    else {	    
	    err = createChilds(FMMV, root, 0, TARGET);
    }	    
    if (err) return -1;
    
    genNeighborRelations(root); 
    if (FMMV->periodicBoundaryConditions) {
        genPeriodicNeighborRelations(root);
    }	    
    
    FMMV_FREE(FMMV, FMMV->particles0, FMM_DIM*FMMV->NParticles*sizeof(_FLOAT_));
    FMMV_FREE(FMMV, FMMV->targets0, FMM_DIM*FMMV->NTargets*sizeof(_FLOAT_));
    return 0;
}

void freeTree(FmmvHandle *FMMV, Box *box)
{
    if (!box) return;

#ifdef USE_PTHREADS
    if(FMMV->useFarfieldNearfieldThreads)  {
    	pthread_mutex_destroy(&(box->mutex));
    }	
#endif    
    freeTree(FMMV, box->child[0]);
    freeTree(FMMV, box->child[1]);
    freeTree(FMMV, box->child[2]);
    freeTree(FMMV, box->child[3]);
    #if (FMM_DIM>=3)
    freeTree(FMMV, box->child[4]);
    freeTree(FMMV, box->child[5]);
    freeTree(FMMV, box->child[6]);
    freeTree(FMMV, box->child[7]);
    #endif
    
    FMMV_FREE(FMMV, box, sizeof(Box));
}	

void freeTree_ST(FmmvHandle *FMMV, Box *box)
{
    if (!box) return;

    if(isTarget(box)) {
#ifdef USE_PTHREADS
   	if(FMMV->useFarfieldNearfieldThreads)  {
		pthread_mutex_destroy(&(box->mutex));
  	}	
#endif
    }	    
    freeTree_ST(FMMV, box->child[0]);
    freeTree_ST(FMMV, box->child[1]);
    freeTree_ST(FMMV, box->child[2]);
    freeTree_ST(FMMV, box->child[3]);
    #if (FMM_DIM>=3)
    freeTree_ST(FMMV, box->child[4]);
    freeTree_ST(FMMV, box->child[5]);
    freeTree_ST(FMMV, box->child[6]);
    freeTree_ST(FMMV, box->child[7]);
    #endif
    
    FMMV_FREE(FMMV, box, sizeof(Box));
}	


