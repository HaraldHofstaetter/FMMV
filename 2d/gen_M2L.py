N   = {    0 : "0",
	  .5 : "1",
	 -.5 : "m1",
	 1.0 : "2",
	-1.0 : "m2",
	 1.5 : "3",
	-1.5 : "m3",
	 2.0 : "4",
	-2.0 : "m4",
	 2.5 : "5",
	-2.5 : "m5",
	 3.0 : "6",
	-3.0 : "m6",	
	 3.5 : "7",
	-3.5 : "m7",	
	 4.0 : "8",
	-4.0 : "m8",	
	 4.5 : "9",
	-4.5 : "m9"}	

def gen_enum():
	X = [
	('N', +0.5, +0.5),
	('N', -0.5, +0.5),
	('N', +0.5, -0.5),
	('N', -0.5, -0.5),
	('S', +0.5, +0.5),
	('S', -0.5, +0.5),
	('S', +0.5, -0.5),
	('S', -0.5, -0.5),
	('W', +0.5, +0.5),
	('W', -0.5, +0.5),
	('W', +0.5, -0.5),
	('W', -0.5, -0.5),
	('E', +0.5, +0.5),
	('E', -0.5, +0.5),
	('E', +0.5, -0.5),
	('E', -0.5, -0.5),
	
	('N', -0.5, +2.5),
	('N', +0.5, +2.5),
	('S', -0.5, -2.5),
	('S', +0.5, -2.5),
	('W', -2.5, -0.5), 
	('W', -2.5, +0.5),
	('E', +2.5, -0.5),
	('E', +2.5, +0.5),

	('W', -2.5, -1.5),
	('S', -1.5, -2.5),
	('S', -2.5, -2.5),
	('W', -2.5, +1.5),
	('N', -1.5, +2.5),
	('W', -2.5, +2.5),
	('N', +1.5, +2.5),
	('E', +2.5, +1.5),
	('N', +2.5, +2.5),
	('S', +1.5, -2.5),
	('E', +2.5, -1.5),
	('E', +2.5, -2.5),

	('N', -1,  2),
	('N',  0,  2),
	('N',  1,  2),
	('N',  2,  2),
	
	('E',  2,  1),
	('E',  2,  0),
	('E',  2, -1),
	('E',  2, -2),
	
	('S',  1, -2),
	('S',  0, -2),
	('S', -1, -2),
	('S', -2, -2),
	
	('W', -2, -1),
	('W', -2,  0),
	('W', -2,  1),
	('W', -2,  2)
	]

	n = 0
	print "enum {"	
	k=0
	for (dir, re, im) in X:
		n += 1
		print "D_%s_%s_%s," % (dir, N[re], N[im]),
		k += 1
		if k==4:
			k=0
			print
	print "};"
	#print "#define N_D_X2X %i" % n
	#	
	#print
	#for (dir, re, im) in X:
	#	print "\t{D_%s_%s_%s, X%s, %.1f, %.1f}," % (dir, N[re], N[im], dir, re, im)
	#print "};"	

	








	
def gen_X2X(): 
	print '#include "fmmv2d.h"'
	print '#include "fmmv.h"'
	print '#include <string.h> /* memset */'
	print '_FLOAT_ *VEC_addmul_c(FmmvHandle *FMMV, int n, _FLOAT_ *a, _FLOAT_ *out, _FLOAT_ *in);'
	print 'void gen_M2L_interaction_list(FmmvHandle *FMMV, Box *box, Box **dummy);'
	print 'void prepare_X2L(FmmvHandle *FMMV, Box *box);'
	print 
	print 'void M2L(FmmvHandle *FMMV, Box *box)' 
	print '{'	
	print "\t_FLOAT_ *D_X2X = FMMV->D_X2X;"
	print "\tBox *nc;"
	print "\tint s_exp = FMMV->s_exp;"
	print "\tint s_exp2 = 2*FMMV->s_exp;"
	print "\tBox *IL[32]; /*6^2-2^2*/"


	#print '\t_FLOAT_ *outE, *outN, *outW, *outS;'
	for i in [ 'E', 'N', 'W', 'S']:
		print "\t_FLOAT_ out%s[2*FMM_S_EPS_MAX];" % i
	for i in [ 'E', 'N', 'W', 'S']:
		print "\t_FLOAT_ s_%s[2*FMM_S_EPS_MAX];" % i
	print
	print '\tif (!isSource(box) || !hasSourceChilds(box)) return;'
	print
	print '\tgen_M2L_interaction_list(FMMV, box, IL);'
	print

	for i in [ 'E', 'N', 'W', 'S']:
		print "\tmemset(s_%s, 0, s_exp2*sizeof(_FLOAT_));" %i

	print
		
	for (child, w_re, w_im, direct_inter) in [
		('SW',  0.5,  0.5,
		       [("N", "NW", "SE", -1,  2),
			("N", "NM", "SW",  0,  2),
			("N", "NM", "SE",  1,  2),
			("N", "NE", "SW",  2,  2),
			("E", "ME" ,"NW",  2,  1),	
			("E", "ME" ,"SW",  2,  0),	
			("E", "SE" ,"NW",  2, -1)]),
		('NW',  0.5, -0.5,
		       [("E", "NE", "SW",  2,  1),
			("E", "ME", "NW",  2,  0),
			("E", "ME", "SW",  2, -1),
			("E", "SE", "NW",  2, -2),#
			("S", "SM", "NE",  1, -2),
			("S", "SM", "NW",  0, -2),
			("S", "SW", "NE", -1, -2)]),
		('NE', -0.5, -0.5, 
		       [("S", "SE", "NW",  1, -2), #
			("S", "SM", "NE",  0, -2),
			("S", "SM", "NW", -1, -2),
			("S", "SW", "NE", -2, -2),
			("W", "MW", "SE", -2, -1),
			("W", "MW", "NE", -2,  0),
			("W", "NW", "SE", -2,  1)]),
		('SE', -0.5,  0.5,
		       [("W", "SW", "NE", -2, -1),
			("W", "MW", "SE", -2,  0),
			("W", "MW", "NE", -2,  1),
			("W", "NW", "SE", -2,  2),#
			("N", "NM", "SW", -1,  2),
			("N", "NM", "SE",  0,  2),
			("N", "NE", "SW",  1,  2)])		
		]:
		print '\n\tif (box->child[%s]) {' % child
		print '\t\tM2X(FMMV, box->child[%s], outN, outS, outW, outE);' % child
		for dir in ['E', 'N','W', 'S']:
			print "\t\tVEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_%s_%s_%s, out%s, s_%s);" % \
				(dir, N[w_re], N[w_im], dir, dir)
		oldneighbor = None
        	for (dir, neighbor, neighborchild, w_re, w_im) in direct_inter:
			if oldneighbor!=neighbor:
				if oldneighbor!=None:
					print '\t\t}' # close neighbor
				oldneighbor = neighbor
				print '\t\tif (box->neighbor[%s_]) {' % neighbor
			print '\t\t\tnc = box->neighbor[%s_]->child[%s];' % (neighbor, neighborchild)
			print '\t\t\tif (nc) {' 
			print "\t\t\t\tnc->X[X%s] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_%s_%s_%s, out%s, nc->X[X%s]);" % \
			(dir, dir, N[w_re], N[w_im], dir, dir)
			print '\t\t\t}'
		print '\t\t}' # close neighbor
		print '\t}' #close child

	for (dir, nb, child1, w_re1, w_im1, child2, w_re2, w_im2) in [
	    ('N', 'NM_', 'NW', -0.5, +2.5, 'NE', +0.5, +2.5),
	    ('S', 'SM_', 'SW', -0.5, -2.5, 'SE', +0.5, -2.5),
	    ('W', 'MW_', 'SW', -2.5, -0.5, 'NW', -2.5, +0.5),
	    ('E', 'ME_', 'SE', +2.5, -0.5, 'NE', +2.5, +0.5),
	    ]:
		print '\n\tif (box->neighbor[%s]) {' % nb
		print '\t\tnc = box->neighbor[%s]->child[%s];' % (nb, child1)
		print '\t\tif (nc) {' 
		print "\t\t\tnc->X[X%s] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_%s_%s_%s, s_%s, nc->X[X%s]);" % \
			(dir, dir, N[w_re1], N[w_im1], dir, dir)
		print '\t\t}'				
		print '\t\tnc = box->neighbor[%s]->child[%s];' % (nb, child2)
		print '\t\tif (nc) {' 
		print "\t\t\tnc->X[X%s] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_%s_%s_%s, s_%s, nc->X[X%s]);" % \
			(dir, dir, N[w_re2], N[w_im2], dir, dir)
		print '\t\t}'				
		print '\t}'			
			

	for (nb, child1, dir1, w_re1, w_im1, 
                 child2, dir2, w_re2, w_im2,
                 child3, dir3, w_re3, w_im3) in [
	    ('SW_', 'NW', 'W', -2.5, -1.5,
	            'SE', 'S', -1.5, -2.5,
		    'SW', 'S', -2.5, -2.5),
	    ('NW_', 'SW', 'W', -2.5, 1.5,
	            'NE', 'N', -1.5, 2.5,
		    'NW', 'W', -2.5, 2.5),
	    ('NE_', 'NW', 'N', +1.5, 2.5,
	            'SE', 'E', +2.5, 1.5,
		    'NE', 'N', +2.5, 2.5),
	    ('SE_', 'SW', 'S', +1.5, -2.5,
	            'NE', 'E', +2.5, -1.5,
		    'SE', 'E', +2.5, -2.5)
            ]:
		print '\n\tif (box->neighbor[%s]) {' % nb
		print '\t\tnc = box->neighbor[%s]->child[%s];' % (nb, child1)
		print '\t\tif (nc) {' 
		print "\t\t\tnc->X[X%s] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_%s_%s_%s, s_%s, nc->X[X%s]);" % \
			(dir1, dir1, N[w_re1], N[w_im1], dir1, dir1)
		print '\t\t}'	
		print '\t\tnc = box->neighbor[%s]->child[%s];' % (nb, child2)
		print '\t\tif (nc) {' 
		print "\t\t\tnc->X[X%s] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_%s_%s_%s, s_%s, nc->X[X%s]);" % \
			(dir2, dir2, N[w_re2], N[w_im2], dir2, dir2)
		print '\t\t}'	
		print '\t\tnc = box->neighbor[%s]->child[%s];' % (nb, child3)
		print '\t\tif (nc) {' 
		print "\t\t\tnc->X[X%s] = VEC_addmul_c(FMMV, s_exp, D_X2X + s_exp2*D_%s_%s_%s, s_%s, nc->X[X%s]);" % \
			(dir3, dir3, N[w_re3], N[w_im3], dir3, dir3)
		print '\t\t}'	
		print '\t}'
		
	print '\tprepare_X2L(FMMV, box);'
	print '}' #close M2L

	
gen_enum()
gen_X2X()	
print """
/*** dummy routines yet... ***/

#include<assert.h>

void M2L_ws2(FmmvHandle *FMMV, Box *box)
{
	assert(0);
}

void M2L_ws2_reduced(FmmvHandle *FMMV, Box *box)
{
	assert(0);
}
"""
	
