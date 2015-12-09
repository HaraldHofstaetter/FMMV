# FMMV - the Fastest Multipole Method of Vienna
# Copyright (c) 2006-2015 Harald Hofstaetter
# http://www.harald-hofstaetter.at
#
# This file is part of FMMV.
#
# FMMV is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# FMMV is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with FMMV; if not, write to the Free Software  Foundation, Inc.,
# 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#

def gen_core_gen_M_L(funname, p_max):
	type = Op.templates["type"] 	
	basetype = Op.templates["basetype"]
	
	begin_block('void %s(int p, %s *scale, %s q, %s *dML, %s *ML)\n' % (funname, basetype, type, basetype, basetype))
	scale = parameter("scale", array=True)
	q = parameter("q")
	dML = parameter("dML", array=True)
	ML = parameter("ML", array=True, horiz_add=-1)
	f = var("f")

	print
	begin_block("switch(p) ")
	for n in range(0, p_max+1):
		print "#if FMM_P_MAX >= %i" % n
	for n in range(p_max, -1, -1):
		print "\tcase %i:" % n
		#print "\t\tf = scale[%i];" % n
                f ^= scale[n]*q;
		for m in range((n+1)*(n+2)-1, n*(n+1)-1, -1):
			ML[m] += f * dML[m]
		print "#endif /* FMM_P_MAX >= %i */" % n
	end_block()

	end_block()
	print


def J(n,m): 
	return (n*(n+1))/2 + m

def Re(n,m): 
	return 2*((n*(n+1))/2 + m)

def Im(n,m): 
	return 2*((n*(n+1))/2 + m) + 1


def gen_core_eval_L_M(funname, p_max):
	type = Op.templates["type"] 	
	basetype = Op.templates["basetype"]
	
	begin_block('%s %s(int p, %s *scale, %s *LM, %s *Y)\n' % (type, funname, basetype, basetype, basetype))
	LM = parameter("LM", array=True, base_array=True)
	Y = parameter("Y", array=True)
	scale = parameter("scale", array=True)
	
	h = var("h")
	c = var("c")
	two = var("two", 2.0)
	
	print
	h ^= 0

	print
	begin_block("switch(p) ")
	for j in range(0, p_max+1):
		print "#if FMM_P_MAX >= %i" % j
	for j in range(p_max,-1,-1):
		print "\tcase %i:" % j
		if j==0:
			c ^= LM[Re(j,0)]*Y[Re(j,0)]
		else:
			c ^= LM[Re(j,1)]*Y[Re(j,1)]-LM[Im(j,1)]*Y[Im(j,1)]
			for k in range(2,j+1):
				c +=  LM[Re(j,k)]*Y[Re(j,k)] - LM[Im(j,k)]*Y[Im(j,k)]
			c ^= two*c + LM[Re(j,0)]*Y[Re(j,0)]
		#h ^= c + scale*h
		h += scale[j]*c 
		print "#endif /* FMM_P_MAX >= %i */" % j
	end_block()
	print '\treturn h;'		
	end_block()
	print


from math import sqrt


def gen_core_eval_L_M_grad_minus(funname, p_max):
	type = Op.templates["type"] 	
	basetype = Op.templates["basetype"]
	
	begin_block('void %s(int p, %s *scale, %s *L, %s *Y, %s *x, %s *y, %s *z)\n' 
                 % (funname, basetype, basetype, basetype, type, type, type))
	L = parameter("L", array=True, base_array=True)
	Y = parameter("Y", array=True)
	scale = parameter("scale", array=True)
	
	hx = var("hx")
	hy = var("hy")
	hz = var("hz")
	cx = var("cx")
	cy = var("cy")
	cz = var("cz")
	f1 = var("f1")
	f2 = var("f2")
	
	print
	hx ^= 0
	hy ^= 0
	hz ^= 0
	
	print
	begin_block("switch(p) ")
	for j in range(1, p_max+1):
		print "#if FMM_P_MAX >= %i" % j
	for j in range(p_max,0,-1):
		print "\tcase %i:" % j
		cx ^= sqrt(float(j*(j-1))) * L[Re(j,0)] * Y[Re(j-1,1)]
		cy ^= sqrt(float(j*(j-1))) * L[Re(j,0)] * Y[Im(j-1,1)]
		cz ^= j*L[Re(j,0)] * Y[Re(j-1,0)]

		for k in range(1,j+1):
			if k<=j-2:
				f1 ^= sqrt(float((j-k)*(j-k-1)))
				f2 ^= sqrt(float((j+k)*(j+k-1)))
				cx ^= cx + f1*(L[Re(j,k)]*Y[Re(j-1,k+1)] - L[Im(j,k)]*Y[Im(j-1,k+1)])\
				     + f2*(L[Im(j,k)]*Y[Im(j-1,k-1)] - L[Re(j,k)]*Y[Re(j-1,k-1)])
				cy ^= cy + f1*(L[Im(j,k)]*Y[Re(j-1,k+1)] + L[Re(j,k)]*Y[Im(j-1,k+1)])\
				     + f2*(L[Im(j,k)]*Y[Re(j-1,k-1)] + L[Re(j,k)]*Y[Im(j-1,k-1)])
			else:
				cx ^= cx - sqrt(float((j+k)*(j+k-1)))*(L[Re(j,k)]*Y[Re(j-1,k-1)] - L[Im(j,k)]*Y[Im(j-1,k-1)])
				cy ^= cy + sqrt(float((j+k)*(j+k-1)))*(L[Im(j,k)]*Y[Re(j-1,k-1)] + L[Re(j,k)]*Y[Im(j-1,k-1)])
				
			cz ^= cz + 2.0*sqrt(float((j+k)*(j-k)))*(L[Re(j,k)]*Y[Re(j-1,k)] - L[Im(j,k)]*Y[Im(j-1,k)])
		hx ^= hx + scale[j]*cx
		hy ^= hy + scale[j]*cy
		hz ^= hz + scale[j]*cz
		print "#endif /* FMM_P_MAX >= %i */" % j
	end_block()

	print '\t*x = hx;'		
	print '\t*y = hy;'		
	print '\t*z = hz;'		

	end_block()
	print


def gen_core_eval_L_M_grad_plus(funname, p_max):
	type = Op.templates["type"] 	
	basetype = Op.templates["basetype"]
	
	begin_block('void %s(int p, %s *scale, %s *M, %s *Y, %s *x, %s *y, %s *z)\n' 
	 	 % (funname, basetype, basetype, basetype, type, type, type))
	M = parameter("M", array=True, base_array=True)
	Y = parameter("Y", array=True)
	scale = parameter("scale", array=True)
	
	hx = var("hx")
	hy = var("hy")
	hz = var("hz")
	cx = var("cx")
	cy = var("cy")
	cz = var("cz")
	f1 = var("f1")
	f2 = var("f2")

	print
	hx ^= 0
	hy ^= 0
	hz ^= 0

	print
	begin_block("switch(p) ")
	for j in range(0, p_max+1):
		print "#if FMM_P_MAX >= %i" % j
	for j in range(p_max,-1,-1):
		print "\tcase %i:" % j
		cx ^= -sqrt(float((j+2)*(j+1)))*M[Re(j,0)]*Y[Re(j+1,1)]
		cy ^= -sqrt(float((j+2)*(j+1)))*M[Re(j,0)]*Y[Im(j+1,1)]
		cz ^= (j+1)*M[Re(j,0)]*Y[Re(j+1,0)]

		for k in range(1,j+1):
			f1 ^= sqrt(float((j+k+2)*(j+k+1)))
			f2 ^= sqrt(float((j-k+2)*(j-k+1)))
			cx ^= cx - f1*(M[Re(j,k)]*Y[Re(j+1,k+1)] - M[Im(j,k)]*Y[Im(j+1,k+1)])\
			      -f2*(M[Im(j,k)]*Y[Im(j+1,k-1)] - M[Re(j,k)]*Y[Re(j+1,k-1)])
			cy ^= cy - f1*(M[Im(j,k)]*Y[Re(j+1,k+1)] + M[Re(j,k)]*Y[Im(j+1,k+1)])\
			      -f2*(M[Im(j,k)]*Y[Re(j+1,k-1)] + M[Re(j,k)]*Y[Im(j+1,k-1)])
		 	cz ^= cz - 2.0*sqrt(float((j+k+1)*(j-k+1)))*(M[Im(j,k)]*Y[Im(j+1,k)] - M[Re(j,k)]*Y[Re(j+1,k)])
		hx ^= hx + scale[j]*cx
		hy ^= hy + scale[j]*cy
		hz ^= hz + scale[j]*cz
		print "#endif /* FMM_P_MAX >= %i */" % j
	end_block()

	print '\t*x = hx;'		
	print '\t*y = hy;'		
	print '\t*z = hz;'		

	end_block()
	print

def gen_core_gen_M_L_dipole_minus(funname, p_max):
	type = Op.templates["type"] 	
	basetype = Op.templates["basetype"]
	
	begin_block('void %s(int p, %s *scale, %s mx, %s my, %s mz, %s *Y, %s *M)\n' % (funname, basetype, type, type, type, basetype, basetype))
	scale = parameter("scale")
	mx = parameter("mx")
	my = parameter("my")
	mz = parameter("mz")
	Y = parameter("Y", array=True)
	M = parameter("M", array=True, horiz_add=-1)
	f = var("f")

	print
	begin_block("switch(p) ")
	for j in range(1, p_max+1):
		print "#if FMM_P_MAX >= %i" % j
	for j in range(p_max, 3-1, -1):
		print "\tcase %i:" % j
		##print "\t\tf = scale[%i];" % j
		f ^= scale[j]

		M[Im(j,j)] += f * ((0.5*sqrt(float(2*j-1)*(2*j)))*(my*Y[Re(j-1,j-1)] - mx*Y[Im(j-1,j-1)]) )
		M[Re(j,j)] += f * ((-0.5*sqrt(float(2*j-1)*(2*j)))*(mx*Y[Re(j-1,j-1)] + my*Y[Im(j-1,j-1)]) )
		if j>1:
			M[Im(j,j-1)] += f * ((0.5*sqrt(float((2*j-2)*(2*j-1))))*(my*Y[Re(j-1,j-2)] - mx*Y[Im(j-1,j-2)]) 
					  + (sqrt(float(2*j-1)))*mz*Y[Im(j-1,j-1)] )
			M[Re(j,j-1)] += f * ((-0.5*sqrt(float((2*j-2)*(2*j-1))))*(mx*Y[Re(j-1,j-2)] + my*Y[Im(j-1,j-2)]) 
					  + (sqrt(float(2*j-1)))*mz*Y[Re(j-1,j-1)] )
		for k in range(j-2+1-1, 2-1, -1): 
			M[Im(j,k)] += f * ( (0.5*sqrt(float((j-k-1)*(j-k))))*(my*Y[Re(j-1,k+1)] + mx*Y[Im(j-1,k+1)]) 
					  + (0.5*sqrt(float((j+k-1)*(j+k))))*(my*Y[Re(j-1,k-1)] - mx*Y[Im(j-1,k-1)]) 
					  + (sqrt(float((j-k)*(j+k))))*mz*Y[Im(j-1,k)] )
			M[Re(j,k)] += f * ( (0.5*sqrt(float((j-k-1)*(j-k))))*(mx*Y[Re(j-1,k+1)] - my*Y[Im(j-1,k+1)]) 
					  - (0.5*sqrt(float((j+k-1)*(j+k))))*(mx*Y[Re(j-1,k-1)] + my*Y[Im(j-1,k-1)]) 
					  + (sqrt(float((j-k)*(j+k))))*mz*Y[Re(j-1,k)] )
		
		M[Im(j,1)] += f * ( (0.5*sqrt(float((j-1-1)*(j-1))))*(my*Y[Re(j-1,1+1)] + mx*Y[Im(j-1,1+1)]) 
					  + (0.5*sqrt(float((j+1-1)*(j+1))))*my*Y[Re(j-1,1-1)] + (sqrt(float((j-1)*(j+1))))*mz*Y[Im(j-1,1)] )
		M[Re(j,1)] += f * ( (0.5*sqrt(float((j-1-1)*(j-1))))*(mx*Y[Re(j-1,1+1)] - my*Y[Im(j-1,1+1)])
					  - (0.5*sqrt(float((j+1-1)*(j+1))))*mx*Y[Re(j-1,1-1)] + (sqrt(float((j-1)*(j+1))))*mz*Y[Re(j-1,1)] )
		M[Im(j,0)] += 0
		M[Re(j,0)] += f * (j*mz*Y[Re(j-1,0)] + sqrt(float((j-1)*j))*(mx*Y[Re(j-1,1)] - my*Y[Im(j-1,1)]))
		print "#endif /* FMM_P_MAX >= %i */" % j

	print "\tcase 2:"
	f ^= scale[2]
	M[11] += f * ( sqrt(float(3))*(my*Y[4] - mx*Y[5]) )
	M[10] += f * ((-sqrt(float(3)))*(mx*Y[4] + my*Y[5]) )
	M[9] += f * ( sqrt(float(1.5))*my*Y[2] + sqrt(float(3))*mz*Y[5] )
	M[8] += f * ((-sqrt(float(1.5)))*mx*Y[2] + sqrt(float(3))*mz*Y[4] )
	M[7] += 0
	M[6] += f * (2.0*mz*Y[2] + sqrt(float(2))*(mx*Y[4]-my*Y[5]))
	print "#endif /* FMM_P_MAX >= 2 */" 
	print "\tcase 1:"
	f ^= scale[1]
	M[5] += f * sqrt(float(0.5))*my*Y[0]
	M[4] += f * (-sqrt(float(0.5)))*mx*Y[0]
	M[3] += 0
	M[2] += f * mz*Y[0]
	print "#endif /* FMM_P_MAX >= 1 */" 
	print "\tcase 0:"
	M[1] += 0
	M[0] += 0

	end_block()
		
	end_block()
	print

def gen_core_gen_M_L_dipole_plus(funname, p_max):
	type = Op.templates["type"] 	
	basetype = Op.templates["basetype"]
	
	begin_block('void %s(int p, %s *scale, %s mx, %s my, %s mz, %s *Y, %s *L)\n' % (funname, basetype, type, type, type, basetype, basetype))
	scale = parameter("scale", array=True)
	mx = parameter("mx")
	my = parameter("my")
	mz = parameter("mz")
	Y = parameter("Y", array=True)
	L = parameter("L", array=True, horiz_add=-1)
	f = var("f")

	print
	begin_block("switch(p) ")
	for j in range(0, p_max+1):
		print "#if FMM_P_MAX >= %i" % j
	for j in range(p_max, -1, -1):
		print "\tcase %i:" % j
		#print "\t\tf = scale[%i];" % j
		f ^= scale[j]
		for k in range(j+1-1,2-1,-1): 
			L[Im(j,k)] += f * ((0.5*sqrt(float((j+k+1)*(j+k+2))))*(my*Y[Re(j+1,k+1)] + mx*Y[Im(j+1,k+1)]) 
					 + (0.5*sqrt(float((j-k+1)*(j-k+2))))*(my*Y[Re(j+1,k-1)] - mx*Y[Im(j+1,k-1)])
					 - (sqrt(float((j-k+1)*(j+k+1))))*mz*Y[Im(j+1,k)] )
			L[Re(j,k)] += f * ((0.5*sqrt(float((j+k+1)*(j+k+2))))*(mx*Y[Re(j+1,k+1)] - my*Y[Im(j+1,k+1)])
					 - (0.5*sqrt(float((j-k+1)*(j-k+2))))*(mx*Y[Re(j+1,k-1)] + my*Y[Im(j+1,k-1)]) 
					 - (sqrt(float((j-k+1)*(j+k+1))))*mz*Y[Re(j+1,k)] )
		if j>0:
			L[Im(j,1)] += f * ((0.5*sqrt(float((j+1+1)*(j+1+2))))*(my*Y[Re(j+1,1+1)] + mx*Y[Im(j+1,1+1)]) 
					 + (0.5*sqrt(float((j-1+1)*(j-1+2))))*my*Y[Re(j+1,1-1)] 
					 - (sqrt(float((j-1+1)*(j+1+1))))*mz*Y[Im(j+1,1)] )
			L[Re(j,1)] += f * ((0.5*sqrt(float((j+1+1)*(j+1+2))))*(mx*Y[Re(j+1,1+1)] - my*Y[Im(j+1,1+1)]) 
				         - (0.5*sqrt(float((j-1+1)*(j-1+2))))*mx*Y[Re(j+1,1-1)] 
					 - (sqrt(float((j-1+1)*(j+1+1))))*mz*Y[Re(j+1,1)] )
		L[Im(j,0)] += 0
		L[Re(j,0)] += f * ((sqrt(float(j+1)*(j+2)))*(mx*Y[Re(j+1,1)] - my*Y[Im(j+1,1)]) - (j+1)*mz*Y[Re(j+1,0)])
		print "#endif /* FMM_P_MAX >= %i */" % j
	end_block()

	end_block()
	print

from gen_straightline_code import Op, var, parameter, templates, begin_block, end_block, fini_horiz_add
p_max = 30

print '''
/* This file is automatically generated by gen__gen_eval_expansions.py */
/* DO NOT EDIT! */

#include"_fmmv.h"
'''
import sys

if len(sys.argv)>1:
	par = sys.argv[1]
	print '#include "simd.h"'
	print "#if (FMM_PRECISION==0)"
        print '   #include"%ss.h"' % par
        print "#else"
        print '   #include"%sd.h"' % par
        print "#endif"
	suf = "_"+par
        Op.templates = templates[par]
else:
	suf = ""
    	Op.templates = templates["generic"]

gen_core_gen_M_L("core_gen_M_L"+suf, p_max);	
gen_core_eval_L_M("core_eval_L_M"+suf, p_max);	
gen_core_eval_L_M_grad_minus("core_eval_L_M_grad_minus"+suf, p_max);	
gen_core_eval_L_M_grad_plus("core_eval_L_M_grad_plus"+suf,p_max);	
gen_core_gen_M_L_dipole_minus("core_gen_M_L_dipole_minus"+suf, p_max);	
gen_core_gen_M_L_dipole_plus("core_gen_M_L_dipole_plus"+suf, p_max);	

