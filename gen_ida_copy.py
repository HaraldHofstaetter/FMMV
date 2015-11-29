from fmmv_setup import dim, charge_kinds, pot_kinds, charge_combs, pot_combs, charge_pot_combs

import sys

simd_vector_length=1
if len(sys.argv)>1:
	if "simd2" in sys.argv[1]:
		simd_vector_length=2
		div_shift = 1
	elif "simd4" in sys.argv[1]:
		simd_vector_length=4
		div_shift = 2

def gen_loop(tt, access, what, source_target, target=False):
	# what ... "zero", "copy", or "backcopy"
	if target:
		nn = "FMMV->NTargets"
		nv = "NCompleteTargetVectors"
		pp = "permTargets"
	else:	
		nn = "FMMV->NParticles"
		nv = "NCompleteParticleVectors"
		pp = "perm"
	xxx = [None, None, "0.0"]
	if what=="zero":
		i1=0; i2=2
	elif what=="copy":
		i1=0; i2=1
	elif what=="backcopy":
		i1=1; i2=0
        if simd_vector_length==1:
		print "%sfor (i = 0; i < %s; i++) {" % (tt, nn)
		if len(access)==1:
			xxx[0] = "access_%s(i)" % access[0] 
			xxx[1] = "%s[%s[i]]" % (source_target, pp)
			print "%s\t%s = %s;" % (tt, xxx[i1], xxx[i2])
		else:
			jj = 0
			for j in access:
				xxx[0] = "access_%s(i)" % j
				xxx[1] = "%s[%s[i]][%i]" % (source_target, pp, jj)
				print "%s\t\t%s = %s;" % (tt, xxx[i1], xxx[i2])
				jj += 1	
		print "%s}" % tt
        else:
		print "%sfor (i = 0, ii=0; i < %s; i++, ii+=%i) {" % (tt, nv, simd_vector_length) 
		if len(access)==1:
			xxx[0] = "access_%s(i)[0]" % access[0] 
			xxx[1] = "%s[%s[ii]]" % (source_target, pp)
			print "%s\t%s = %s;" % (tt, xxx[i1], xxx[i2])
			xxx[0] = "access_%s(i)[1]" %  access[0]
			xxx[1] = "%s[%s[ii+1]]" % (source_target, pp)
			print "%s\t%s = %s;" % (tt, xxx[i1], xxx[i2])
			if simd_vector_length==4:
				xxx[0] = "access_%s(i)[2]" %  access[0]
				xxx[1] = "%s[%s[ii+2]]" % (source_target, pp)
				print "%s\t%s = %s;" % (tt, xxx[i1], xxx[i2])
				xxx[0] = "access_%s(i)[3]" %  access[0]
				xxx[1] = "%s[%s[ii+3]]" % (source_target, pp)
				print "%s\t%s = %s;" % (tt, xxx[i1], xxx[i2])
		else:
			jj = 0
			for j in access:
				xxx[0] = "access_%s(i)[0]" % j
				xxx[1] = "%s[%s[ii]][%i]" % (source_target, pp, jj)
				print "%s\t%s = %s;" % (tt, xxx[i1], xxx[i2])
				xxx[0] = "access_%s(i)[1]" % j
				xxx[1] = "%s[%s[ii+1]][%i]" % (source_target, pp, jj)
				print "%s\t%s = %s;" % (tt, xxx[i1], xxx[i2])
				if simd_vector_length==4:
					xxx[0] = "access_%s(i)[2]" % j
					xxx[1] = "%s[%s[ii+2]][%i]" % (source_target, pp, jj)
					print "%s\t%s = %s;" % (tt, xxx[i1], xxx[i2])
					xxx[0] = "access_%s(i)[3]" % j
					xxx[1] = "%s[%s[ii+3]][%i]" % (source_target, pp, jj)
					print "%s\t%s = %s;" % (tt, xxx[i1], xxx[i2])
				jj += 1	
		print "%s}" % tt
		if simd_vector_length==2:
			print "%sif (%s&1) {" % (tt, nn)
			if what != "zero":
				print "%s\tii = %s - 1;" % (tt, nn)
			if len(access)==1:
				xxx[0] = "access_%s(%s)[0]" % (access[0], nv)
				xxx[1] = "%s[%s[ii]]" % (source_target, pp)
				print "%s\t%s = %s;" % (tt, xxx[i1], xxx[i2])
			else:				
				jj = 0
				for j in access:
					xxx[0] = "access_%s(%s)[0]" % (j, nv) 
					xxx[1] = "%s[%s[ii]][%i]" % (source_target, pp, jj)
					print "%s\t%s = %s;" % (tt, xxx[i1], xxx[i2])
					jj += 1
			print "%s}" % tt
		elif simd_vector_length==4:
			if what != "zero":
				print "%sii = (%s<<2);" % (tt,nv)
			print "%sfor (i=0; i < (%s&3); i++) {" % (tt, nn)
			if len(access)==1:
				xxx[0] = "access_%s(%s)[i]" % (access[0], nv)
				xxx[1] = "%s[%s[ii+i]]" % (source_target, pp)
				print "%s\t%s = %s;" % (tt, xxx[i1], xxx[i2])
			else:
				jj = 0
				for j in access:
					xxx[0] = "access_%s(%s)[i]" % (j, nv)
					xxx[1] = "%s[%s[ii+i]][%i]" % (source_target, pp, jj)
					print "%s\t%s = %s;" % (tt, xxx[i1], xxx[i2])
					jj += 1
			print "%s}" % tt


	

def gen_copy_particles(suf, st):
	print "static void copy_particles%s(FmmvHandle *FMMV)" % (suf+st)
	print "{"
	print "\tDEFINE_IDA_LOCAL_ALIASES(FMMV)"
	print "\t_FLOAT_ (*particles)[%i] = FMMV->particles;" % dim
	print "\tint *perm = FMMV->perm;"
	if st == "_st":
		print "\t_FLOAT_ (*targets)[%i] = FMMV->targets;" % dim
		print "\tint *permTargets = FMMV->permTargets;"
	print "\tint i;"
        if simd_vector_length>1:
		print "\tint ii;"
		print "\tint NCompleteParticleVectors = FMMV->NParticles>>%i;" % div_shift
		if st == "_st":
			print "\tint NCompleteTargetVectors = FMMV->NTargets>>%i;" % div_shift
        access = ["x", "y", "z"][0:dim]
	gen_loop("\t", access, "copy", "particles", target=False)
	if st == "_st":
                access = ["tx", "ty", "tz"][0:dim]
		gen_loop("\t", access, "copy", "targets", target=True)
	print "}"
	print




def gen_copy_charges(suf, st):
	print "static void copy_charges%s(FmmvHandle *FMMV)" % (suf+st)
	print "{"
	print "\tDEFINE_IDA_LOCAL_ALIASES(FMMV)"
	print "\tint *perm = FMMV->perm;"
	charge_comb = charge_combs[charge_pot_combs[suf][0]]
	for i in charge_comb:
		d = len(charge_kinds[i][1])
		n = charge_kinds[i][0]
		if d==1:
			print "\t_FLOAT_ *%s = FMMV->%s;" % (n, n)
		else:
			print "\t_FLOAT_ (*%s)[%i] =  FMMV->%s;" % (n, d, n)
	print "\tint i;"
        if simd_vector_length>1:
		print "\tint ii;"
		print "\tint NCompleteParticleVectors = FMMV->NParticles>>%i;" % div_shift
	for i in charge_comb:
		n = charge_kinds[i][0]
		v = charge_kinds[i][1]
		print
		print "\tif (%s==0) {" % n
		gen_loop("\t\t", v, "zero", n, target=False)
		print "\t}"
		print "\telse {"
		gen_loop("\t\t", v, "copy", n, target=False)
		print "\t}"	
	print "}"
	print


def gen_backcopy_pot(suf, st):
	print "static void backcopy_pot%s(FmmvHandle *FMMV)" % (suf+st)
	print "{"
	print "\tDEFINE_IDA_LOCAL_ALIASES(FMMV)"
	print "\tint *permTargets = FMMV->permTargets;"
	pot_comb = pot_combs[charge_pot_combs[suf][1]]
	for i in pot_comb:
		d = len(pot_kinds[i][1])
		n = pot_kinds[i][0]
		if d==1:
			print "\t_FLOAT_ *%s = FMMV->%s;" % (n, n)
		else:
			print "\t_FLOAT_ (*%s)[%i] =  FMMV->%s;" % (n, d, n)
	print "\tint i;"
        if simd_vector_length>1:
		print "\tint ii;"
		print "\tint NCompleteTargetVectors = FMMV->NTargets>>%i;" % div_shift
	print

	for i in pot_comb:
		n = pot_kinds[i][0]
		v = pot_kinds[i][1]
		print
		print "\tif (%s!=0) {" % n
		gen_loop("\t\t", v, "backcopy", n, target=True)
		print "\t}"
	print "}"
	print

def gen_zero_pot(suf, st):
	print "static void zero_pot%s(FmmvHandle *FMMV)" % (suf+st)
	print "{"
	print "\tDEFINE_IDA_LOCAL_ALIASES(FMMV)"
	#print "\tint *permTargets = FMMV->permTargets;"
	pot_comb = pot_combs[charge_pot_combs[suf][1]]
	print "\tint i;"
        if simd_vector_length>1:
		print "\tint ii;"
		print "\tint NCompleteTargetVectors = FMMV->NTargets>>%i;" % div_shift
	for i in pot_comb:
		v = pot_kinds[i][1]
		print
		gen_loop("\t", v, "zero", None, target=True)
	print "}"
	print

def gen_dispatch(rout_name, charge_pot_combs):
	print "void %s(FmmvHandle *FMMV)" % rout_name
	print "{"
        print "\tswitch(FMMV->dataKind) {"
	for suf in charge_pot_combs: 	
	       	print "\tcase FMM%s:" % suf.upper()
                print "\t\t%s(FMMV);" % (rout_name + suf)
		print "\t\tbreak;"
	for suf in charge_pot_combs: 	
	       	print "\tcase FMM_ST%s:" % suf.upper()
	        print "\t\t%s_st(FMMV);" % (rout_name + suf )
		print "\t\tbreak;"
	print "\t}"		
	print "}"
	print



print "/* This file is automatically generated by gen_gen_neighbor2_list.py */"
print "/* DO NOT EDIT! */"
print
print '#include"_fmmv.h"'
print
for suf in charge_pot_combs: 	
	print "#undef FMM_KIND"
	print "#define FMM_KIND FMM%s" % suf.upper()
	print '#include"fmmv_access.h"'
	print
	gen_copy_particles(suf, "")
	gen_copy_charges(suf, "")
	gen_backcopy_pot(suf, "")
	gen_zero_pot(suf, "")
	print "#undef FMM_KIND"
	print "#define FMM_KIND FMM_ST%s" % suf.upper()
	print '#include"fmmv_access.h"'
	print
	gen_copy_particles(suf, "_st")
	gen_copy_charges(suf, "_st")
	gen_backcopy_pot(suf, "_st")
	gen_zero_pot(suf, "_st")
	
	
for rout_name in ["copy_particles", "copy_charges", "backcopy_pot", "zero_pot"]:
	gen_dispatch(rout_name, charge_pot_combs)




