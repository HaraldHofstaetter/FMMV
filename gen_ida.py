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

def gen_struct(name, suffix, a, number_type, simd_vector_length=1):
	vc = 0
	for i in a:
		if type(i)==list:
			print "\ntypedef struct {"
			for j in i:
				if simd_vector_length>1:
					print "\t%s %s[%i];" % (number_type, j, simd_vector_length)
				else:
					print "\t%s %s;" % (number_type, j)
			print "} _%s%s_V%i_t;" % (name, suffix, vc)
			vc += 1
		
	print "\ntypedef struct {"
	vc = 0
	for i in a:
		if type(i)==list:
			print  ("\t_%s%s_V%i_t *_%i; /* %s" % (name,suffix, vc, vc, i[0])) ,
			for j in i[1:]:
				print  ", %s" % j ,
			print  " */"
			vc += 1
		else:	
			if simd_vector_length>1:
				print  "\t%s (*%s)[%i];" % (number_type, i, simd_vector_length)
			else:
				print  "\t%s *%s;" % (number_type, i)
	print  "} %s%s_t;" % (name, suffix)	

def gen_access(name, suffix, a, alias):
	vc = 0
	r = {}
	for i in a:
		if type(i)==list:
			for j in i:
				r[j] = "_%i[%%s].%s" % (vc, j)
			vc += 1	
		else:
			r[i] = "%s[%%s]" % i
			
	for i in a:
		if type(i)==list:
			for j in i:
				print "#define access_%s(i) (((%s%s_t*)FMMV->DATA)->%s)" % (j, name, suffix, r[j] % "i")
		else:
			print "#define access_%s(i) (((%s%s_t*)FMMV->DATA)->%s)" % (i, name, suffix, r[i] % "i")
	for (i,j) in alias:
		print "#define access_%s(i) (((%s%s_t*)FMMV->DATA)->%s)" % (i, name, suffix, r[j] % "i")
		
		
def gen_define_local_aliases(name, suffix, a, short_vectort_length=1):
	print "#define DEFINE_IDA_LOCAL_ALIASES(FMMV) \\"
	vc = 0
	for i in a:
		if type(i)==list:
			print "\t_%s%s_V%i_t *_%i = ((%s%s_t*)FMMV->DATA)->_%i; \\" % (name, suffix, vc, vc, name, suffix, vc)
			vc += 1
		else:
			if simd_vector_length>1:
				print "\t_FLOAT_ (*_%s)[%i] = ((%s%s_t*)FMMV->DATA)->%s; \\" % (i, simd_vector_length, name, suffix, i)
			else:
				print "\t_FLOAT_ *_%s = ((%s%s_t*)FMMV->DATA)->%s; \\" % (i, name, suffix, i)
	print		
	print

def gen_access_local_aliases(name, suffix, a, alias):
	vc = 0
	r = {}
	for i in a:
		if type(i)==list:
			for j in i:
				r[j] = "_%i[%%s].%s" % (vc, j)
			vc += 1	
		else:
			r[i] = "_%s[%%s]" % i
			
	for i in a:
		if type(i)==list:
			for j in i:
				print "#define access_%s(i) (%s)" % (j, r[j] % "i")
		else:
			print "#define access_%s(i) (%s)" % (i, r[i] % "i")
	for (i,j) in alias:
		print "#define access_%s(i) (%s)" % (i, r[j] % "i")
		
			

def gen_allocate(name, suffix, a, targets, number_type, simd_vector_length=1):
	print "static int ida_allocate%s(FmmvHandle *FMMV)" % suffix
	print "{"
	if simd_vector_length==2:
		print "\tint NPairs = FMMV->NParticles/2 + 1;"
		if len(targets)>0:
			print "\tint NTargetPairs = FMMV->NTargets/2 + 1;"
	elif simd_vector_length==4:
		print "\tint NQuadruples = FMMV->NParticles/4 + 1;"
		if len(targets)>0:
			print "\tint NTargetQuadruples = FMMV->NTargets/4 + 1;"
	print			
	print "\tFMMV->DATA = (void*) FMMV_MALLOC(FMMV, sizeof(%s%s_t));" % (name, suffix)		
	print "\tif (!FMMV->DATA) goto _err;" 
	vc = 0
	for i in a:
		if type(i)==list:
			if len(i)!=0 and not i[0] in targets:
				if simd_vector_length==2:
					print ("\t((%s%s_t*)FMMV->DATA)->_%i = (_%s%s_V%i_t *) FMMV_MALLOC(FMMV, NPairs*sizeof(_%s%s_V%i_t)); " %
						(name, suffix, vc, name, suffix, vc, name, suffix,vc))
				elif simd_vector_length==4:
					print ("\t((%s%s_t*)FMMV->DATA)->_%i = (_%s%s_V%i_t *) FMMV_MALLOC(FMMV, NQuadruples*sizeof(_%s%s_V%i_t)); " %
						(name, suffix, vc, name, suffix, vc, name, suffix,vc))
				else:		
					print ("\t((%s%s_t*)FMMV->DATA)->_%i = (_%s%s_V%i_t *) FMMV_MALLOC(FMMV, FMMV->NParticles*sizeof(_%s%s_V%i_t)); " %
						(name, suffix, vc, name, suffix, vc, name, suffix,vc))
			else:		
				if simd_vector_length==2:
					print ("\t((%s%s_t*)FMMV->DATA)->_%i = (_%s%s_V%i_t *) FMMV_MALLOC(FMMV, NTargetPairs*sizeof(_%s%s_V%i_t)); " %
						(name, suffix, vc, name, suffix, vc, name, suffix,vc))
				elif simd_vector_length==4:
					print ("\t((%s%s_t*)FMMV->DATA)->_%i = (_%s%s_V%i_t *) FMMV_MALLOC(FMMV, NTargetQuadruples*sizeof(_%s%s_V%i_t)); " %
						(name, suffix, vc, name, suffix, vc, name, suffix,vc))
				else:
					print ("\t((%s%s_t*)FMMV->DATA)->_%i = (_%s%s_V%i_t *) FMMV_MALLOC(FMMV, FMMV->NTargets*sizeof(_%s%s_V%i_t)); " %
						(name, suffix, vc, name, suffix, vc, name, suffix,vc))
			print "\tif (!((%s%s_t*)FMMV->DATA)->_%i) goto _err;" % (name, suffix, vc)
			vc +=1 
		else:
			if not i in targets:
				if simd_vector_length==2:
					print ("\t((%s%s_t*)FMMV->DATA)->%s = (%s (*)[2]) FMMV_MALLOC(FMMV, 2*NPairs*sizeof(%s)); " %
						(name, suffix, i, number_type, number_type))
				elif simd_vector_length==4:
					print ("\t((%s%s_t*)FMMV->DATA)->%s = (%s (*)[4]) FMMV_MALLOC(FMMV, 4*NQuadruples*sizeof(%s)); " %
						(name, suffix, i, number_type, number_type))
				else:
					print ("\t((%s%s_t*)FMMV->DATA)->%s = (%s *) FMMV_MALLOC(FMMV, FMMV->NParticles*sizeof(%s)); " %
						(name, suffix, i, number_type, number_type))
			else:		
				if simd_vector_length==2:
					print ("\t((%s%s_t*)FMMV->DATA)->%s = (%s (*)[2]) FMMV_MALLOC(FMMV, 2*NTargetPairs*sizeof(%s)); " %
						(name, suffix, i, number_type, number_type))
				elif simd_vector_length==4:
					print ("\t((%s%s_t*)FMMV->DATA)->%s = (%s (*)[4]) FMMV_MALLOC(FMMV, 4*NTargetQuadruples*sizeof(%s)); " %
						(name, suffix, i, number_type, number_type))
				else:
					print ("\t((%s%s_t*)FMMV->DATA)->%s = (%s *) FMMV_MALLOC(FMMV, FMMV->NTargets*sizeof(%s)); " %
						(name, suffix, i, number_type, number_type))
			print "\tif (!((%s%s_t*)FMMV->DATA)->%s) goto _err;" % (name, suffix, i)
	print 
	print "\treturn 0;"
	print 
	print "_err:"
	print "\treturn 1;"
	print "}"
	print 
	


def gen_free(name, suffix, a, targets, simd_vector_length=1):
	print "static void ida_free%s(FmmvHandle *FMMV)" % suffix
	print "{"
	if simd_vector_length==2:
		print "\tint NPairs = FMMV->NParticles/2 + 1;"
		if len(targets)>0:
			print "\tint NTargetPairs = FMMV->NTargets/2 + 1;"
	elif simd_vector_length==4:
		print "\tint NQuadruples = FMMV->NParticles/4 + 1;"
		if len(targets)>0:
			print "\tint NTargetQuadruples = FMMV->NTargets/4 + 1;"
	vc = 0


	for i in a:
		if type(i)==list:
			if len(i)!=0 and not i[0] in targets:
				if simd_vector_length==2:
					print ("\tFMMV_FREE(FMMV, ((%s%s_t*)FMMV->DATA)->_%i, NPairs*sizeof(_%s%s_V%i_t)); " %
						(name, suffix, vc, name, suffix,vc))
				elif simd_vector_length==4:
					print ("\tFMMV_FREE(FMMV, ((%s%s_t*)FMMV->DATA)->_%i, NQuadruples*sizeof(_%s%s_V%i_t)); " %
						(name, suffix, vc, name, suffix,vc))
				else:		
					print ("\tFMMV_FREE(FMMV, ((%s%s_t*)FMMV->DATA)->_%i, FMMV->NParticles*sizeof(_%s%s_V%i_t)); " %
						(name, suffix, vc, name, suffix,vc))
			else:		
				if simd_vector_length==2:
					print ("\tFMMV_FREE(FMMV, ((%s%s_t*)FMMV->DATA)->_%i, NTargetPairs*sizeof(_%s%s_V%i_t)); " %
						(name, suffix, vc, name, suffix,vc))
				elif simd_vector_length==4:
					print ("\tFMMV_FREE(FMMV, ((%s%s_t*)FMMV->DATA)->_%i, NTargetQuadruples*sizeof(_%s%s_V%i_t)); " %
						(name, suffix, vc, name, suffix,vc))
				else:
					print ("\tFMMV_FREE(FMMV, ((%s%s_t*)FMMV->DATA)->_%i, FMMV->NTargets*sizeof(_%s%s_V%i_t)); " %
						(name, suffix, vc, name, suffix,vc))
			vc +=1 
		else:
			if not i in targets:
				if simd_vector_length==2:
					print ("\tFMMV_FREE(FMMV, ((%s%s_t*)FMMV->DATA)->%s, 2*NPairs*sizeof(%s)); " % (name, suffix, i, number_type)) 
				elif simd_vector_length==4:
					print ("\tFMMV_FREE(FMMV, ((%s%s_t*)FMMV->DATA)->%s, 4*NQuadruples*sizeof(%s)); " % (name, suffix, i, number_type))
				else:
					print ("\tFMMV_FREE(FMMV, ((%s%s_t*)FMMV->DATA)->%s, FMMV->NParticles*sizeof(%s)); " % (name, suffix, i, number_type))
			else:		
				if simd_vector_length==2:
					print ("\tFMMV_FREE(FMMV, ((%s%s_t*)FMMV->DATA)->%s, 2*NTargetPairs*sizeof(%s)); " % (name, suffix, i, number_type))
				elif simd_vector_length==4:
					print ("\tFMMV_FREE(FMMV, ((%s%s_t*)FMMV->DATA)->%s, 4*NTargetQuadruples*sizeof(%s)); " % (name, suffix, i, number_type))
				else:
					print ("\tFMMV_FREE(FMMV, ((%s%s_t*)FMMV->DATA)->%s, FMMV->NTargets*sizeof(%s)); " % (name, suffix, i, number_type))
	print "\tFMMV_FREE(FMMV, FMMV->DATA, sizeof(%s%s_t));" % (name, suffix)		
	print "}"
	print 
	
###############################################################################################




import sys		
from fmmv_setup import simd_vector_length, ida, ida_st, charge_pot_combs, pot_kinds, float_types, dim, name, dir, charge_kinds, additional_data_in_box, additional_global_data, type_name, p_max, s_eps_max, s_exp_max, eval_direct_number_of_accuracies

targets = ['tx', 'ty', 'tz']
for i in pot_kinds:
	targets = targets + pot_kinds[i][1]

print "/* This file is automatically generated by gen_ida.py */"
print "/* DO NOT EDIT! */"
print

if "header" in sys.argv[1]:
	print '#ifndef _FMMV_H_'
	print '#define _FMMV_H_'
	print
	print '#define FMM_DIM  %i' % dim
	print
	i=0
	for float_type in float_types:
		if i==0:
			print "#if (FMM_PRECISION==0)"
		else:	
			print "#elif (FMM_PRECISION==%i)" % i
		i += 1	
		print '\ttypedef %s _FLOAT_;' % type_name[float_type]
		print "\t#define FMM_P_MAX %i" % p_max[float_type]
		print "\t#define FMM_S_EPS_MAX %i" % s_eps_max[float_type]
		print "\t#define FMM_S_EXP_MAX %i" % s_exp_max[float_type]
		print "\t#define FMM_SIMD_VECTOR_LENGTH %i" % simd_vector_length[float_type]
	print "#endif"	
	print 
elif "access" in sys.argv[1]:
	print "#undef DEFINE_IDA_LOCAL_ALIASES"
	for j in ["x","y","z"][0:dim]:
		print "#undef access_%s" % j
	for j in ["tx","ty","tz"][0:dim]:
		print "#undef access_%s" % j
	for i in charge_kinds:
		for j in charge_kinds[i][1]:
			print "#undef access_%s" % j
	for i in pot_kinds:
		for j in pot_kinds[i][1]:
			print "#undef access_%s" % j
	print		
elif "source" in sys.argv[1]:
	print "#undef USE_LOCAL_ALIASES"
	print '#include"_fmmv.h"'

kind_nr = 0
for ida1 in [ida, ida_st]:
	if ida1==ida_st:
		alias = []
		suf1 = "_ST"
	else:
		alias =  [("tx", "x"), ("ty", "y"), ("tz", "z")][0:dim]
		suf1 = ""
   	for suffix in charge_pot_combs:
		if "header" in sys.argv[1]:
			print "#define FMM%s %i" % ((suf1+suffix).upper(), kind_nr)
			kind_nr = kind_nr + 1
			i=0
			for float_type in float_types:
				if i==0:
					print "\n#if (FMM_PRECISION==0)"
				else:	
					print "\n#elif (FMM_PRECISION==%i)" % i
				i += 1	
				gen_struct("DATA", suffix+suf1, ida1[float_type][suffix], float_type, simd_vector_length[float_type])
			print "#endif"
			print
		elif "access" in sys.argv[1]:
                        print "#if (FMM_KIND == FMM%s)" % ((suf1+suffix).upper())
			i=0
			for float_type in float_types:
				if i==0:
					print "\n#if (FMM_PRECISION==0)"
				else:	
					print "\n#elif (FMM_PRECISION==%i)" % i
				i += 1	
				gen_define_local_aliases("DATA", suffix+suf1, ida1[float_type][suffix], simd_vector_length[float_type])
				gen_access_local_aliases("DATA", suffix+suf1, ida1[float_type][suffix], alias)
			print "\n#endif /* FMM_PRECISION */"
                	print "\n#endif /* FMM_KIND == FMM%s */" % ((suf1+suffix).upper())
		elif "source" in sys.argv[1]:	
			i=0
			for float_type in float_types:
				if i==0:
					print "\n#if (FMM_PRECISION==0)"
				else:	
					print "\n#elif (FMM_PRECISION==%i)" % i
				i += 1	
				if "ST" in suf1:
					gen_allocate("DATA", suffix+suf1, ida1[float_type][suffix], targets, float_type, simd_vector_length[float_type])
					gen_free("DATA", suffix+suf1, ida1[float_type][suffix], targets, simd_vector_length[float_type])
				else:
					gen_allocate("DATA", suffix+suf1, ida1[float_type][suffix], [], float_type, simd_vector_length[float_type])
					gen_free("DATA", suffix+suf1, ida1[float_type][suffix], [], simd_vector_length[float_type])	
			print "\n#endif /* FMM_PRECISION */"
		print
		print


if "header" in sys.argv[1]:
	print "#define FMM_ADDITIONAL_DATA_IN_BOX\\"
	for data in additional_data_in_box:
		print "\t%s;\t\\" % data
	print 	
	print "#define FMM_ADDITIONAL_GLOBAL_DATA\\"
	for charge_kind in charge_kinds:
		d = len(charge_kinds[charge_kind][1])
		n = charge_kinds[charge_kind][0]
		if d==1:
			print "\t_FLOAT_ *%s;\\" % n
		else:	
			print "\t_FLOAT_ (*%s)[%i];\\" % (n, d)
	for pot_kind in pot_kinds:
		d = len(pot_kinds[pot_kind][1])
		n = pot_kinds[pot_kind][0]
		if d==1:
			print "\t_FLOAT_ *%s;\\" % n
		else:	
			print "\t_FLOAT_ (*%s)[%i];\\" % (n, d)
	for data in additional_global_data:
		print "\t%s;\t\\" % data
	print	

	print '#include"fmmv%s.h"' % name["generic"]
	print '#include"fmmv_common.h"'
	i=0
	for float_type in float_types:
		if i==0:
			print "#if (FMM_PRECISION==0)"
		else:	
			print "#elif (FMM_PRECISION==%i)" % i
		i += 1	
		if simd_vector_length[float_type] > 1:
			print '  #include"simd.h"'
	print "#endif /* FMM_PRECISION */"
	print
	# gen prototypes
	protos = [
		"void gen_M%s(FmmvHandle *FMMV, Box *box);",
		"void gen_M_ST%s(FmmvHandle *FMMV, Box *box);",
		"void eval_L%s(FmmvHandle *FMMV, Box *box);",
		"void eval_L_ST%s(FmmvHandle *FMMV, Box *box);",
		"void eval_M%s(FmmvHandle *FMMV, Box *target, Box *source);",
		"void eval_M_ST%s(FmmvHandle *FMMV, Box *target, Box *source);",
		"void gen_L%s(FmmvHandle *FMMV, Box *target, Box *source);",
		"void gen_L_ST%s(FmmvHandle *FMMV, Box *target, Box *source);",
		"void gen_L_eval_M%s(FmmvHandle *FMMV, Box *list3, Box *list4);",
		"void extrinsic_correction%s(FmmvHandle *FMMV);",
		"void extrinsic_correction_ST%s(FmmvHandle *FMMV);",
	]
	    
	for proto in protos: 
		for i in charge_pot_combs:
			print proto % i
			
	protos = [
		"void eval_direct%s_acc%i(FmmvHandle *FMMV, Box *target, Box *source);",
		"void eval_direct_ST%s_acc%i(FmmvHandle *FMMV, Box *target, Box *source);",
	]
	if dim==2:
		protos += [
		"void eval_direct_periodic%s_acc%i(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy);",
		"void eval_direct_periodic_ST%s_acc%i(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy);",
		]
	else:	
		protos += [
		"void eval_direct_periodic%s_acc%i(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);",
		"void eval_direct_periodic_ST%s_acc%i(FmmvHandle *FMMV, Box *target, Box *source, _FLOAT_ dx, _FLOAT_ dy, _FLOAT_ dz);",
		]
	for proto in protos: 
		for i in charge_pot_combs:
			for j in range(eval_direct_number_of_accuracies["generic"]):
				print proto % (i, j)
		
	
	print '#endif /* _FMMV_H_ */' 
elif "source" in sys.argv[1]:
	print 
	print "int ida_allocate(FmmvHandle *FMMV)"
	print "{"
	print "\tswitch(FMMV->dataKind) {"
	for st in ["", "_ST"]:
   	    for suffix in charge_pot_combs:	        
	        print "\tcase FMM%s:" % (st+suffix).upper()
		print "\t\treturn ida_allocate%s(FMMV);" %(suffix+st)
	print "\t}"
	print "\treturn -1;"
	print "}"
	print 
	print "void ida_free(FmmvHandle *FMMV)"
	print "{"
	print "\tswitch(FMMV->dataKind) {"
	for st in ["", "_ST"]:
   	    for suffix in charge_pot_combs:	        
	        print "\tcase FMM%s:" % (st+suffix).upper()
		print "\t\tida_free%s(FMMV);" %(suffix+st)
		print "\t\tbreak;"
	print "\t}"
	print "}"
	
        
