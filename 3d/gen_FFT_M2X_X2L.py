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

from pyparsing import Literal,Word,Combine,Optional,Group,alphas,nums,alphanums,\
Forward,ZeroOrMore,cStyleComment,ParseResults,CaselessLiteral,CharsNotIn

fun_name = ""
constants = [] 
indices = []
vars = []	
arrays = []
extern_arrays = []

def parse_spiral_output(input_string): 
	global fun_name, indices, vars, arrays, extern_arrays, constants, statements

	fun_name = ""
	constants = [] 
	indices = []
	vars = []	
	arrays = []
	extern_arrays = []

	ident = Word(alphas, alphas+nums+"_")
	fnumber = Combine( Word( "+-"+nums, nums ) +
			   Optional( "." + Optional( Word( nums ) ) ) +
			   Optional( CaselessLiteral("E") + Word( "+-"+nums, nums ) ) )
	def fnumber_hndl(s,l,t): 
		if not t[0] in constants:
			constants.append(t[0])
	fnumber.setParseAction(fnumber_hndl)
			   
	inumber = Word( "+-"+nums, nums ) 		   


	fun_decl = (Literal("void") + ident + Literal("(double *y, double *x)"))

	def fun_decl_hndl(s,l,t): 
		global fun_name
		fun_name = t[1]

	fun_decl.setParseAction(fun_decl_hndl)

	int_decl = Literal("int") + ident + Literal(";")
	def int_decl_hndl(s,l,t): 
		indices.append(t[1])

	int_decl.setParseAction(int_decl_hndl)

	double_decl = Literal("double") + ident + Literal(";") 
	def double_decl_hndl(s,l,t): 
		vars.append(t[1])

	double_decl.setParseAction(double_decl_hndl)
		


	array_decl = (Literal("static") + Literal("double") + 
		ident + Literal("[") + inumber + Literal("]")+Literal(";")) 
	def array_decl_hndl(s,l,t):
		arrays.append((t[2], int(t[4])))

	array_decl.setParseAction(array_decl_hndl)


	extern_array_decl = (Literal("extern") + Literal("double") + 
		ident + Literal("[") + inumber +
		Literal("]")+Literal(";")) 

	def extern_array_decl_hndl(s,l,t):
		extern_arrays.append((t[2], int(t[4])))

	extern_array_decl.setParseAction(extern_array_decl_hndl)

		
	op = Word("+-*/", max=1)
	def op_hndl(s,l,t):
		if t[0]=="+": return ["add"]
		if t[0]=="-": return ["sub"]
		if t[0]=="*": return ["mul"]
		if t[0]=="/": return ["div"]
	op.setParseAction(op_hndl)

	unary_neg = Word("-", max=1)
	def unary_neg_hndl(s,l,t):
		return ["neg"]
	unary_neg.setParseAction(unary_neg_hndl)
		
	arrayitem = Group(ident + Literal("[").suppress() + ZeroOrMore( CharsNotIn("]")) +Literal("]").suppress())
	atom = (arrayitem ^ ident ^ fnumber)
	simple_expr = Group( atom ^ (atom + op + atom) ^ (unary_neg + atom) ^ (unary_neg + atom + op + atom))
	assign_stmt = Group((ident ^ arrayitem) + "=" + simple_expr  + Literal(";").suppress())
	decl_list = ZeroOrMore( int_decl | double_decl | array_decl | extern_array_decl )
	for_stmt = Forward()
	stmt_list = ZeroOrMore( assign_stmt| for_stmt )

	for_stmt << Group(Literal("for") + 
		Group(Literal("(").suppress() +
		ZeroOrMore( CharsNotIn(";"))+
		Literal(";").suppress() + 
		ZeroOrMore( CharsNotIn(";"))+
		Literal(";").suppress() + 
		ZeroOrMore( CharsNotIn(")"))+
		Literal(")").suppress()) + 
		Literal("{").suppress() +
			Group(stmt_list) +
		Literal("}").suppress())	

	fun = fun_decl.suppress() + Literal("{").suppress() + decl_list.suppress() + stmt_list + Literal("}").suppress()
	fun.ignore(cStyleComment)
	
	
	statements = fun.parseString(input_string)
	#return (fun_name, indices, vars, arrays, extern_arrays, constants, statements)

def gen_code(name):
	global indices, vars, arrays, extern_arrays, constants, statements, array_names, extern_array_names
	print "void %s(%s *y, %s *x)" % (name, templates["basetype"], templates["basetype"])
	print "{"
	indent = 1
	tabs = "\t"*indent
	for s in indices:
		print "%sint %s;" % (tabs, s)
	for s in vars:	
		print "%s%s %s;" % (tabs, templates["type"], s)
	i=0
	for s in constants:
		print "%s%s c%i = %s;" % (tabs, templates["type"], i, templates["set1"] % s)
		i += 1
		
	for (s, n) in arrays:	
		print "%sstatic %s %s[%i];" % (tabs, templates["basetype"], s, n*templates["short_vector_length"])
	for (s, n) in extern_arrays:	
		print "%sextern %s %s[%i];" % (tabs, templates["basetype"], s, n)
	
	print	

	array_names = dict(arrays).keys() + ["x", "y"]
	extern_array_names = dict(extern_arrays).keys() 

	def get_atom(x):
		global vars, constants, array_names, extern_array_names
		if x in vars:
			return x
		if x in constants:
			return "c%i" % constants.index(x)
		if x[0] in array_names:
			if templates["short_vector_length"]>1:
				return templates["getitem"] % (x[0], "%i*(%s)"%(templates["short_vector_length"], x[1]))
			else:
				return templates["getitem"] % (x[0], "%s" % x[1])
		if x[0] in extern_array_names:
			return templates["getbaseitem"] % (x[0], "%s" % x[1])

	def get_rhs(x):		
		global vars, constants, arrays, extern_arrays
		if type(x)==ParseResults: 
			if len(x)==3: # A ? B
				return templates[x[1]] % (get_atom(x[0]), get_atom(x[2]))
			elif len(x)==2: # -A
				return templates[x[0]] % (get_atom(x[1]))
			elif len(x)==1: # A
				return get_atom(x[0])
			elif len(x)==4:  # -A ? B
				return templates[x[2]] % (templates[x[0]] % get_atom(x[1]), get_atom(x[3]))
		#else:
		#	return get_atom(x[0])
			
		
	def print_statement_sequence(statements, indent):
		global vars, arrays, extern_arrays, constants
		tabs = "\t"*indent
		for stmnt in statements:
			if stmnt[0] == "for":
				print "%sfor(%s; %s; %s) {" % (tabs, stmnt[1][0], stmnt[1][1], stmnt[1][2])
				print_statement_sequence(stmnt[2], indent+1)
				print "%s}" % tabs
			elif len(stmnt)==3 and stmnt[1]=="=":
				if stmnt[0] in vars:
					print "%s%s = %s;" % (tabs, stmnt[0], get_rhs(stmnt[2]))
				else: 
					if templates["short_vector_length"]>1:
						index = "%i*(%s)" % (templates["short_vector_length"],stmnt[0][1])
					else:
						index = stmnt[0][1]
					print "%s%s;" % (tabs, templates["setitem"] % (stmnt[0][0], index, get_rhs(stmnt[2])))
				
	
	print_statement_sequence(statements, indent)

	print "}"

from sys import stdin, argv

from gen_straightline_code import templates

M_list = [2, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64]
dirname = "SPIRAL_FFT_M2X_X2L"
M2X_template = "FFT_M2X_%i.c"
X2L_template = "FFT_X2L_%i.c"

print '''
/* This file is automatically generated by gen_FFT_M2X_X2L_aux.py */
/* DO NOT EDIT! */

#include"_fmmv.h"
'''


if len(argv)>1 and argv[1]=="cores_init":
	ident = Word(alphas, alphas+nums+"_")
	fun_decl = (Literal("void") + ident + Literal("(double *y, double *x)"))
	group = Forward() 
	group << "{" + ZeroOrMore(group^CharsNotIn("{}")) + "}" 
	tail = fun_decl.suppress() + group.suppress() + ZeroOrMore(CharsNotIn(""))
	tail.ignore(cStyleComment)
	
	for M in M_list:
		file = open(dirname+"/"+(M2X_template %  M))
		t = tail.parseString("".join(file.readlines()))	
		print t[0].replace("double", "_FLOAT_")
		file.close()
		file = open(dirname+"/"+(X2L_template %  M))
		t = tail.parseString("".join(file.readlines()))	
		print t[0].replace("double", "_FLOAT_")
		file.close()

elif len(argv)>1 and "cores" in argv[1]: ####################################################
	if "simd2" in argv[1]:
		par = "simd2"
	        suf = "_simd2"
        	templates = templates[par]
	elif "simd4" in argv[1]:
		par = "simd4"
	        suf = "_simd4"
        	templates = templates[par]
	else:
		suf = ""
		templates = templates["generic"]
	if "simd" in argv[1]:
		print '#include "simd.h"'
		print "#if (FMM_PRECISION==0)"
		print '   #include"%ss.h"' % par
		print "#else"
		print '   #include"%sd.h"' % par
		print "#endif"
	for M in M_list:
		file = open(dirname+"/"+(M2X_template % M))
		parse_spiral_output("".join( file.readlines() )) 
		gen_code(fun_name+suf)	
		file.close()
		file = open(dirname+"/"+(X2L_template % M))
		parse_spiral_output("".join( file.readlines() )) 
		gen_code(fun_name+suf)	
		file.close()
