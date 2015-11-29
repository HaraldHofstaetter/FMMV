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

templates = {}

templates["single"] = {
	"add": "(%s+%s)",
	"sub": "(%s-%s)",
	"mul": "(%s*%s)",
	"div": "(%s/%s)",
	"neg": "(-%s)",
	"set": "%s",
	"set1": "%s",
	"getitem": "%s[%s]",
	"setitem": "%s[%s] = %s",
	"type": "float",
	"sqrt": "sqrtf(%s)",
	"zero": "0",
	"basetype": "float",
	"getbaseitem": "%s[%s]",
	"short_vector_length": 1,
	"init": "",
}	

templates["double"] = {
	"add": "(%s+%s)",
	"sub": "(%s-%s)",
	"mul": "(%s*%s)",
	"div": "(%s/%s)",
	"neg": "(-%s)",
	"set": "%s",
	"set1": "%s",
	"getitem": "%s[%s]",
	"setitem": "%s[%s] = %s",
	"type": "double",
	"sqrt": "sqrt(%s)",
	"zero": "0",
	"basetype": "double",
	"getbaseitem": "%s[%s]",
	"short_vector_length": 1,
	"init": ""
}	

templates["generic"] = {
        "add": "(%s+%s)",
        "sub": "(%s-%s)",
        "mul": "(%s*%s)",
        "div": "(%s/%s)",
        "neg": "(-%s)",
	"fma": "(%s*%s+%s)",
	"fms": "(%s*%s-%s)",
	"fnms": "(-(%s*%s-%s))",
        "set": "%s",
        "set1": "%s",
        "getitem": "%s[%s]",
        "setitem": "%s[%s] = %s",
        "type": "_FLOAT_",
        "sqrt": "SQRT(%s)",
	"zero": "0",
	"basetype": "_FLOAT_",
	"getbaseitem": "%s[%s]",
	"short_vector_length": 1,
        "init": ""
}

templates["sse1"] = {
	"add": "_mm_add_ps(%s, %s)",
	"sub": "_mm_sub_ps(%s, %s)",
	"mul": "_mm_mul_ps(%s, %s)",
	"div": "_mm_div_ps(%s, %s)",
	"xor": "_mm_xor_ps(%s, %s)",
	#"neg": "_mm_xor_ps(%s, _mm_set1_ps(-0.0))",
	"neg": "_mm_xor_ps(%s, __negative_zero__)",
	"set": "_mm_set_ps(%s, %s, %s, %s)",
	"set1": "_mm_set1_ps(%s)",
	"getitem": "_mm_load_ps(%s+%s)",
	"setitem": "_mm_store_ps(%s+%s, %s)",
	"type": "__m128",
	"sqrt": "_mm_sqrt_ps(%s)",
	"zero": "_mm_setzero_ps()",
	"basetype": "float",
	"getbaseitem": "_mm_load_ps1(%s+%s)",
	"short_vector_length": 4,
	"init": " \t__m128 __negative_zero__ = _mm_set1_ps(-0.0);",
	"horiz_add_init": (
		"\t__m128 _HA_aux_0, _HA_aux_1;\n"
		"#ifndef HORIZ_ADD\n"
		"#define HORIZ_ADD(v0, v1, v2, v3) \\\n"
		"\t(_HA_aux_0 = _mm_add_ps(_mm_movelh_ps(v0, v1), _mm_movehl_ps(v1, v0)), \\\n"
		"\t_HA_aux_1 = _mm_add_ps(_mm_movelh_ps(v2, v3), _mm_movehl_ps(v3, v2)), \\\n"
		"\t_mm_add_ps(_mm_shuffle_ps(_HA_aux_0, _HA_aux_1, 0xDD), _mm_shuffle_ps(_HA_aux_0, _HA_aux_1, 0x88)))\n"
		"#endif"),
}	

templates["sse2"] = {
	"add": "_mm_add_pd(%s, %s)",
	"sub": "_mm_sub_pd(%s, %s)",
	"mul": "_mm_mul_pd(%s, %s)",
	"div": "_mm_div_pd(%s, %s)",
	"xor": "_mm_xor_pd(%s, %s)",
	#"neg": "_mm_xor_pd(%s, _mm_set1_pd(-0.0))",
	"neg": "_mm_xor_pd(%s, __negative_zero__)",
	"set": "_mm_set_pd(%s, %s)",
	"set1": "_mm_set1_pd(%s)",
	"getitem": "_mm_load_pd(%s+%s)",
	"setitem": "_mm_store_pd(%s+%s, %s)",
	"type": "__m128d",
	"sqrt": "_mm_sqrt_pd(%s)",
	"zero": "_mm_setzero_pd()",
	"basetype": "double",
	"short_vector_length": 2,
	"getbaseitem": "_mm_load1_pd(%s+%s)",
	"init": "\t__m128d __negative_zero__ = _mm_set1_pd(-0.0);",
	"horiz_add_init": ("#ifndef HORIZ_ADD\n"
		"#define HORIZ_ADD(v0, v1) _mm_add_pd(_mm_shuffle_pd(v0, v1, _MM_SHUFFLE2(1, 1)), _mm_shuffle_pd(v0, v1, _MM_SHUFFLE2(0, 0)))\n"
		"#endif"),
}	


templates["simd4"] = {
	"add": "V4_ADD(%s, %s)",
	"sub": "V4_SUB(%s, %s)",
	"mul": "V4_MUL(%s, %s)",
	"div": "V4_DIV(%s, %s)",
	"xor": "V4_XOR(%s, %s)",
	"neg": "V4_NEG(%s)",
	"fma": "V4_FMA(%s, %s, %s)",
	"fms": "V4_FMS(%s, %s, %s)",
	"fnms": "V4_FNMS(%s, %s, %s)",
	"set": "V4_SET(%s, %s, %s, %s)",
	"set1": "V4_SET1(%s)",
	"getitem": "V4_GETITEM(%s, %s)",
	"setitem": "V4_SETITEM(%s, %s, %s)",
	"type": "V4_TYPE",
	"sqrt": "V4_SQRT(%s)",
	"zero": "V4_ZERO",
	"basetype": "V4_BASETYPE",
	"getbaseitem": "V4_GETBASEITEM(%s, %s)",
	"short_vector_length": 4,
	"init": "V4_INIT",
	"horizadd": "V4_HORIZADD(%s, %s, %s, %s)",
}	

templates["simd2"] = {
	"add": "V2_ADD(%s, %s)",
	"sub": "V2_SUB(%s, %s)",
	"mul": "V2_MUL(%s, %s)",
	"div": "V2_DIV(%s, %s)",
	"xor": "V2_XOR(%s, %s)",
	"neg": "V2_NEG(%s)",
	"fma": "V2_FMA(%s, %s, %s)",
	"fms": "V2_FMS(%s, %s, %s)",
	"fnms": "V2_FNMS(%s, %s, %s)",
	"set": "V2_SET(%s, %s)",
	"set1": "V2_SET1(%s)",
	"getitem": "V2_GETITEM(%s, %s)",
	"setitem": "V2_SETITEM(%s, %s, %s)",
	"type": "V2_TYPE",
	"sqrt": "V2_SQRT(%s)",
	"zero": "V2_ZERO",
	"basetype": "V2_BASETYPE",
	"getbaseitem": "V2_GETBASEITEM(%s, %s)",
	"short_vector_length": 2,
	"init": "V2_INIT",
	"horizadd": "V2_HORIZADD(%s, %s)",
}	

from math import ceil

indent = 0

def begin_block(s):
	global indent
	print "%s%s{" % (indent * '\t', s)
	if indent==0:
		Op.ops = 0
		print "\t"+Op.templates["init"]
	indent += 1

def end_block():
	global indent
	indent -= 1
	if indent==0:
		print "/* #ops = %i */" % Op.ops
	print "%s}" % (indent * '\t')

class Op:
	def __init__(self, s, array_flag=False, base_array_flag=False, horiz_add_flag=0):
		self.array_flag = array_flag
		self.base_array_flag = base_array_flag
		self.s = s
		self.horiz_add_flag = horiz_add_flag
		if horiz_add_flag!=0:
			if horiz_add_flag < 0:
				self.horiz_add_index = -1 
			else:
				self.horiz_add_index = 0
			self.horiz_add_modulo_index = 0
			h = ()
			for i in range(Op.templates["short_vector_length"]):
				h += ("_HA_%s_%i" % (str(s), i),)
			self.horiz_add_template = Op.templates["setitem"] % (s, "%s", 
				Op.templates["add"] % (Op.templates["getitem"] % (s, "%s"), 
						       Op.templates["horizadd"] % h))
		self.fused_mul = None						       


	def __repr__(self):
		return self.s

	def binary(self, other, which, ops=1, reverse=False):	
		if not isinstance(other, Op): 
			if other==0:
				other = Op(Op.templates["zero"])
			else:	
				other = Op(Op.templates["set1"] % repr(other))
		Op.ops += ops
		if reverse:
			(self, other) = (other, self)
		if not Op.templates.has_key("USE_FMA"):
			return Op(Op.templates[which] % (repr(self), repr(other)))
		
		if which=="add" and self.fused_mul!=None:
			return  Op(Op.templates["fma"] % (self.fused_mul[0], self.fused_mul[1], repr(other)))
		elif which=="add" and other.fused_mul!=None:
			return  Op(Op.templates["fma"] % (other.fused_mul[0], other.fused_mul[1], repr(self)))
		elif which=="sub" and other.fused_mul!=None:
			return  Op(Op.templates["fnms"] % (other.fused_mul[0], other.fused_mul[1], repr(self)))
		elif which=="sub" and self.fused_mul!=None:
			return  Op(Op.templates["fms"] % (self.fused_mul[0], self.fused_mul[1], repr(other)))
		else:
			ret = Op(Op.templates[which] % (repr(self), repr(other)))
			if which=="mul":
				ret.fused_mul = (repr(self), repr(other))
			return ret	

	def unary(self, which, ops=1):
		Op.ops += ops
		return Op(Op.templates[which] % repr(self))

	#Unfortunately, Python doesn't allow to overload the assignment operator '='.
	#We use '^=' instead of '='
	def __ixor__(self, other):
		if not isinstance(other, Op): 
			if other==0:
				other = Op(Op.templates["zero"])
			else:
				other = Op(Op.templates["set1"] % repr(other))
		if self.array_flag:	
			print indent*'\t' + (Op.templates["setitem"] % (repr(self.array), self.array_index, repr(other))) + ';'
		else:	
			print "%s%s = %s;" % (indent*'\t', repr(self), repr(other))
		return self

	def ibinary(self, other, which, ops=1):
		if not isinstance(other, Op): 
			if other==0:
				other = Op(Op.templates["zero"])
			else:	
				other = Op(Op.templates["set1"] % repr(other))
		print "%s%s = %s;" % (indent*'\t', repr(self), Op.templates[which] % (repr(self), repr(other)))
		Op.ops += ops
		return self
		
	def __add__(self, other):
		if other==0.0: return self	
		return Op.binary(self, other, "add")

	def __sub__(self, other):
		if other==0.0: return self	
		return Op.binary(self, other, "sub")
		
	def __mul__(self, other):
		if other==0.0: return 0
		if other==1.0: return self
		return Op.binary(self, other, "mul")
		
	def __div__(self, other):
		#if other==0.0: ***ERROR***
		if other==1.0: return self
		return Op.binary(self, other, "div")

	def __xor__(self, other):
		return Op.binary(self, other, "xor")

	def __neg__(self):
		return Op.unary(self, "neg")
	
	def __radd__(self, other):
		if other==0.0: return self	
		return Op.binary(self, other, "add", reverse=True)

	def __rsub__(self, other):
		if other==0.0: return Op.unary(self, "neg")	
		return Op.binary(self, other, "sub", reverse=True)

	def __rmul__(self, other):
		if other==0.0: return 0
		if other==1.0: return self
		return Op.binary(self, other, "mul", reverse=True)

	def __rdiv__(self, other):
		if other==0.0: return 0
		return Op.binary(self, other, "div", reverse=True)

	def __iadd__(self, other):
		if  self.array_flag and self.array.horiz_add_flag==+1:
			if self.array.horiz_add_index>int(self.array_index):
				raise ValueError, "access to array has to be monotone"
			for i in range (self.array.horiz_add_index+1, int(self.array_index)+1):
				dummy = Op("_HA_%s_%i" % (repr(self.array), self.array.horiz_add_modulo_index)).__ixor__(0)
				self.array.horiz_add_index += 1
				self.array.horiz_add_modulo_index += 1
				if self.array.horiz_add_modulo_index == Op.templates["short_vector_length"]:
					i = self.array.horiz_add_index - Op.templates["short_vector_length"]
					print indent*'\t' + (self.array.horiz_add_template % (i, i )) + ';'
					self.array.horiz_add_modulo_index = 0
			dummy = Op("_HA_%s_%i" % (repr(self.array), self.array.horiz_add_modulo_index)).__ixor__(other)
			self.array.horiz_add_index += 1
			self.array.horiz_add_modulo_index += 1
			if self.array.horiz_add_modulo_index == Op.templates["short_vector_length"]:
				i = self.array.horiz_add_index - Op.templates["short_vector_length"]
				print indent*'\t' + (self.array.horiz_add_template % (i, i )) + ';'

				self.array.horiz_add_modulo_index = 0
			return None
		if  self.array_flag and self.array.horiz_add_flag==-1:
			if self.array.horiz_add_index <0:
				self.array.horiz_add_index = Op.templates["short_vector_length"]*int(ceil(int(self.array_index)/(Op.templates["short_vector_length"]+0.0)))
				self.array.horiz_add_modulo_index = int(self.array_index)%Op.templates["short_vector_length"]+1 
			elif self.array.horiz_add_index<int(self.array_index):
				raise ValueError, "access to array has to be monotone"
			for i in range (self.array.horiz_add_index-1, int(self.array_index)-1):
				dummy = Op("_HA_%s_%i" % (repr(self.array), self.array.horiz_add_modulo_index)).__ixor__(0)
				self.array.horiz_add_index -= 1
				self.array.horiz_add_modulo_index -= 1
				if self.array.horiz_add_modulo_index == 0:
					i = self.array.horiz_add_index #- Op.templates["short_vector_length"]
					print indent*'\t' + (self.array.horiz_add_template % (i, i )) + ';'
					self.array.horiz_add_modulo_index = Op.templates["short_vector_length"]
			dummy = Op("_HA_%s_%i" % (repr(self.array), self.array.horiz_add_modulo_index-1)).__ixor__(other)
			self.array.horiz_add_index -= 1
			self.array.horiz_add_modulo_index -= 1
			if self.array.horiz_add_modulo_index == 0:
				i = self.array.horiz_add_index #- Op.templates["short_vector_length"]
				print indent*'\t' + (self.array.horiz_add_template % (i, i )) + ';'

				self.array.horiz_add_modulo_index = Op.templates["short_vector_length"]
			return None

		return self.__ixor__(self.__add__(other))
	
	def __isub__(self, other):
		return self.__ixor__(self.__sub__(other))
	
	def __imul__(self, other):
		return self.__ixor__(self.__mul__(other))
	
	def __idiv__(self, other):
		return self.__ixor__(self.__div__(other))

	def __getitem__(self, index):
		if self.base_array_flag:
			ret = Op(Op.templates["getbaseitem"] % (repr(self), repr(index)),
				base_array_flag = True)
		else:
			if self.horiz_add_flag==0:
				index *= Op.templates["short_vector_length"]
			ret = Op(Op.templates["getitem"] % (repr(self), repr(index)),
				array_flag = True)
		ret.array = self
		ret.array_index = repr(index)
		return ret

	def __setitem__(self, index, value):
		return self

def set_vec(x0, x1=None, x2=None, x3=None):
	if Op.templates["short_vector_length"] == 1:
		return Op(Op.templates["set"] % repr(x0))
	if Op.templates["short_vector_length"] == 2:
		return Op(Op.templates["set"] % (repr(x0), repr(x1)))
	if Op.templates["short_vector_length"] == 4:
		return Op(Op.templates["set"] % (repr(x0), repr(x1), repr(x2), repr(x3)))

def sqrt(x):
	return Op(Op.templates["sqrt"] % repr(x)) 
	
def parameter(name, array=False, base_array=False, horiz_add=0):
	if Op.templates["short_vector_length"] == 1:
		horiz_add = False
	if horiz_add:
		base_array = False
		array = True
		for i in range(Op.templates["short_vector_length"]):
			print "%s%s _HA_%s_%i = %s;" %(indent*'\t', Op.templates["type"], name, i, Op.templates["zero"])
	return Op(name, array_flag=array, base_array_flag=base_array, horiz_add_flag=horiz_add)

def var(name, init=None, array_length=0, base_array=False):
	if array_length>0:
		print "%s%s %s[%s];" % (indent*'\t', Op.templates["type"], name, array_length)
	elif init!=None:
		if not isinstance(init, Op): 
			if init==0:
				init = Op(Op.templates["zero"])
			else:	
				init = Op(Op.templates["set1"] % repr(init))
		print "%s%s %s = %s;" % (indent*'\t', Op.templates["type"], name, init)
	else:
		print "%s%s %s;" % (indent*'\t', Op.templates["type"], name)
	return  Op(name, array_flag=(array_length>0), base_array_flag=base_array)

def fini_horiz_add(M):
	if M.horiz_add_flag != +1:
		return
	s = "h[0]"
	for i in range(1, Op.templates["short_vector_length"]):
		s += " + h[%i]" % i
	if M.horiz_add_modulo_index == 1:
		print "%s{%s *h = (%s *) &_HA_%s_0; %s[%i] += %s;}" \
		% (indent*'\t', Op.templates["basetype"], Op.templates["basetype"], M.s, M.s, M.horiz_add_index-1, s)
	elif M.horiz_add_modulo_index == 2:
		print "%s{%s *h = (%s *) &_HA_%s_0; %s[%i] += %s;" \
		% (indent*'\t', Op.templates["basetype"], Op.templates["basetype"], M.s, M.s, M.horiz_add_index-2, s)
		print "%s h = (%s *) &_HA_%s_1; %s[%i] += %s;}" \
		% (indent*'\t', Op.templates["basetype"], M.s, M.s, M.horiz_add_index-1, s)
	elif M.horiz_add_modulo_index == 3:
		print "%s{%s *h = (%s *) &_HA_%s_0; %s[%i] += %s;" \
		% (indent*'\t', Op.templates["basetype"], Op.templates["basetype"], M.s, M.s, M.horiz_add_index-3, s)
		print "%s h = (%s *) &_HA_%s_1; %s[%i] += %s;" \
		% (indent*'\t', Op.templates["basetype"], M.s, M.s, M.horiz_add_index-2, s)
		print "%s h = (%s *) &_HA_%s_2; %s[%i] += %s;}" \
		% (indent*'\t', Op.templates["basetype"], M.s, M.s, M.horiz_add_index-1, s)
	

