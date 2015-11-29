# FMMV - the Fastest Multipole Method of Vienna
# Copyright (c) 2006 Harald Hofstaetter
# Institute for Analysis and Scientific Computing
# Vienna University of Technology
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

##############################################################################


include ../Makefile.inc

INSTALLDIR = $(PWD)/../lib
INC_INSTALLDIR = $(PWD)/../include
LIBNAME = libfmmv2d
LIBFNAME = libfmmv2df

all: lib 

single: Makefile2 $(ODIR)/$(LIBFNAME).la

double: Makefile2 $(ODIR)/$(LIBNAME).la

lib: single double fmmv2d.h
	if test ! -d $(LIB_INSTALLDIR)  ; then mkdir $(LIB_INSTALLDIR) ; fi
	$(LIBTOOL_INSTALL) $(ODIR)/$(LIBFNAME).la $(INSTALLDIR)
	$(LIBTOOL_INSTALL) $(ODIR)/$(LIBNAME).la $(INSTALLDIR)
	$(LIBTOOL_FINISH) $(INSTALLDIR)
	if test ! -d $(INC_INSTALLDIR)  ; then mkdir $(INC_INSTALLDIR) ; fi
	$(LIBTOOL_INSTALL) fmmv2d.h $(INC_INSTALLDIR)



#dll: $(ODIR)/fmmvcoul2d.dll $(ODIR)/fmmvcoul2df.dll 
dll: $(ODIR)/fmmvcoul2d.dll 

clean:
	rm -f $(ODIR)/*.o
	rm -f $(ODIR)/*.a
	rm -f $(ODIR)/*.lo
	rm -f $(ODIR)/*.la
	rm -f -r $(ODIR)/.libs

odir:
	if test ! -d $(ODIR)  ; then mkdir $(ODIR) ; fi



Makefile2: ../gen_Makefile2.py fmmv_setup.py
	PYTHONPATH=".:..:$(PYTHONPATH)" python ../gen_Makefile2.py  > Makefile2

include Makefile2



### generated source files: ####


fmmv2d.h: ../gen_frontend.py fmmv_setup.py
	PYTHONPATH=".:..:$(PYTHONPATH)" python ../gen_frontend.py header > fmmv2d.h

fmmv2d.c: ../gen_frontend.py fmmv_setup.py fmmv2d.h
	PYTHONPATH=".:..:$(PYTHONPATH)" python ../gen_frontend.py source > fmmv2d.c

fmmv_sys.h: ../gen_ida.py fmmv_setup.py
	PYTHONPATH=".:..:$(PYTHONPATH)" python ../gen_ida.py header > fmmv_sys.h

fmmv_access.h: ../gen_ida.py fmmv_setup.py
	PYTHONPATH=".:..:$(PYTHONPATH)" python ../gen_ida.py access > fmmv_access.h

ida.c: ../gen_ida.py fmmv_setup.py
	PYTHONPATH=".:..:$(PYTHONPATH)" python ../gen_ida.py source > ida.c

ida_copy.c: ../gen_ida_copy.py fmmv_setup.py
	PYTHONPATH=".:..:$(PYTHONPATH)" python ../gen_ida_copy.py > ida_copy.c

ida_copy_simd2.c: ../gen_ida_copy.py fmmv_setup.py
	PYTHONPATH=".:..:$(PYTHONPATH)" python ../gen_ida_copy.py simd2 > ida_copy_simd2.c

ida_copy_simd4.c: ../gen_ida_copy.py fmmv_setup.py
	PYTHONPATH=".:..:$(PYTHONPATH)" python ../gen_ida_copy.py simd4 > ida_copy_simd4.c

gen_neighbor2_list.c: ../gen_gen_neighbor2_list.py 
	PYTHONPATH=".:..:$(PYTHONPATH)" python ../gen_gen_neighbor2_list.py > gen_neighbor2_list.c

M2L.c: gen_M2L.py
	PYTHONPATH=".:..:$(PYTHONPATH)" python gen_M2L.py > M2L.c

#M2L.c: M2L.h


GENERATED_FILES = Makefile2 

NAMES = \
	fmmv2d	\
        direct_method_xp \
	ida		\
	fmmv_common 	\
	build_tree 	\
	create_box 	\
	gen_adaptive_lists 	\
	gen_adaptive_lists_st 	\
	generic_fmm generic_fmm_ws2 	\
	gen_neighbor_relations	\
	gen_neighbor2_list	\
	M2L		\
	statistics 	\
	periodic_lattice_M2L	\
	M2L_init	\
	init_all	\
	quad_coeffs	\
	M2L_aux		\
	bessel		\
	gen_eval


NAMES_SIMD = \
	ida_copy \
	vec_ops \
	M2M	\
	L2L	\
	M2X	\
	X2L	\


fmmv_common.c build_tree.c create_box.c gen_adaptive_lists.c \
   gen_adaptive_lists_st.c generic_fmm.c generic_fmm_ws2.c \
   gen_neighbor_relations.c statistics.c \
   vec_ops.h vec_ops.c vec_ops_simd2.c vec_ops_simd4.c : 
	ln -s -f ../$@ $@

ida_copy.c ida_copy_simd2.c ida_copy_simd4.c: fmmv_access.h

$(ODIR)/%_float.o: %.c
	$(LIBTOOL_COMPILE) -c -DFMM_PRECISION=0 $(CFLAGS) -I..  $< -o $@

$(ODIR)/%_double.o: %.c
	$(LIBTOOL_COMPILE) -c -DFMM_PRECISION=1 $(CFLAGS) -I..  $< -o $@


#OBJS_FLOAT = $(addprefix $(ODIR)/, $(addsuffix _float.o, $(NAMES))) \
#	     $(addprefix $(ODIR)/, $(addsuffix _simd4_float.o, $(NAMES_ACCESS) $(NAMES_SIMD)))
#
#LOBJS_FLOAT = $(addprefix $(ODIR)/, $(addsuffix _float.lo, $(NAMES))) \
#	      $(addprefix $(ODIR)/, $(addsuffix _simd4_float.lo, $(NAMES_ACCESS) $(NAMES_SIMD)))
#
#OBJS_DOUBLE = $(addprefix $(ODIR)/, $(addsuffix _double.o, $(NAMES))) \
#	      $(addprefix $(ODIR)/, $(addsuffix _simd2_double.o, $(NAMES_ACCESS) $(NAMES_SIMD)))
#
#LOBJS_DOUBLE = $(addprefix $(ODIR)/, $(addsuffix _double.lo, $(NAMES))) \
#	       $(addprefix $(ODIR)/, $(addsuffix _simd2_double.lo, $(NAMES_ACCESS) $(NAMES_SIMD)))
#

#NON_SIMD-VERSIONS:

OBJS_FLOAT = $(addprefix $(ODIR)/, $(addsuffix _float.o, $(NAMES) $(NAMES_ACCESS) $(NAMES_SIMD)))

LOBJS_FLOAT = $(addprefix $(ODIR)/, $(addsuffix _float.lo, $(NAMES) $(NAMES_ACCESS) $(NAMES_SIMD)))

OBJS_DOUBLE = $(addprefix $(ODIR)/, $(addsuffix _double.o, $(NAMES) $(NAMES_ACCESS) $(NAMES_SIMD)))

LOBJS_DOUBLE = $(addprefix $(ODIR)/, $(addsuffix _double.lo, $(NAMES) $(NAMES_ACCESS) $(NAMES_SIMD)))

$(OBJS_FLOAT) $(OBJS_DOUBLE): fmmv_sys.h fmmv_access.h 

$(ODIR)/$(LIBNAME).la: odir $(OBJS_DOUBLE) 
	$(LIBTOOL_LINK) -o \
	$(ODIR)/$(LIBNAME).la $(LOBJS_DOUBLE)  \
	-rpath $(INSTALLDIR) $(LNFLAGS) 

$(ODIR)/$(LIBFNAME).la: odir $(OBJS_FLOAT) 
	$(LIBTOOL_LINK) -o \
	$(ODIR)/$(LIBFNAME).la $(LOBJS_FLOAT) \
	-rpath $(INSTALLDIR) $(LNFLAGS)


