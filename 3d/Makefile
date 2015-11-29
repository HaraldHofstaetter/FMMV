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

LIB_INSTALLDIR = $(PWD)/../lib
INC_INSTALLDIR = $(PWD)/../include
LIBNAME = libfmmv3d
LIBFNAME = libfmmv3df

all: lib 

single: Makefile2 $(ODIR)/$(LIBFNAME).la

double: Makefile2 $(ODIR)/$(LIBNAME).la

lib: single double fmmv3d.h
	if test ! -d $(LIB_INSTALLDIR)  ; then mkdir $(LIB_INSTALLDIR) ; fi
	$(LIBTOOL_INSTALL) $(ODIR)/$(LIBFNAME).la $(LIB_INSTALLDIR)
	$(LIBTOOL_INSTALL) $(ODIR)/$(LIBNAME).la $(LIB_INSTALLDIR)
	$(LIBTOOL_FINISH) $(LIB_INSTALLDIR)
	if test ! -d $(INC_INSTALLDIR)  ; then mkdir $(INC_INSTALLDIR) ; fi
	$(LIBTOOL_INSTALL) fmmv3d.h $(INC_INSTALLDIR)

#dll: $(ODIR)/fmmvcoul3d.dll $(ODIR)/fmmvcoul3df.dll 
dll: $(ODIR)/fmmvcoul3d.dll 


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

#$(ODIR)/$(LIBNAME).la: odir $(OBJS_FMMV)
#	$(LIBTOOL_LINK) -o \
#	$(ODIR)/$(LIBNAME).la $(LOBJS_FMMV) \
#	-rpath $(LIB_INSTALLDIR) $(LIBS) 
#
#$(ODIR)/$(LIBFNAME).la: odir $(OBJS_FMMV_SINGLE) 
#	$(LIBTOOL_LINK) -o \
#	$(ODIR)/$(LIBFNAME).la $(LOBJS_FMMV_SINGLE) \
#	-rpath $(LIB_INSTALLDIR) $(LIBS)



### generated source files: ####


fmmv3d.h: ../gen_frontend.py fmmv_setup.py
	PYTHONPATH=".:..:$(PYTHONPATH)" python ../gen_frontend.py header > fmmv3d.h

fmmv3d.c: ../gen_frontend.py fmmv_setup.py fmmv3d.h
	PYTHONPATH=".:..:$(PYTHONPATH)" python ../gen_frontend.py source > fmmv3d.c

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

M2L.h: ../gen_M2L_h.py
	PYTHONPATH=".:..:$(PYTHONPATH)" python ../gen_M2L_h.py > M2L.h

M2L.c: M2L.h
	ln -s -f ../M2L.c M2L.c
	
core_gen_eval.c: gen_core_gen_eval.py ../gen_straightline_code.py
	PYTHONPATH="..:.:$(PYTHONPATH)" python gen_core_gen_eval.py > core_gen_eval.c
	
core_gen_eval_simd2.c: gen_core_gen_eval.py ../gen_straightline_code.py
	PYTHONPATH="..:.:$(PYTHONPATH)" python gen_core_gen_eval.py simd2 > core_gen_eval_simd2.c

core_gen_eval_simd4.c: gen_core_gen_eval.py ../gen_straightline_code.py
	PYTHONPATH="..:.:$(PYTHONPATH)" python gen_core_gen_eval.py simd4 > core_gen_eval_simd4.c

spherical_harmonics.c: gen_spherical_harmonics.py ../gen_straightline_code.py
	PYTHONPATH="..:.:$(PYTHONPATH)" python gen_spherical_harmonics.py > spherical_harmonics.c

spherical_harmonics_simd2.c: gen_spherical_harmonics.py ../gen_straightline_code.py
	PYTHONPATH="..:.:$(PYTHONPATH)" python gen_spherical_harmonics.py simd2 > spherical_harmonics_simd2.c

spherical_harmonics_simd4.c: gen_spherical_harmonics.py ../gen_straightline_code.py
	PYTHONPATH="..:.:$(PYTHONPATH)" python gen_spherical_harmonics.py simd4 > spherical_harmonics_simd4.c

FFT_M2X_X2L_cores_init.c: gen_FFT_M2X_X2L.py ../gen_straightline_code.py
	PYTHONPATH="..:.:$(PYTHONPATH)" python gen_FFT_M2X_X2L.py cores_init > FFT_M2X_X2L_cores_init.c

FFT_M2X_X2L_cores.c: gen_FFT_M2X_X2L.py ../gen_straightline_code.py
	PYTHONPATH="..:.:$(PYTHONPATH)" python gen_FFT_M2X_X2L.py cores > FFT_M2X_X2L_cores.c

FFT_M2X_X2L_cores_simd4.c: gen_FFT_M2X_X2L.py ../gen_straightline_code.py
	PYTHONPATH="..:.:$(PYTHONPATH)" python gen_FFT_M2X_X2L.py cores_simd4 > FFT_M2X_X2L_cores_simd4.c

FFT_M2X_X2L_cores_simd2.c: gen_FFT_M2X_X2L.py ../gen_straightline_code.py
	PYTHONPATH="..:.:$(PYTHONPATH)" python gen_FFT_M2X_X2L.py cores_simd2 > FFT_M2X_X2L_cores_simd2.c

FFT_M2X_X2L_aux.c: gen_FFT_M2X_X2L_aux.py ../gen_straightline_code.py
	PYTHONPATH="..:.:$(PYTHONPATH)" python  gen_FFT_M2X_X2L_aux.py > FFT_M2X_X2L_aux.c

FFT_M2X_X2L_aux_simd4.c: gen_FFT_M2X_X2L_aux.py ../gen_straightline_code.py
	PYTHONPATH="..:.:$(PYTHONPATH)" python gen_FFT_M2X_X2L_aux.py simd4 > FFT_M2X_X2L_aux_simd4.c

FFT_M2X_X2L_aux_simd2.c: gen_FFT_M2X_X2L_aux.py ../gen_straightline_code.py
	PYTHONPATH="..:.:$(PYTHONPATH)" python gen_FFT_M2X_X2L_aux.py simd2 > FFT_M2X_X2L_aux_simd2.c

FFT_M2X_X2L_tables.h: gen_FFT_M2X_X2L_tables.py ../gen_straightline_code.py
	PYTHONPATH="..:.:$(PYTHONPATH)" python gen_FFT_M2X_X2L_tables.py > FFT_M2X_X2L_tables.h

FFT_M2X_X2L_tables_simd4.h: gen_FFT_M2X_X2L_tables.py ../gen_straightline_code.py
	PYTHONPATH="..:.:$(PYTHONPATH)" python gen_FFT_M2X_X2L_tables.py simd4 > FFT_M2X_X2L_tables_simd4.h

FFT_M2X_X2L_tables_simd2.h: gen_FFT_M2X_X2L_tables.py ../gen_straightline_code.py
	PYTHONPATH="..:.:$(PYTHONPATH)" python gen_FFT_M2X_X2L_tables.py simd2 > FFT_M2X_X2L_tables_simd2.h

FFT_M2X_X2L.c: FFT_M2X_X2L_tables.h

FFT_M2X_X2L_simd2.c: FFT_M2X_X2L_tables_simd2.h

FFT_M2X_X2L_simd4.c: FFT_M2X_X2L_tables_simd4.h

GENERATED_FILES = Makefile2 \
	spherical_harmonics.c \
	spherical_harmonics_simd2.c \
	spherical_harmonics_simd4.c \
	core_gen_eval.c \
	core_gen_eval_simd2.c \
	core_gen_eval_simd4.c \
	M2L.h \
	FFT_M2X_X2L_cores_init.c \
	FFT_M2X_X2L_cores.c \
	FFT_M2X_X2L_cores_simd4.c \
	FFT_M2X_X2L_cores_simd2.c \
	FFT_M2X_X2L_aux.c \
	FFT_M2X_X2L_aux_simd4.c \
	FFT_M2X_X2L_aux_simd2.c \
	FFT_M2X_X2L_tables.h \
	FFT_M2X_X2L_tables_simd4.h \
	FFT_M2X_X2L_tables_simd2.h 

NAMES = \
	fmmv3d	\
        direct_method_xp \
	ida \
	fmmv_common 	\
	build_tree 	\
	create_box 	\
	gen_adaptive_lists 	\
	gen_adaptive_lists_st 	\
	generic_fmm generic_fmm_ws2 	\
	gen_neighbor_relations	\
	gen_neighbor2_list	\
	M2L statistics 	\
	periodic_lattice_M2L	\
	gen_eval	\
	core_gen_eval	\
	bessel \
	spherical_harmonics	\
	spherical_harmonics_generic	\
	init_all	\
	init_coeffs	\
	quad_coeffs	\
	Ry_init	\
	M2M_L2L_init	\
	M2L_init	\
	perm_init	\
	FFT_M2X_X2L_init	\
	FFT_M2X_X2L_cores_init

NAMES_SIMD = \
	ida_copy \
	vec_ops \
	gen_eval	\
	core_gen_eval	\
	bessel \
        math \
	spherical_harmonics	\
	spherical_harmonics_generic	\
	M2M_L2L	\
	M2X_X2L	\
	Ry	\
	Rz	\
	FFT_M2X_X2L_cores	\
	FFT_M2X_X2L	\
	FFT_M2X_X2L_aux	

NAMES_NOSIMD = \
	ida_copy \
	vec_ops \
	M2M_L2L	\
	M2X_X2L	\
	Ry	\
	Rz	\
	FFT_M2X_X2L_cores	\
	FFT_M2X_X2L	\
	FFT_M2X_X2L_aux	



fmmv_common.c build_tree.c create_box.c gen_adaptive_lists.c \
   gen_adaptive_lists_st.c generic_fmm.c generic_fmm_ws2.c \
   gen_neighbor_relations.c statistics.c \
   vec_ops.h vec_ops.c vec_ops_simd2.c vec_ops_simd4.c : 
	ln -s -f ../$@ $@

$(ODIR)/%_float.o: %.c
	$(LIBTOOL_COMPILE) -c -DFMM_PRECISION=0 $(CFLAGS) -I..  $< -o $@

$(ODIR)/%_double.o: %.c
	$(LIBTOOL_COMPILE) -c -DFMM_PRECISION=1 $(CFLAGS) -I..  $< -o $@


OBJS_FLOAT = $(addprefix $(ODIR)/, $(addsuffix _float.o, $(NAMES))) \
	     $(addprefix $(ODIR)/, $(addsuffix _simd4_float.o, $(NAMES_ACCESS) $(NAMES_SIMD)))

LOBJS_FLOAT = $(addprefix $(ODIR)/, $(addsuffix _float.lo, $(NAMES))) \
	      $(addprefix $(ODIR)/, $(addsuffix _simd4_float.lo, $(NAMES_ACCESS) $(NAMES_SIMD)))

OBJS_DOUBLE = $(addprefix $(ODIR)/, $(addsuffix _double.o, $(NAMES))) \
	      $(addprefix $(ODIR)/, $(addsuffix _simd2_double.o, $(NAMES_ACCESS) $(NAMES_SIMD)))

LOBJS_DOUBLE = $(addprefix $(ODIR)/, $(addsuffix _double.lo, $(NAMES))) \
	       $(addprefix $(ODIR)/, $(addsuffix _simd2_double.lo, $(NAMES_ACCESS) $(NAMES_SIMD)))

#OBJS_FLOAT = $(addprefix $(ODIR)/, $(addsuffix _float.o, $(NAMES))) \
#	     $(addprefix $(ODIR)/, $(addsuffix _float.o, $(NAMES_ACCESS) $(NAMES_NOSIMD)))
#
#LOBJS_FLOAT = $(addprefix $(ODIR)/, $(addsuffix _float.lo, $(NAMES))) \
#	      $(addprefix $(ODIR)/, $(addsuffix _float.lo, $(NAMES_ACCESS) $(NAMES_NOSIMD)))
#
#OBJS_DOUBLE = $(addprefix $(ODIR)/, $(addsuffix _double.o, $(NAMES))) \
#	      $(addprefix $(ODIR)/, $(addsuffix _double.o, $(NAMES_ACCESS) $(NAMES_NOSIMD)))
#
#LOBJS_DOUBLE = $(addprefix $(ODIR)/, $(addsuffix _double.lo, $(NAMES))) \
#	       $(addprefix $(ODIR)/, $(addsuffix _double.lo, $(NAMES_ACCESS) $(NAMES_NOSIMD)))


$(OBJS_FLOAT) $(OBJS_DOUBLE): fmmv_sys.h fmmv_access.h 

$(ODIR)/$(LIBNAME).la: odir $(OBJS_DOUBLE) 
	$(LIBTOOL_LINK) -o \
	$(ODIR)/$(LIBNAME).la $(LOBJS_DOUBLE)  \
	-rpath $(LIB_INSTALLDIR) $(LIBS) 

$(ODIR)/$(LIBFNAME).la: odir $(OBJS_FLOAT) 
	$(LIBTOOL_LINK) -o \
	$(ODIR)/$(LIBFNAME).la $(LOBJS_FLOAT) \
	-rpath $(LIB_INSTALLDIR) $(LIBS)


