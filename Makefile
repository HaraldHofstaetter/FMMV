all:
	cd 2d; make 
	cd 3d; make

install: all
	cd 2d; make install
	cd 3d; make install

clean:
	cd 2d; make clean
	cd 3d; make clean
	cd test; make clean
	cd python; make clean
	cd matlab; make clean
	cd fortran; make clean

tgz: fmmv.tgz

TGZ_FILES = \
    include/fmmv2d.h \
    include/fmmv3d.h \
    build_tree.c \
    create_box.c \
    determine_direct_eval_treshold.c \
    fmmv_common.c \
    fmmv_common.h \
    gen_adaptive_lists.c \
    gen_adaptive_lists_st.c \
    generic_fmm.c \
    generic_fmm_ws2.c \
    gen_frontend.py \
    gen_gen_neighbor2_list.py \
    gen_ida_copy.py \
    gen_ida.py \
    gen_M2L_h_2d.py \
    gen_M2L_h.py \
    gen_Makefile2.py \
    gen_neighbor_relations.c \
    gen_straightline_code.py \
    M2L.c \
    Makefile \
    Makefile.inc \
    Makefile.gcc \
    Makefile.osx \
    pyparsing.py \
    simd2d.h \
    simd2d_sse2.h \
    simd4s_altivec.h \
    simd4s.h \
    simd4s_sse.h \
    simd_altivec.h \
    simd.h \
    simd_sse.h \
    statistics.c \
    timeval.h \
    vec_ops.c \
    vec_ops_simd2.c \
    vec_ops_simd4.c \
    2d/build_tree.c \
    2d/create_box.c \
    2d/direct_method_xp.c \
    2d/M2L_init.c \
    2d/gen_eval_expansions_base.c \
    2d/gen_eval_expansions.c \
    2d/gen_eval_expansions_base_simd2.c \
    2d/gen_eval_expansions_simd2.c \
    2d/gen_eval_expansions_base_simd4.c \
    2d/gen_eval_expansions_simd4.c \
    2d/eval_direct.c \
    2d/eval_direct_simd2.c \
    2d/eval_direct_simd4.c \
    2d/extrinsic_correction.c \
    2d/_fmmv.h \
    2d/fmmv_access.h \
    2d/fmmv_common.c \
    2d/fmmv2d.c \
    2d/fmmv2d.h \
    2d/fmmv_sys.h \
    2d/fmmv_setup.py \
    2d/gen_adaptive_lists.c \
    2d/gen_adaptive_lists_st.c \
    2d/generic_fmm.c \
    2d/generic_fmm_ws2.c \
    2d/gen_M2L.py \
    2d/gen_neighbor2_list.c \
    2d/gen_neighbor_relations.c \
    2d/ida.c \
    2d/ida_copy.c \
    2d/init_all.c \
    2d/Makefile \
    2d/Makefile2 \
    2d/L2L.c \
    2d/M2L_aux.c \
    2d/M2L.c \
    2d/M2M.c \
    2d/M2X.c \
    2d/periodic_lattice_M2L.c \
    2d/quad_coeffs.c \
    2d/bessel.c \
    2d/statistics.c \
    2d/vec_ops.c \
    2d/X2L.c \
    3d/build_tree.c \
    3d/create_box.c \
    3d/direct_method_xp.c \
    3d/gen_eval_expansions_base.c \
    3d/gen_eval_expansions.c \
    3d/gen_eval_expansions_base_simd2.c \
    3d/gen_eval_expansions_simd2.c \
    3d/gen_eval_expansions_base_simd4.c \
    3d/gen_eval_expansions_simd4.c \
    3d/eval_direct.c \
    3d/eval_direct_simd2.c \
    3d/eval_direct_simd4.c \
    3d/extrinsic_correction.c \
    3d/extrinsic_correction_simd2.c \
    3d/extrinsic_correction_simd4.c \
    3d/FFT_M2X_X2L.c \
    3d/FFT_M2X_X2L_init.c \
    3d/FFT_M2X_X2L_simd2.c \
    3d/FFT_M2X_X2L_simd4.c \
    3d/_fmmv.h \
    3d/fmmv3d.h \
    3d/fmmv_sys.h \
    3d/fmmv_common.c \
    3d/fmmv_setup.py \
    3d/gen_adaptive_lists.c \
    3d/gen_adaptive_lists_st.c \
    3d/gen_gen_eval_core.py \
    3d/generic_fmm.c \
    3d/generic_fmm_ws2.c \
    3d/gen_FFT_M2X_X2L_aux.py \
    3d/gen_FFT_M2X_X2L.py \
    3d/gen_FFT_M2X_X2L_tables.py \
    3d/gen_neighbor_relations.c \
    3d/gen_spherical_harmonics.py \
    3d/init_all.c \
    3d/init_coeffs.c \
    3d/Makefile \
    3d/Makefile2 \
    3d/M2L.c \
    3d/M2M_L2L.c \
    3d/M2M_L2L_simd2.c \
    3d/M2M_L2L_simd4.c \
    3d/M2X_X2L.c \
    3d/M2X_X2L_simd2.c \
    3d/M2X_X2L_simd4.c \
    3d/periodic_lattice_M2L.c \
    3d/perm_init.c \
    3d/quad_coeffs.c \
    3d/Ry.c \
    3d/Ry_init.c \
    3d/Ry_simd2.c \
    3d/Ry_simd4.c \
    3d/Rz.c \
    3d/Rz_simd2.c \
    3d/Rz_simd4.c \
    3d/spherical_harmonics_generic.c \
    3d/spherical_harmonics_generic.c \
    3d/spherical_harmonics_generic_simd2.c \
    3d/spherical_harmonics_generic_simd4.c \
    3d/math_simd2.c \
    3d/math_simd4.c \
    3d/bessel.c \
    3d/bessel_simd2.c \
    3d/bessel_simd4.c \
    3d/statistics.c \
    3d/M2M_L2L_init.c \
    3d/M2L_init.c \
    3d/vec_ops_simd2.c \
    3d/vec_ops_simd4.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_M2X_12.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_M2X_16.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_M2X_20.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_M2X_24.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_M2X_28.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_M2X_2.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_M2X_32.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_M2X_36.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_M2X_40.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_M2X_44.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_M2X_48.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_M2X_4.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_M2X_52.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_M2X_56.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_M2X_60.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_M2X_64.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_M2X_8.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_M2X_X2L.g \
    3d/SPIRAL_FFT_M2X_X2L/FFT_X2L_12.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_X2L_16.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_X2L_20.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_X2L_24.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_X2L_28.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_X2L_2.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_X2L_32.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_X2L_36.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_X2L_40.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_X2L_44.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_X2L_48.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_X2L_4.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_X2L_52.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_X2L_56.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_X2L_60.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_X2L_64.c \
    3d/SPIRAL_FFT_M2X_X2L/FFT_X2L_8.c \
    test/Makefile \
    test/test_fmmv.c \
    test/test_utilities.c \
    test/test_utilities.h \
    test/test_M_L_3d.c \
    test/test_M_L_2d.c \
    test/test_X_3d.c \
    test/test_X_2d.c \
    test/test_grad.c \
    python/Makefile \
    python/fmmv_demo.py \
    python/fmmvmodule.c \
    python/setup.py \
    python/test_fmmv.py  \
    python/aurora.gif \
    fortran/fmmv2d.f90  \
    fortran/fmmv3d.f90  \
    fortran/Makefile  \
    fortran/test_fmmv.F90  \
    fortran/test_utilities.F90\
    matlab/Makefile\
    matlab/fmmv.c           \
    matlab/fmmv_finalize.c    \
    matlab/fmmv_evaluate.c  \
    matlab/fmmv_initialize.c  \
    matlab/statistics.c\
    affmmlap/adaplapdriver.f  \
    affmmlap/lapadap.f       \
    affmmlap/laptable.f  \
    affmmlap/parm-alap.h  \
    affmmlap/second.f        \
    affmmlap/fmmadaplap.f     \
    affmmlap/functions.f  \
    affmmlap/lapoperators.f  \
    affmmlap/Makefile    \
    affmmlap/prini.f      \
    affmmlap/test_huang.f90  \
    affmmlap/treeadap.f \
    lib/dummy



fmmv.tgz: $(TGZ_FILES)
	tar czvf fmmv.tgz $(TGZ_FILES)

