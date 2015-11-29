/* This file is automatically generated by gen_ida.py */
/* DO NOT EDIT! */

#undef DEFINE_IDA_LOCAL_ALIASES
#undef access_x
#undef access_y
#undef access_tx
#undef access_ty
#undef access_mx
#undef access_my
#undef access_q
#undef access_pot
#undef access_gradx
#undef access_grady

#if (FMM_KIND == FMM_STANDARD)

#if (FMM_PRECISION==0)
#define DEFINE_IDA_LOCAL_ALIASES(FMMV) \
	_DATA_standard_V0_t *_0 = ((DATA_standard_t*)FMMV->DATA)->_0; \


#define access_x(i) (_0[i].x)
#define access_y(i) (_0[i].y)
#define access_q(i) (_0[i].q)
#define access_pot(i) (_0[i].pot)
#define access_tx(i) (_0[i].x)
#define access_ty(i) (_0[i].y)

#elif (FMM_PRECISION==1)
#define DEFINE_IDA_LOCAL_ALIASES(FMMV) \
	_DATA_standard_V0_t *_0 = ((DATA_standard_t*)FMMV->DATA)->_0; \


#define access_x(i) (_0[i].x)
#define access_y(i) (_0[i].y)
#define access_q(i) (_0[i].q)
#define access_pot(i) (_0[i].pot)
#define access_tx(i) (_0[i].x)
#define access_ty(i) (_0[i].y)

#endif /* FMM_PRECISION */

#endif /* FMM_KIND == FMM_STANDARD */


#if (FMM_KIND == FMM_DIPOLE)

#if (FMM_PRECISION==0)
#define DEFINE_IDA_LOCAL_ALIASES(FMMV) \
	_DATA_dipole_V0_t *_0 = ((DATA_dipole_t*)FMMV->DATA)->_0; \


#define access_x(i) (_0[i].x)
#define access_y(i) (_0[i].y)
#define access_q(i) (_0[i].q)
#define access_mx(i) (_0[i].mx)
#define access_my(i) (_0[i].my)
#define access_pot(i) (_0[i].pot)
#define access_tx(i) (_0[i].x)
#define access_ty(i) (_0[i].y)

#elif (FMM_PRECISION==1)
#define DEFINE_IDA_LOCAL_ALIASES(FMMV) \
	_DATA_dipole_V0_t *_0 = ((DATA_dipole_t*)FMMV->DATA)->_0; \


#define access_x(i) (_0[i].x)
#define access_y(i) (_0[i].y)
#define access_q(i) (_0[i].q)
#define access_mx(i) (_0[i].mx)
#define access_my(i) (_0[i].my)
#define access_pot(i) (_0[i].pot)
#define access_tx(i) (_0[i].x)
#define access_ty(i) (_0[i].y)

#endif /* FMM_PRECISION */

#endif /* FMM_KIND == FMM_DIPOLE */


#if (FMM_KIND == FMM_GRAD)

#if (FMM_PRECISION==0)
#define DEFINE_IDA_LOCAL_ALIASES(FMMV) \
	_DATA_grad_V0_t *_0 = ((DATA_grad_t*)FMMV->DATA)->_0; \


#define access_x(i) (_0[i].x)
#define access_y(i) (_0[i].y)
#define access_q(i) (_0[i].q)
#define access_pot(i) (_0[i].pot)
#define access_gradx(i) (_0[i].gradx)
#define access_grady(i) (_0[i].grady)
#define access_tx(i) (_0[i].x)
#define access_ty(i) (_0[i].y)

#elif (FMM_PRECISION==1)
#define DEFINE_IDA_LOCAL_ALIASES(FMMV) \
	_DATA_grad_V0_t *_0 = ((DATA_grad_t*)FMMV->DATA)->_0; \


#define access_x(i) (_0[i].x)
#define access_y(i) (_0[i].y)
#define access_q(i) (_0[i].q)
#define access_pot(i) (_0[i].pot)
#define access_gradx(i) (_0[i].gradx)
#define access_grady(i) (_0[i].grady)
#define access_tx(i) (_0[i].x)
#define access_ty(i) (_0[i].y)

#endif /* FMM_PRECISION */

#endif /* FMM_KIND == FMM_GRAD */


#if (FMM_KIND == FMM_DIPOLE_GRAD)

#if (FMM_PRECISION==0)
#define DEFINE_IDA_LOCAL_ALIASES(FMMV) \
	_DATA_dipole_grad_V0_t *_0 = ((DATA_dipole_grad_t*)FMMV->DATA)->_0; \


#define access_x(i) (_0[i].x)
#define access_y(i) (_0[i].y)
#define access_q(i) (_0[i].q)
#define access_mx(i) (_0[i].mx)
#define access_my(i) (_0[i].my)
#define access_pot(i) (_0[i].pot)
#define access_gradx(i) (_0[i].gradx)
#define access_grady(i) (_0[i].grady)
#define access_tx(i) (_0[i].x)
#define access_ty(i) (_0[i].y)

#elif (FMM_PRECISION==1)
#define DEFINE_IDA_LOCAL_ALIASES(FMMV) \
	_DATA_dipole_grad_V0_t *_0 = ((DATA_dipole_grad_t*)FMMV->DATA)->_0; \


#define access_x(i) (_0[i].x)
#define access_y(i) (_0[i].y)
#define access_q(i) (_0[i].q)
#define access_mx(i) (_0[i].mx)
#define access_my(i) (_0[i].my)
#define access_pot(i) (_0[i].pot)
#define access_gradx(i) (_0[i].gradx)
#define access_grady(i) (_0[i].grady)
#define access_tx(i) (_0[i].x)
#define access_ty(i) (_0[i].y)

#endif /* FMM_PRECISION */

#endif /* FMM_KIND == FMM_DIPOLE_GRAD */


#if (FMM_KIND == FMM_ST_STANDARD)

#if (FMM_PRECISION==0)
#define DEFINE_IDA_LOCAL_ALIASES(FMMV) \
	_DATA_standard_ST_V0_t *_0 = ((DATA_standard_ST_t*)FMMV->DATA)->_0; \
	_DATA_standard_ST_V1_t *_1 = ((DATA_standard_ST_t*)FMMV->DATA)->_1; \


#define access_x(i) (_0[i].x)
#define access_y(i) (_0[i].y)
#define access_q(i) (_0[i].q)
#define access_tx(i) (_1[i].tx)
#define access_ty(i) (_1[i].ty)
#define access_pot(i) (_1[i].pot)

#elif (FMM_PRECISION==1)
#define DEFINE_IDA_LOCAL_ALIASES(FMMV) \
	_DATA_standard_ST_V0_t *_0 = ((DATA_standard_ST_t*)FMMV->DATA)->_0; \
	_DATA_standard_ST_V1_t *_1 = ((DATA_standard_ST_t*)FMMV->DATA)->_1; \


#define access_x(i) (_0[i].x)
#define access_y(i) (_0[i].y)
#define access_q(i) (_0[i].q)
#define access_tx(i) (_1[i].tx)
#define access_ty(i) (_1[i].ty)
#define access_pot(i) (_1[i].pot)

#endif /* FMM_PRECISION */

#endif /* FMM_KIND == FMM_ST_STANDARD */


#if (FMM_KIND == FMM_ST_DIPOLE)

#if (FMM_PRECISION==0)
#define DEFINE_IDA_LOCAL_ALIASES(FMMV) \
	_DATA_dipole_ST_V0_t *_0 = ((DATA_dipole_ST_t*)FMMV->DATA)->_0; \
	_DATA_dipole_ST_V1_t *_1 = ((DATA_dipole_ST_t*)FMMV->DATA)->_1; \


#define access_x(i) (_0[i].x)
#define access_y(i) (_0[i].y)
#define access_q(i) (_0[i].q)
#define access_mx(i) (_0[i].mx)
#define access_my(i) (_0[i].my)
#define access_tx(i) (_1[i].tx)
#define access_ty(i) (_1[i].ty)
#define access_pot(i) (_1[i].pot)

#elif (FMM_PRECISION==1)
#define DEFINE_IDA_LOCAL_ALIASES(FMMV) \
	_DATA_dipole_ST_V0_t *_0 = ((DATA_dipole_ST_t*)FMMV->DATA)->_0; \
	_DATA_dipole_ST_V1_t *_1 = ((DATA_dipole_ST_t*)FMMV->DATA)->_1; \


#define access_x(i) (_0[i].x)
#define access_y(i) (_0[i].y)
#define access_q(i) (_0[i].q)
#define access_mx(i) (_0[i].mx)
#define access_my(i) (_0[i].my)
#define access_tx(i) (_1[i].tx)
#define access_ty(i) (_1[i].ty)
#define access_pot(i) (_1[i].pot)

#endif /* FMM_PRECISION */

#endif /* FMM_KIND == FMM_ST_DIPOLE */


#if (FMM_KIND == FMM_ST_GRAD)

#if (FMM_PRECISION==0)
#define DEFINE_IDA_LOCAL_ALIASES(FMMV) \
	_DATA_grad_ST_V0_t *_0 = ((DATA_grad_ST_t*)FMMV->DATA)->_0; \
	_DATA_grad_ST_V1_t *_1 = ((DATA_grad_ST_t*)FMMV->DATA)->_1; \


#define access_x(i) (_0[i].x)
#define access_y(i) (_0[i].y)
#define access_q(i) (_0[i].q)
#define access_tx(i) (_1[i].tx)
#define access_ty(i) (_1[i].ty)
#define access_pot(i) (_1[i].pot)
#define access_gradx(i) (_1[i].gradx)
#define access_grady(i) (_1[i].grady)

#elif (FMM_PRECISION==1)
#define DEFINE_IDA_LOCAL_ALIASES(FMMV) \
	_DATA_grad_ST_V0_t *_0 = ((DATA_grad_ST_t*)FMMV->DATA)->_0; \
	_DATA_grad_ST_V1_t *_1 = ((DATA_grad_ST_t*)FMMV->DATA)->_1; \


#define access_x(i) (_0[i].x)
#define access_y(i) (_0[i].y)
#define access_q(i) (_0[i].q)
#define access_tx(i) (_1[i].tx)
#define access_ty(i) (_1[i].ty)
#define access_pot(i) (_1[i].pot)
#define access_gradx(i) (_1[i].gradx)
#define access_grady(i) (_1[i].grady)

#endif /* FMM_PRECISION */

#endif /* FMM_KIND == FMM_ST_GRAD */


#if (FMM_KIND == FMM_ST_DIPOLE_GRAD)

#if (FMM_PRECISION==0)
#define DEFINE_IDA_LOCAL_ALIASES(FMMV) \
	_DATA_dipole_grad_ST_V0_t *_0 = ((DATA_dipole_grad_ST_t*)FMMV->DATA)->_0; \
	_DATA_dipole_grad_ST_V1_t *_1 = ((DATA_dipole_grad_ST_t*)FMMV->DATA)->_1; \


#define access_x(i) (_0[i].x)
#define access_y(i) (_0[i].y)
#define access_q(i) (_0[i].q)
#define access_mx(i) (_0[i].mx)
#define access_my(i) (_0[i].my)
#define access_tx(i) (_1[i].tx)
#define access_ty(i) (_1[i].ty)
#define access_pot(i) (_1[i].pot)
#define access_gradx(i) (_1[i].gradx)
#define access_grady(i) (_1[i].grady)

#elif (FMM_PRECISION==1)
#define DEFINE_IDA_LOCAL_ALIASES(FMMV) \
	_DATA_dipole_grad_ST_V0_t *_0 = ((DATA_dipole_grad_ST_t*)FMMV->DATA)->_0; \
	_DATA_dipole_grad_ST_V1_t *_1 = ((DATA_dipole_grad_ST_t*)FMMV->DATA)->_1; \


#define access_x(i) (_0[i].x)
#define access_y(i) (_0[i].y)
#define access_q(i) (_0[i].q)
#define access_mx(i) (_0[i].mx)
#define access_my(i) (_0[i].my)
#define access_tx(i) (_1[i].tx)
#define access_ty(i) (_1[i].ty)
#define access_pot(i) (_1[i].pot)
#define access_gradx(i) (_1[i].gradx)
#define access_grady(i) (_1[i].grady)

#endif /* FMM_PRECISION */

#endif /* FMM_KIND == FMM_ST_DIPOLE_GRAD */


