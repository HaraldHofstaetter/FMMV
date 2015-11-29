/*
 * FMMV - the Fastest Multipole Method of Vienna
 * Copyright (c) 2006 Harald Hofstaetter
 * Institute for Analysis and Scientific Computing
 * Vienna University of Technology
 * 
 * This file is part of FMMV.
 * 
 * FMMV is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * FMMV is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with FMMV; if not, write to the Free Software  Foundation, Inc.,
 * 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 * 
 */

#include "Python.h"

#if (FMM_DIM==2)
#include "fmmv2d.h"
#else
#include "fmmv3d.h"
#endif
#include "arrayobject.h"
#include <stdlib.h>

#ifdef USE_PAPI
	#include<papi.h>
#endif

#ifdef USE_SINGLE_PRECISION
	#define PyArray_MYTYPE PyArray_FLOAT
	//#define PyArray_MYTYPE PyArray_FLOAT32
	typedef float _FLOAT_;
#else
	#define PyArray_MYTYPE PyArray_DOUBLE
	//#define PyArray_MYTYPE PyArray_FLOAT64
	typedef double _FLOAT_;
#endif

static PyObject *ErrorObject;

static PyObject *
genStatistics(struct FmmvStatistics *statistics) 
{
	static char *stat_names[_STAT_LAST_] = {
		"total",
		"buildTree",
		"genM",
		"M2M",
		"M2L",
		"L2L",
		"evalL",
		"list1",
		"list3",
		"list4",
		"list34",
		"farfield",
		"nearfield",
		"initialize",
		"evaluate",
		"finalize",
	};	
	
	int i;
	PyObject* stat;
	PyObject* obj;

	stat = Py_BuildValue("{sisisisisisisisisisdsdsi}", 
			"noOfParticles", statistics->noOfParticles,
			"noOfTargets", statistics->noOfTargets,
			"maxNoOfStoredXin", statistics->maxNoOfStoredXin,	
			"noOfSourceLevels", statistics->noOfSourceLevels,
			"noOfTargetLevels", statistics->noOfTargetLevels,
			"noOfSourceBoxes", statistics->noOfSourceBoxes,
			"noOfTargetBoxes", statistics->noOfTargetBoxes,
			"noOfSourceLeafBoxes", statistics->noOfSourceLeafBoxes,
			"noOfTargetLeafBoxes", statistics->noOfTargetLeafBoxes,
			"averageNoOfParticlesPerLeafBox", (double) statistics->averageNoOfParticlesPerLeafBox,
			"averageNoOfTargetsPerLeafBox", (double) statistics->averageNoOfTargetsPerLeafBox,
			"maxAllocatedMemory", statistics->maxAllocatedMemory);

	obj = PyDict_New();
	for (i=0; i<_STAT_LAST_; i++) {
		if (statistics->time[i]>=-0.5) {
			PyDict_SetItemString(obj, stat_names[i], 
				PyFloat_FromDouble(statistics->time[i]));
		}	
  	}
	PyDict_SetItemString(stat, "time", obj);

	obj = PyDict_New();
	for (i=0; i<_STAT_LAST_; i++) {
		if (statistics->time[i]>=-0.5) {
			PyDict_SetItemString(obj, stat_names[i], 
				PyFloat_FromDouble(statistics->etime[i]));
		}	
  	}
	PyDict_SetItemString(stat, "etime", obj);
#ifdef USE_PAPI
	if (statistics->PAPIeventSet!=PAPI_NULL) {
		int ret, j;
		int events[MAX_NUM_PAPI_EVENTS];
		int NEvents=MAX_NUM_PAPI_EVENTS;
		char eventName[128];

		ret = PAPI_list_events(statistics->PAPIeventSet, events, &NEvents);
		for (i=0; i<NEvents; i++) {
			obj = PyDict_New();
			PAPI_event_code_to_name(events[i], eventName);
			for (j=0; j<_STAT_LAST_; j++) {
				if (statistics->time[j]>=-0.5) {
					PyDict_SetItemString(obj, stat_names[j], 
						PyLong_FromLongLong(statistics->PAPIvalues[j][i]));
				}	
			}	
			PyDict_SetItemString(stat, eventName, obj);
  		}
	}
#endif	
	
	
	return stat;
}	

#define FMM_S_EPS_MAX 100

static double user_defined_x[FMM_S_EPS_MAX];
static double user_defined_w[FMM_S_EPS_MAX];
static int user_defined_M[FMM_S_EPS_MAX];


static char py_fmmv__doc__[] ="";

static PyObject *
py_fmmv(PyObject *self, PyObject *args, PyObject *kwds)
{

	char *err = NULL;

	static char* argnames[] = {
		"particles", 
		"charges", 
		"dipoleMoments",
		"targets", 

		"computeGradient",
		"getStatistics",
		
		"precision",
		"beta",
		"scale",
		"splitThreshold",
		"splitTargetThreshold",
		"levels",
		"directEvalThreshold",
		"periodicBoundaryConditions",
		"extrinsicCorrection",
		"useHilbertOrder",
		"directEvalAccuracy",
		"useFarfieldNearfieldThreads",

		"p",
		"pM",
		"pL",
		"s",
		"ws",
		"reducedScheme",
		"x",
		"w",
		"M",
		
		"papi",
       		NULL };

	PyArrayObject *particles = NULL;
	PyArrayObject *charges = NULL;
	PyArrayObject *targets = NULL;
	PyArrayObject *dipoleMoments = NULL;

	struct FmmvOptions options;
	struct FmmvStatistics statistics;
	struct FmmvStatistics *statistics_p;

	int computeGradient = 0;
	int getStatistics = 0;
	
	int NParticles;
	PyArrayObject *particles1 = NULL;
	PyArrayObject *charges1 = NULL;
	PyArrayObject *dipoleMoments1 = NULL;
	
	int NTargets;
	PyArrayObject *targets1 = NULL;
	
	PyArrayObject *pot = NULL;
	PyArrayObject *grad = NULL;
	void *grad_data = NULL;
	_FLOAT_ *charges_data = NULL;
#if (FMM_DIM==2)
	_FLOAT_ (*targets_data)[2] = NULL;
	_FLOAT_ (*dipoleMoments_data)[2] = NULL;
#else
	_FLOAT_ (*targets_data)[3] = NULL;
	_FLOAT_ (*dipoleMoments_data)[3] = NULL;
#endif
	PyObject *stat = NULL;
	npy_intp dims[2];

	PyObject *papi = 0;
	int precision = -1;
	int p = -1;

	PyObject *x = 0;
	PyObject *w = 0;
	PyObject *M = 0;

#ifdef USE_PAPI
	int ret;
	int events[MAX_NUM_PAPI_EVENTS];
	int NEvents=0;
#endif
	

	options = fmmvGetDefaultOptions();

	if (!PyArg_ParseTupleAndKeywords(
		args, kwds, "O|OOO" "ii" "iddiiiiiiiii" "iiiiiiOOO" "O", argnames,
		&particles, 
		&charges,
		&dipoleMoments,
		&targets,

		&computeGradient,
		&getStatistics,
		
		&precision,
		&options.beta,
		&options.scale,
		&options.splitThreshold,
		&options.splitTargetThreshold,
		&options.levels,
		&options.directEvalThreshold,
		&options.periodicBoundaryConditions,
		&options.extrinsicCorrection,
		&options.useHilbertOrder,
		&options.directEvalAccuracy,
		&options.useFarfieldNearfieldThreads,

		&p,
		&options.pM,
		&options.pL,
		&options.s,
		&options.ws,
		&options.reducedScheme,
		&x,
		&w,
		&M,

		&papi
		)) goto _fail;

	if ((PyObject*)charges==Py_None) {
		Py_XDECREF(charges);
		charges = 0;
	}
	if ((PyObject*)dipoleMoments==Py_None) {
		Py_XDECREF(dipoleMoments);
		dipoleMoments = 0;
	}
	if ((PyObject*)targets==Py_None) {
		Py_XDECREF(targets);
		targets = 0;
	}
#ifdef USE_PAPI
	if ((PyObject*)papi==Py_None) {
		Py_XDECREF(papi);
		papi = 0;
	}
#endif	


	if(!PyArray_Check((PyObject*)particles)) {
		PyErr_Format(ErrorObject,"'particles' must be a numpy array.");
		goto _fail;
	}
	if (particles->nd != 2) {
		PyErr_Format(ErrorObject,"'particles' must be 2-dimensional.");
		goto _fail;
	}
	NParticles = particles->dimensions[0];
	
	if (targets) {
		if(!PyArray_Check((PyObject*)targets)) {
			PyErr_Format(ErrorObject,"'targets' must be a numpy array.");
			goto _fail;
		}
#if (FMM_DIM==2)
		if (particles->dimensions[1] != 2) {
			PyErr_Format(ErrorObject,"second dimension of 'particles' must be 2.");
			goto _fail;
		}
#else
		if (particles->dimensions[1] != 3) {
			PyErr_Format(ErrorObject,"second dimension of 'particles' must be 3.");
			goto _fail;
		}
#endif
		NTargets = targets->dimensions[0];
		targets1 = (PyArrayObject *)PyArray_ContiguousFromObject(
			(PyObject*)targets, PyArray_MYTYPE, 0, 0);
		if (targets1 == NULL) goto _fail;
#if (FMM_DIM==2)
		targets_data = (_FLOAT_ (*)[2])((PyArrayObject*)targets1)->data;
#else
		targets_data = (_FLOAT_ (*)[3])((PyArrayObject*)targets1)->data;
#endif
	}	
	else {
		NTargets = NParticles;
		targets_data = NULL;
	}
	
	if (charges) {
		if(!PyArray_Check((PyObject*)charges)) {
			PyErr_Format(ErrorObject,"'charges' must be a numpy array.");
			goto _fail;
		}
	       	if (charges->nd != 1) {
			PyErr_Format(ErrorObject,"'charges' must be 1-dimensional.");
			goto _fail;
		}
		if (charges->dimensions[0] != NParticles) {
			PyErr_Format(ErrorObject,"wrong length of 'charges'.");
			goto _fail;
		}	
	}

	if (dipoleMoments) {
		if(!PyArray_Check((PyObject*)dipoleMoments)) {
			PyErr_Format(ErrorObject,"'dipoleMoments' must be a numpy array.");
			goto _fail;
		}
		if (dipoleMoments->nd != 2) {
			PyErr_Format(ErrorObject,"'dipoleMoments' must be 2-dimensional.");
			goto _fail;
		}
		if (dipoleMoments->dimensions[0] != NParticles) {
			PyErr_Format(ErrorObject,"wrong length of 'dipoleMoments'.");
			goto _fail;
		}	
#if (FMM_DIM==2)
		if (dipoleMoments->dimensions[1] != 2) {
			PyErr_Format(ErrorObject,"second dimension of 'dipoleMoments' must be 2.");
			goto _fail;
		}	
#else
		if (dipoleMoments->dimensions[1] != 3) {
			PyErr_Format(ErrorObject,"second dimension of 'dipoleMoments' must be 3.");
			goto _fail;
		}	
#endif
	}	
	if (!charges && !dipoleMoments) {
		PyErr_Format(ErrorObject,"at least one of 'charges' or 'dipoleMoments' must be present.");
		goto _fail;
	}	

	particles1 = (PyArrayObject *)PyArray_ContiguousFromObject(
		(PyObject*)particles, PyArray_MYTYPE, 0, 0);
	if (particles1 == NULL) goto _fail;

	if (charges) {
		charges1 = (PyArrayObject *)PyArray_ContiguousFromObject(
			(PyObject*)charges, PyArray_MYTYPE, 0, 0);
		if (charges1 == NULL) goto _fail;	
		charges_data = (_FLOAT_ *)((PyArrayObject*)charges1)->data;
	}	

	if (dipoleMoments) {
		dipoleMoments1 = (PyArrayObject *)PyArray_ContiguousFromObject(
			(PyObject*)dipoleMoments, PyArray_MYTYPE, 0, 0);
		if (dipoleMoments1 == NULL) goto _fail;
#if (FMM_DIM==2)
		dipoleMoments_data = (_FLOAT_ (*)[2])((PyArrayObject*)dipoleMoments1)->data;
#else
		dipoleMoments_data = (_FLOAT_ (*)[3])((PyArrayObject*)dipoleMoments1)->data;
#endif
	}	
        
        dims[0] = NTargets;
	pot = (PyArrayObject *)PyArray_SimpleNew(1, dims, PyArray_MYTYPE);
	if (pot == NULL) goto _fail;

	grad_data = NULL;
	if (computeGradient) {
		dims[0] = NTargets;
#if (FMM_DIM==2)
		dims[1] = 2;
#else
		dims[1] = 3;
#endif
		grad = (PyArrayObject *)PyArray_SimpleNew(2, dims, PyArray_MYTYPE);
		if (grad == NULL) goto _fail;
#if (FMM_DIM==2)
		grad_data = (_FLOAT_ (*)[2])((PyArrayObject*)grad)->data;
#else
		grad_data = (_FLOAT_ (*)[3])((PyArrayObject*)grad)->data;
#endif
	}	
	
	if (getStatistics) {
		statistics_p = &statistics;
#ifdef USE_PAPI
		if (papi){
			PyObject *ev;
			char eventName[PAPI_MAX_STR_LEN];
			int eventCode;
			int i;

			ret =PyList_Check(papi);
			if (!ret) {
				PyErr_Format(ErrorObject,"'papi' must be an array of strings.");
				goto _fail;
			}	
			ret = PAPI_library_init(PAPI_VER_CURRENT);
			if (ret!=PAPI_VER_CURRENT) {
				fprintf(stderr, "PAPI error %d: %s\n",ret,PAPI_strerror(ret));
				PyErr_Format(ErrorObject,"PAPI_library_init failed.");
				goto _fail;
			}	
			
			options.PAPIeventSet = PAPI_NULL;	
			ret = PAPI_create_eventset(&options.PAPIeventSet);
			if (ret!=PAPI_OK) {
				fprintf(stderr, "PAPI error %d: %s\n",ret,PAPI_strerror(ret));
				PyErr_Format(ErrorObject,"PAPI_create_eventset failed.");
				goto _fail;
			}	
			NEvents=PyList_Size(papi);
			if (NEvents>MAX_NUM_PAPI_EVENTS) NEvents=MAX_NUM_PAPI_EVENTS;
			for (i=0; i<NEvents; i++) {
				ev = PyList_GetItem(papi, i);
				if (!PyString_Check(ev)) {
					PyErr_Format(ErrorObject,"'papi' must be an array of strings.");
					goto _fail;
				}
				strncpy(eventName, PyString_AsString(ev), PAPI_MAX_STR_LEN);
				ret =PAPI_event_name_to_code(eventName,&eventCode);
				if (ret!=PAPI_OK) {
					fprintf(stderr, "PAPI error %d: %s\n",ret,PAPI_strerror(ret));
					PyErr_Format(ErrorObject,"Unknown papi event name.");
					goto _fail;
				}	

				ret = PAPI_add_event(options.PAPIeventSet, eventCode);
				if (ret!=PAPI_OK) {
					fprintf(stderr, "PAPI error %d: %s\n",ret,PAPI_strerror(ret));
					PyErr_Format(ErrorObject,"PAPI_add_event failed.");
					goto _fail;
				}	
				events[i] = eventCode;
			}	
		}	
#endif	
	}
	else {
		statistics_p = NULL;
	}


	if (precision>=0) {
		switch(precision) {
		case 0:
			options.pM = 6;
			options.pL = 6;
			options.s = 8;
			options.ws = 1;
			options.directEvalAccuracy=0;
			break;
		case 1:
			options.pM = 16;
			options.pL = 16;
			options.s = 17;
			options.ws = 1;
			options.directEvalAccuracy=0;
			break;
		case 2:
			options.pM = 23;
			options.pL = 23;
			options.s = 26;
			options.ws = 1;
			options.directEvalAccuracy=2;
			break;
		default:
			PyErr_Format(ErrorObject,"'precision' must be one of 0, 1, 2.");
			goto _fail;
		}	
	}	
	else {
		if ((options.ws!=1)&&(options.ws!=2)) {
			PyErr_Format(ErrorObject,"'ws' must be 1 or 2.");
			goto _fail;
		}	
		if (p>=0) {
			options.pM = p;
			options.pL = p;
		}
		if (x||w||M) {
			PyObject *item;
			int i;

			if (!(x&&w&&M)) {
				PyErr_Format(ErrorObject,"If one of 'x', 'w', and 'M' is specifed, all of them must be specified.");
				goto _fail;
			}	
			if(!PyList_Check((PyObject*)x)) {
				PyErr_Format(ErrorObject,"'x' must be a list of numbers.");
				goto _fail;
			}
			if(!PyList_Check((PyObject*)w)) {
				PyErr_Format(ErrorObject,"'w' must be a list of numbers.");
				goto _fail;
			}
			if(!PyList_Check((PyObject*)M)) {
				PyErr_Format(ErrorObject,"'M' must be a list of numbers.");
				goto _fail;
			}
			options.s = PyList_Size(x);
			if ((PyList_Size(w)!=options.s)||PyList_Size(M)!=options.s) {
				PyErr_Format(ErrorObject,"Lengths of lists 'x', 'w', and 'M' must match.");
				goto _fail;
			}
			if (options.s > FMM_S_EPS_MAX) {
				PyErr_Format(ErrorObject,"Lists 'x', 'w', and 'M' too long.");
				goto _fail;
			}
			for (i=0; i<options.s; i++) {
				item = PyList_GetItem(x, i);
				if (!PyFloat_Check(item)) {
					PyErr_Format(ErrorObject,"'x' must be a list of numbers.");
					goto _fail;
				}
				user_defined_x[i] = PyFloat_AsDouble(item);
				
				item = PyList_GetItem(w, i);
				if (!PyFloat_Check(item)) {
					PyErr_Format(ErrorObject,"'w' must be a list of numbers.");
					goto _fail;
				}
				user_defined_w[i] = PyFloat_AsDouble(item);
				
				item = PyList_GetItem(M, i);
				if (!PyInt_Check(item)) {
					PyErr_Format(ErrorObject,"'M' must be a list of positive integers.");
					goto _fail;
				}
				user_defined_M[i] = PyInt_AsLong(item);
				if (user_defined_M[i]<=0) {
					PyErr_Format(ErrorObject,"'M' must be a list of positive integers.");
					goto _fail;
				}
			}
#if (FMM_DIM==3)
			options.x = user_defined_x;
			options.w = user_defined_w;
			options.M = user_defined_M;
#endif		
		}	
	}
	if ((options.directEvalAccuracy<0)||(options.directEvalAccuracy>2)) {
		PyErr_Format(ErrorObject,"'directEvalAccuracy' must be >= 0 and <= 2.");
		goto _fail;
	}
	
		
#if (FMM_DIM==2)
#ifdef USE_SINGLE_PRECISION
        fmmv2df
#else	
        fmmv2d
#endif		
#else
#ifdef USE_SINGLE_PRECISION
        fmmv3df
#else	
        fmmv3d
#endif		
#endif		
	    (	NParticles, 
#if (FMM_DIM==2)
		(_FLOAT_ (*)[2])((PyArrayObject*)particles1)->data,
#else
		(_FLOAT_ (*)[3])((PyArrayObject*)particles1)->data,
#endif
		charges_data,  
		dipoleMoments_data,
		NTargets,		
		targets_data,
		(_FLOAT_ *)((PyArrayObject*)pot)->data,
		grad_data,
	        &options,
		statistics_p,
                &err
	    );
	

	if (err) {
		PyErr_Format(ErrorObject,"%s", err);
		goto _fail;
	}	

	if (getStatistics) {
		stat = genStatistics(&statistics);
	}	
#ifdef USE_PAPI
	if (statistics.PAPIeventSet!=PAPI_NULL) {
		ret = PAPI_cleanup_eventset(statistics.PAPIeventSet);
		if (ret!=PAPI_OK) {
			fprintf(stderr, "PAPI error %d: %s\n",ret,PAPI_strerror(ret));
			PyErr_Format(ErrorObject,"PAPI_cleanup_eventset failed.");
			goto _fail;
		}	
	
		ret = PAPI_destroy_eventset(&statistics.PAPIeventSet);
		if (ret!=PAPI_OK) {
			fprintf(stderr, "PAPI error %d: %s\n",ret,PAPI_strerror(ret));
			PyErr_Format(ErrorObject,"PAPI_destroy_eventset failed.");
			goto _fail;
		}	
	}
#endif	
	

	Py_XDECREF(particles1);
	Py_XDECREF(charges1);
	Py_XDECREF(dipoleMoments1);
	Py_XDECREF(targets1);

	if (computeGradient) {
		if (getStatistics) {
			return Py_BuildValue("(NNN)", pot, grad, stat);
		}
		else {	
			return Py_BuildValue("(NN)", pot, grad); 
		}	
	}	
	else {
		if (getStatistics) {
			return Py_BuildValue("(NN)", pot, stat);
		}
		else {	
			return PyArray_Return(pot);
		}	
	}		
_fail:
	Py_XDECREF(particles1);
	Py_XDECREF(charges1);
	Py_XDECREF(dipoleMoments1);
	Py_XDECREF(targets1);
	return NULL;	
}



typedef struct {
    PyObject_HEAD
	    
    void *fmmvHandle;
    int NParticles;
    int NTargets;
    int dipoleSources;
    int computeGradient;
} fmmv_HandleObject;

static PyTypeObject fmmv_HandleType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "fmmv.Handle",             /*tp_name*/
    sizeof(fmmv_HandleObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    0,                         /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/
    "fmmv handle objects",           /* tp_doc */
};

static char py_fmmv_initialize__doc__[] ="";

static PyObject *
py_fmmv_initialize(PyObject *self, PyObject *args, PyObject *kwds)
{

	char *err = NULL;

	static char* argnames[] = {
		"particles", 
		"targets", 

		"dipoleSources",
		"computeGradient",
		"getStatistics",
		
		"precision",
		"beta",
		"scale",
		"splitThreshold",
		"splitTargetThreshold",
		"levels",
		"directEvalThreshold",
		"periodicBoundaryConditions",
		"extrinsicCorrection",
		"useHilbertOrder",
		"directEvalAccuracy",
		"useFarfieldNearfieldThreads",

		"p",
		"pM",
		"pL",
		"s",
		"ws",
		"reducedScheme",
		"x",
		"w",
		"M",
		
		"papi",
       		NULL };

	PyArrayObject *particles = NULL;
	PyArrayObject *targets = NULL;

	struct FmmvOptions options;
	struct FmmvStatistics statistics;
	struct FmmvStatistics *statistics_p;

	int typeSources = 0;
	int typeTargets = 0;
	int dipoleSources = 0;
	int computeGradient = 0;
	int getStatistics = 0;
	
	int NParticles;
	PyArrayObject *particles1 = NULL;
	
	int NTargets;
	PyArrayObject *targets1 = NULL;
#if (FMM_DIM==2)
	_FLOAT_ (*targets_data)[2] = NULL;
#else
	_FLOAT_ (*targets_data)[3] = NULL;
#endif
	
	PyObject *stat = NULL;
	fmmv_HandleObject *handle = NULL;

	void* fmmvHandle;

	PyObject *papi = 0;
	int precision = -1;
	int p = -1;
	
	PyObject *x = 0;
	PyObject *w = 0;
	PyObject *M = 0;
#ifdef USE_PAPI
	int ret;
	int events[MAX_NUM_PAPI_EVENTS];
	int NEvents=0;
#endif

	options = fmmvGetDefaultOptions();

	if (!PyArg_ParseTupleAndKeywords(
		args, kwds, "O|O" "iii" "iddiiiiiiiiii" "iiiiiiOOO" "O", argnames,
		&particles, 
		&targets,

		&dipoleSources,
		&computeGradient,
		&getStatistics,
		
		&precision,
		&options.beta,
		&options.scale,
		&options.splitThreshold,
		&options.splitTargetThreshold,
		&options.levels,
		&options.directEvalThreshold,
		&options.periodicBoundaryConditions,
		&options.extrinsicCorrection,
		&options.useHilbertOrder,
		&options.directEvalAccuracy,
		&options.useFarfieldNearfieldThreads,

		&p,
		&options.pM,
		&options.pL,
		&options.s,
		&options.ws,
		&options.reducedScheme,
		&x,
		&w,
		&M,

		&papi
		)) goto _fail;

	if ((PyObject*)targets==Py_None) {
		Py_XDECREF(targets);
		targets = 0;
	}
#ifdef USE_PAPI
	if ((PyObject*)papi==Py_None) {
		Py_XDECREF(papi);
		papi = 0;
	}
#endif	

	if(!PyArray_Check((PyObject*)particles)) {
		PyErr_Format(ErrorObject,"'particles' must be a numpy array.");
		goto _fail;
	}
	if (particles->nd != 2) {
		PyErr_Format(ErrorObject,"'particles' must be 2-dimensional.");
		goto _fail;
	}
	NParticles = particles->dimensions[0];
	particles1 = (PyArrayObject *)PyArray_ContiguousFromObject(
		(PyObject*)particles, PyArray_MYTYPE, 0, 0);
	if (particles1 == NULL) goto _fail;
	
	
	if (targets) {
		if(!PyArray_Check((PyObject*)targets)) {
			PyErr_Format(ErrorObject,"'targets' must be a numpy array.");
			goto _fail;
		}
#if (FMM_DIM==2)
		if (particles->dimensions[1] != 2) {
			PyErr_Format(ErrorObject,"second dimension of 'particles' must be 2.");
			goto _fail;
		}
#else
		if (particles->dimensions[1] != 3) {
			PyErr_Format(ErrorObject,"second dimension of 'particles' must be 3.");
			goto _fail;
		}
#endif
		NTargets = targets->dimensions[0];
		targets1 = (PyArrayObject *)PyArray_ContiguousFromObject(
			(PyObject*)targets, PyArray_MYTYPE, 0, 0);
		if (targets1 == NULL) goto _fail;
#if (FMM_DIM==2)
		targets_data = (_FLOAT_ (*)[2])((PyArrayObject*)targets1)->data;
#else
		targets_data = (_FLOAT_ (*)[3])((PyArrayObject*)targets1)->data;
#endif
	}	
	else {
		NTargets = NParticles;
		targets_data = NULL;
	}
	
	
	if (getStatistics) {
		statistics_p = &statistics;
#ifdef USE_PAPI
		if (papi){
			PyObject *ev;
			char eventName[PAPI_MAX_STR_LEN];
			int eventCode;
			int i;

			ret =PyList_Check(papi);
			if (!ret) {
				PyErr_Format(ErrorObject,"'papi' must be an array of strings.");
				goto _fail;
			}	
			ret = PAPI_library_init(PAPI_VER_CURRENT);
			if (ret!=PAPI_VER_CURRENT) {
				fprintf(stderr, "PAPI error %d: %s\n",ret,PAPI_strerror(ret));
				PyErr_Format(ErrorObject,"PAPI_library_init failed.");
				goto _fail;
			}	
			
			options.PAPIeventSet = PAPI_NULL;	
			ret = PAPI_create_eventset(&options.PAPIeventSet);
			if (ret!=PAPI_OK) {
				fprintf(stderr, "PAPI error %d: %s\n",ret,PAPI_strerror(ret));
				PyErr_Format(ErrorObject,"PAPI_create_eventset failed.");
				goto _fail;
			}	
			NEvents=PyList_Size(papi);
			if (NEvents>MAX_NUM_PAPI_EVENTS) NEvents=MAX_NUM_PAPI_EVENTS;
			for (i=0; i<NEvents; i++) {
				ev = PyList_GetItem(papi, i);
				if (!PyString_Check(ev)) {
					PyErr_Format(ErrorObject,"'papi' must be an array of strings.");
					goto _fail;
				}
				strncpy(eventName, PyString_AsString(ev), PAPI_MAX_STR_LEN);
				ret =PAPI_event_name_to_code(eventName,&eventCode);
				if (ret!=PAPI_OK) {
					fprintf(stderr, "PAPI error %d: %s\n",ret,PAPI_strerror(ret));
					PyErr_Format(ErrorObject,"Unknown papi event name.");
					goto _fail;
				}	

				ret = PAPI_add_event(options.PAPIeventSet, eventCode);
				if (ret!=PAPI_OK) {
					fprintf(stderr, "PAPI error %d: %s\n",ret,PAPI_strerror(ret));
					PyErr_Format(ErrorObject,"PAPI_add_event failed.");
					goto _fail;
				}	
				events[i] = eventCode;
			}	
		}	
#endif	
		
	}
	else {
		statistics_p = NULL;
	}	

	if (precision>=0) {
		switch(precision) {
		case 0:
			options.pM = 6;
			options.pL = 6;
			options.s = 8;
			options.ws = 1;
			options.directEvalAccuracy=0;
			break;
		case 1:
			options.pM = 16;
			options.pL = 16;
			options.s = 17;
			options.ws = 1;
			options.directEvalAccuracy=1;
			break;
		case 2:
			options.pM = 23;
			options.pL = 23;
			options.s = 26;
			options.ws = 1;
			options.directEvalAccuracy=2;
			break;
		default:
			PyErr_Format(ErrorObject,"'precision' must be one of 0, 1, 2.");
			goto _fail;
		}	
	}	
	else {
		if ((options.ws!=1)&&(options.ws!=2)) {
			PyErr_Format(ErrorObject,"'ws' must be 1 or 2.");
			goto _fail;
		}	
		if (p>=0) {
			options.pM = p;
			options.pL = p;
		}
		if (x||w||M) {
			PyObject *item;
			int i;

			if (!(x&&w&&M)) {
				PyErr_Format(ErrorObject,"If one of 'x', 'w', and 'M' is specifed, all of them must be specified.");
				goto _fail;
			}	
			if(!PyList_Check((PyObject*)x)) {
				PyErr_Format(ErrorObject,"'x' must be a list of numbers.");
				goto _fail;
			}
			if(!PyList_Check((PyObject*)w)) {
				PyErr_Format(ErrorObject,"'w' must be a list of numbers.");
				goto _fail;
			}
			if(!PyList_Check((PyObject*)M)) {
				PyErr_Format(ErrorObject,"'M' must be a list of numbers.");
				goto _fail;
			}
			options.s = PyList_Size(x);
			if ((PyList_Size(w)!=options.s)||PyList_Size(M)!=options.s) {
				PyErr_Format(ErrorObject,"Lengths of lists 'x', 'w', and 'M' must match.");
				goto _fail;
			}
			if (options.s > FMM_S_EPS_MAX) {
				PyErr_Format(ErrorObject,"Lists 'x', 'w', and 'M' too long.");
				goto _fail;
			}
			for (i=0; i<options.s; i++) {
				item = PyList_GetItem(x, i);
				if (!PyFloat_Check(item)) {
					PyErr_Format(ErrorObject,"'x' must be a list of numbers.");
					goto _fail;
				}
				user_defined_x[i] = PyFloat_AsDouble(item);
				
				item = PyList_GetItem(w, i);
				if (!PyFloat_Check(item)) {
					PyErr_Format(ErrorObject,"'w' must be a list of numbers.");
					goto _fail;
				}
				user_defined_w[i] = PyFloat_AsDouble(item);
				
				item = PyList_GetItem(M, i);
				if (!PyInt_Check(item)) {
					PyErr_Format(ErrorObject,"'M' must be a list of positive integers.");
					goto _fail;
				}
				user_defined_M[i] = PyInt_AsLong(item);
				if (user_defined_M[i]<=0) {
					PyErr_Format(ErrorObject,"'M' must be a list of positive integers.");
					goto _fail;
				}
			}
#if (FMM_DIM==3)
			options.x = user_defined_x;
			options.w = user_defined_w;
			options.M = user_defined_M;
#endif
		}	
	}

	if ((options.directEvalAccuracy<0)||(options.directEvalAccuracy>2)) {
		PyErr_Format(ErrorObject,"'directEvalAccuracy' must be >= 0 and <= 2.");
		goto _fail;
		
	}
         
        if (computeGradient) {
            typeTargets = POTENTIALS || GRADIENTS;
        }
        else {
            typeTargets = POTENTIALS;
        }    

        if (dipoleSources) {
            typeSources = CHARGES || DIPOLEMOMENTS;
        }
        else {
            typeSources = CHARGES;
        }    
	
		
#if (FMM_DIM==2)
#ifdef USE_SINGLE_PRECISION
        fmmv2df_initialize
#else	
        fmmv2d_initialize
#endif		
#else	
#ifdef USE_SINGLE_PRECISION
        fmmv3df_initialize
#else	
        fmmv3d_initialize
#endif		
#endif		
	    (
		&fmmvHandle,
		NParticles, 
#if (FMM_DIM==2)
		(_FLOAT_ (*)[2])((PyArrayObject*)particles1)->data,
#else
		(_FLOAT_ (*)[3])((PyArrayObject*)particles1)->data,
#endif
		typeSources,
		NTargets,		
		targets_data,
		typeTargets,
	        &options,
		statistics_p,
                &err 
	    );
	

	if (err) {
		PyErr_Format(ErrorObject,"%s", err);
		goto _fail;
	}	

	if (getStatistics) {
		stat = genStatistics(&statistics);
	}	

	Py_XDECREF(particles1);
	Py_XDECREF(targets1);

	handle = (fmmv_HandleObject *) PyObject_MALLOC(sizeof(fmmv_HandleObject));
	if (handle == NULL) return PyErr_NoMemory();
	PyObject_INIT(handle, &fmmv_HandleType);
	handle->fmmvHandle = fmmvHandle;
	handle->NParticles = NParticles;
	handle->NTargets = NTargets;
	handle->dipoleSources = dipoleSources;
	handle->computeGradient = computeGradient;

	if (getStatistics) {
		return Py_BuildValue("(NN)", handle, stat);
	}
	else {	
		return (PyObject*) handle; 
	}	

_fail:
	Py_XDECREF(particles1);
	Py_XDECREF(targets1);
	return NULL;	
}



static char py_fmmv_evaluate__doc__[] ="";

static PyObject *
py_fmmv_evaluate(PyObject *self, PyObject *args, PyObject *kwds)
{

	char *err = NULL;

	static char* argnames[] = {
		"handle", 
		"charges", 
		"dipoleMoments",
		"getStatistics",
       		NULL };

	fmmv_HandleObject *handle = NULL;

	PyArrayObject *charges = NULL;
	PyArrayObject *dipoleMoments = NULL;

	struct FmmvStatistics statistics;
	struct FmmvStatistics *statistics_p;

	PyArrayObject *charges1 = NULL;
	PyArrayObject *dipoleMoments1 = NULL;

	int getStatistics = 0;
	
	PyArrayObject *pot = NULL;
	PyArrayObject *grad = NULL;
	void *grad_data = NULL;
	_FLOAT_ *charges_data = NULL;
#if (FMM_DIM==2)
	_FLOAT_ (*dipoleMoments_data)[2] = NULL;
#else
	_FLOAT_ (*dipoleMoments_data)[3] = NULL;
#endif
	PyObject *stat = NULL;
	npy_intp dims[2];

	if (!PyArg_ParseTupleAndKeywords(
		args, kwds, "O|OO" "i" , argnames,
		&handle, 
		&charges,
		&dipoleMoments,
		&getStatistics
		)) goto _fail;


	// if (((handle->typeSources>=0)||(handle->typeSources>=2))&&(!charges)) {
	if (!charges) {
			PyErr_Format(ErrorObject,"no 'charges' specified.");
			goto _fail;
	}			
	
	if ((handle->dipoleSources)&&(!dipoleMoments)) {
			PyErr_Format(ErrorObject,"no 'dipoleMoments' specified.");
			goto _fail;
	}			
	
	if (charges) {
		if(!PyArray_Check((PyObject*)charges)) {
			PyErr_Format(ErrorObject,"'charges' must be a numpy array.");
			goto _fail;
		}
	       	if (charges->nd != 1) {
			PyErr_Format(ErrorObject,"'charges' must be 1-dimensional.");
			goto _fail;
		}
		if (charges->dimensions[0] != handle->NParticles) {
			PyErr_Format(ErrorObject,"wrong length of 'charges'.");
			goto _fail;
		}	
	}

	if (handle->dipoleSources) {
		if(!PyArray_Check((PyObject*)dipoleMoments)) {
			PyErr_Format(ErrorObject,"'dipoleMoments' must be a numpy array.");
			goto _fail;
		}
		if (dipoleMoments->nd != 2) {
			PyErr_Format(ErrorObject,"'dipoleMoments' must be 2-dimensional.");
			goto _fail;
		}
		if (dipoleMoments->dimensions[0] != handle->NParticles) {
			PyErr_Format(ErrorObject,"wrong length of 'dipoleMoments'.");
			goto _fail;
		}	
#if (FMM_DIM==2)
		if (dipoleMoments->dimensions[1] != 2) {
			PyErr_Format(ErrorObject,"second dimension of 'dipoleMoments' must be 2.");
			goto _fail;
		}	
#else
		if (dipoleMoments->dimensions[1] != 3) {
			PyErr_Format(ErrorObject,"second dimension of 'dipoleMoments' must be 3.");
			goto _fail;
		}	
#endif
	}	


	if (charges) {
		charges1 = (PyArrayObject *)PyArray_ContiguousFromObject(
			(PyObject*)charges, PyArray_MYTYPE, 0, 0);
		if (charges1 == NULL) goto _fail;	
		charges_data = (_FLOAT_ *)((PyArrayObject*)charges1)->data;
	}	

	if (handle->dipoleSources) {
		dipoleMoments1 = (PyArrayObject *)PyArray_ContiguousFromObject(
			(PyObject*)dipoleMoments, PyArray_MYTYPE, 0, 0);
		if (dipoleMoments1 == NULL) goto _fail;
#if (FMM_DIM==2)
		dipoleMoments_data = (_FLOAT_ (*)[2])((PyArrayObject*)dipoleMoments1)->data;
#else
		dipoleMoments_data = (_FLOAT_ (*)[3])((PyArrayObject*)dipoleMoments1)->data;
#endif
	}	
	
        dims[0] = handle->NTargets;
	pot = (PyArrayObject *)PyArray_SimpleNew(1, dims, PyArray_MYTYPE);
	if (pot == NULL) goto _fail;

	grad_data = NULL;
	if (handle->computeGradient) {
		dims[0] = handle->NTargets;
#if (FMM_DIM==2)
		dims[1] = 2;
#else
		dims[1] = 3;
#endif
		grad = (PyArrayObject *)PyArray_SimpleNew(2, dims, PyArray_MYTYPE);
		if (grad == NULL) goto _fail;
#if (FMM_DIM==2)
		grad_data = (_FLOAT_ (*)[2])((PyArrayObject*)grad)->data;
#else
		grad_data = (_FLOAT_ (*)[3])((PyArrayObject*)grad)->data;
#endif
	}	
	
	if (getStatistics) {
		statistics_p = &statistics;
	}
	else {
		statistics_p = NULL;
	}	
		
		
#if (FMM_DIM==2)
#ifdef USE_SINGLE_PRECISION
        fmmv2df_evaluate
#else
        fmmv2d_evaluate
#endif		
#else
#ifdef USE_SINGLE_PRECISION
        fmmv3df_evaluate
#else
        fmmv3d_evaluate
#endif		
#endif		
	    (
		handle->fmmvHandle, 
		charges_data,  
		dipoleMoments_data,
		(_FLOAT_ *)((PyArrayObject*)pot)->data,
		grad_data,
		statistics_p,
                &err
	    );
	

	if (err) {
		PyErr_Format(ErrorObject,"%s", err);
		goto _fail;
	}	

	if (getStatistics) {
		stat = genStatistics(&statistics);
	}	

	Py_XDECREF(charges1);
	Py_XDECREF(dipoleMoments1);

	if (handle->computeGradient) {
		if (getStatistics) {
			return Py_BuildValue("(NNN)", pot, grad, stat);
		}
		else {	
			return Py_BuildValue("(NN)", pot, grad); 
		}	
	}	
	else {
		if (getStatistics) {
			return Py_BuildValue("(NN)", pot, stat);
		}
		else {	
			return PyArray_Return(pot);
		}	
	}		
_fail:
	Py_XDECREF(charges1);
	Py_XDECREF(dipoleMoments1);
	return NULL;	
}


static char py_fmmv_finalize__doc__[] ="";

static PyObject *
py_fmmv_finalize(PyObject *self, PyObject *args, PyObject *kwds)
{

	char *err = NULL;

	static char* argnames[] = {
		"handle", 
		"getStatistics",
       		NULL };

	fmmv_HandleObject *handle = NULL;
	PyObject *stat = NULL;

	struct FmmvStatistics statistics;
	struct FmmvStatistics *statistics_p;

	int getStatistics = 0;
	
	if (!PyArg_ParseTupleAndKeywords(
		args, kwds, "O|i" , argnames,
		&handle, 
		&getStatistics
		)) goto _fail;

	if (getStatistics) {
		statistics_p = &statistics;
	}
	else {
		statistics_p = NULL;
	}	
		
#if (FMM_DIM==2)
#ifdef USE_SINGLE_PRECISION
        fmmv2df_finalize
#else	
        fmmv2d_finalize
#endif		
#else	
#ifdef USE_SINGLE_PRECISION
        fmmv3df_finalize
#else	
        fmmv3d_finalize
#endif		
#endif		
	    (
		handle->fmmvHandle, 
		statistics_p, 
                &err
	    );
	

	if (err) {
		PyErr_Format(ErrorObject,"%s", err);
		goto _fail;
	}	

	if (getStatistics) {
		stat = genStatistics(&statistics);
	}	

#ifdef USE_PAPI
	if (statistics.PAPIeventSet!=PAPI_NULL) {
		int ret;
		ret = PAPI_cleanup_eventset(statistics.PAPIeventSet);
		if (ret!=PAPI_OK) {
			fprintf(stderr, "PAPI error %d: %s\n",ret,PAPI_strerror(ret));
			PyErr_Format(ErrorObject,"PAPI_cleanup_eventset failed.");
			goto _fail;
		}	
	
		ret = PAPI_destroy_eventset(&statistics.PAPIeventSet);
		if (ret!=PAPI_OK) {
			fprintf(stderr, "PAPI error %d: %s\n",ret,PAPI_strerror(ret));
			PyErr_Format(ErrorObject,"PAPI_destroy_eventset failed.");
			goto _fail;
		}	
	}
#endif	

	if (getStatistics) {
		return stat;
	}
	else {	
		//Py_RETURN_NONE;
		return Py_INCREF(Py_None), Py_None;
	}		
_fail:
	return NULL;	
}






static char py_fmmv_direct__doc__[] ="";

static PyObject *
py_fmmv_direct(PyObject *self, PyObject *args, PyObject *kwds)
{

	char *err = NULL;

	static char* argnames[] = {
		"particles", 
		"charges", 
		"dipoleMoments",
		"targets", 
		"NTargets", 
		"computeGradient",
		"directEvalAccuracy",
		"extendedPrecision",
		"beta",
       		NULL };

	PyArrayObject *particles = NULL;
	PyArrayObject *charges = NULL;
	PyArrayObject *targets = NULL;
	PyArrayObject *dipoleMoments = NULL;

	double time;
	double beta;
	int computeGradient = 0;
#ifdef USE_SINGLE_PRECISION
	int directEvalAccuracy = 1;
#else        
	int directEvalAccuracy = 2;
#endif        
	int extendedPrecision = 0;
	
	int NParticles;
	PyArrayObject *particles1 = NULL;
	PyArrayObject *charges1 = NULL;
	PyArrayObject *dipoleMoments1 = NULL;
	
	int NTargets=-1234567;
	int NTargets0;
	PyArrayObject *targets1 = NULL;
	
	PyArrayObject *pot = NULL;
	PyArrayObject *grad = NULL;
	void *grad_data = NULL;
	_FLOAT_ *charges_data = NULL;
#if (FMM_DIM==2)
	_FLOAT_ (*targets_data)[2] = NULL;
	_FLOAT_ (*dipoleMoments_data)[2] = NULL;
#else
	_FLOAT_ (*targets_data)[3] = NULL;
	_FLOAT_ (*dipoleMoments_data)[3] = NULL;
#endif
	npy_intp dims[2];


	if (!PyArg_ParseTupleAndKeywords(
		args, kwds, "O|OOOiiiid", argnames,
		&particles, 
		&charges,
		&dipoleMoments,
		&targets,
		&NTargets,
		&computeGradient,
		&directEvalAccuracy,
		&extendedPrecision,
		&beta
		)) goto _fail;

	if ((PyObject*)charges==Py_None) {
		Py_XDECREF(charges);
		charges = 0;
	}
	if ((PyObject*)dipoleMoments==Py_None) {
		Py_XDECREF(dipoleMoments);
		dipoleMoments = 0;
	}
	if ((PyObject*)targets==Py_None) {
		Py_XDECREF(targets);
		targets = 0;
	}
#ifdef USE_PAPI
	if ((PyObject*)papi==Py_None) {
		Py_XDECREF(papi);
		papi = 0;
	}
#endif	


	if(!PyArray_Check((PyObject*)particles)) {
		PyErr_Format(ErrorObject,"'particles' must be a numpy array.");
		goto _fail;
	}
	if (particles->nd != 2) {
		PyErr_Format(ErrorObject,"'particles' must be 2-dimensional.");
		goto _fail;
	}
	NParticles = particles->dimensions[0];
	
	if (targets) {
		if(!PyArray_Check((PyObject*)targets)) {
			PyErr_Format(ErrorObject,"'targets' must be a numpy array.");
			goto _fail;
		}
#if (FMM_DIM==2)
		if (particles->dimensions[1] != 2) {
			PyErr_Format(ErrorObject,"second dimension of 'particles' must be 2.");
			goto _fail;
		}
#else
		if (particles->dimensions[1] != 3) {
			PyErr_Format(ErrorObject,"second dimension of 'particles' must be 3.");
			goto _fail;
		}
#endif
		NTargets0 = targets->dimensions[0];
		targets1 = (PyArrayObject *)PyArray_ContiguousFromObject(
			(PyObject*)targets, PyArray_MYTYPE, 0, 0);
		if (targets1 == NULL) goto _fail;
#if (FMM_DIM==2)
		targets_data = (_FLOAT_ (*)[2])((PyArrayObject*)targets1)->data;
#else
		targets_data = (_FLOAT_ (*)[3])((PyArrayObject*)targets1)->data;
#endif
	}	
	else {
		NTargets0 = NParticles;
		targets_data = NULL;
	}

	if (NTargets==-1234567) {
		NTargets = NTargets0;
	}
	else if (NTargets>NTargets0) {
		NTargets = NTargets0;
	}
	else if (NTargets<1) {	
		PyErr_Format(ErrorObject,"'NTargets' must be an integer >=1.");
		goto _fail;
	}

	
	if (charges) {
		if(!PyArray_Check((PyObject*)charges)) {
			PyErr_Format(ErrorObject,"'charges' must be a numpy array.");
			goto _fail;
		}
	       	if (charges->nd != 1) {
			PyErr_Format(ErrorObject,"'charges' must be 1-dimensional.");
			goto _fail;
		}
		if (charges->dimensions[0] != NParticles) {
			PyErr_Format(ErrorObject,"wrong length of 'charges'.");
			goto _fail;
		}	
	}

	if (dipoleMoments) {
		if(!PyArray_Check((PyObject*)dipoleMoments)) {
			PyErr_Format(ErrorObject,"'dipoleMoments' must be a numpy array.");
			goto _fail;
		}
		if (dipoleMoments->nd != 2) {
			PyErr_Format(ErrorObject,"'dipoleMoments' must be 2-dimensional.");
			goto _fail;
		}
		if (dipoleMoments->dimensions[0] != NParticles) {
			PyErr_Format(ErrorObject,"wrong length of 'dipoleMoments'.");
			goto _fail;
		}	
#if (FMM_DIM==2)
		if (dipoleMoments->dimensions[1] != 2) {
			PyErr_Format(ErrorObject,"second dimension of 'dipoleMoments' must be 2.");
			goto _fail;
		}	
#else
		if (dipoleMoments->dimensions[1] != 3) {
			PyErr_Format(ErrorObject,"second dimension of 'dipoleMoments' must be 3.");
			goto _fail;
		}	
#endif
	}	
	if (!charges && !dipoleMoments) {
		PyErr_Format(ErrorObject,"at least one of 'charges' or 'dipoleMoments' must be present.");
		goto _fail;
	}	
        if (!extendedPrecision) {
	    if ((directEvalAccuracy<0)||(directEvalAccuracy>2)) {
		PyErr_Format(ErrorObject,"'directEvalAccuracy' must be >= 0 and <= 2.");
		goto _fail;
	    }
        }    


	particles1 = (PyArrayObject *)PyArray_ContiguousFromObject(
		(PyObject*)particles, PyArray_MYTYPE, 0, 0);
	if (particles1 == NULL) goto _fail;

	if (charges) {
		charges1 = (PyArrayObject *)PyArray_ContiguousFromObject(
			(PyObject*)charges, PyArray_MYTYPE, 0, 0);
		if (charges1 == NULL) goto _fail;	
		charges_data = (_FLOAT_ *)((PyArrayObject*)charges1)->data;
	}	

	if (dipoleMoments) {
		dipoleMoments1 = (PyArrayObject *)PyArray_ContiguousFromObject(
			(PyObject*)dipoleMoments, PyArray_MYTYPE, 0, 0);
		if (dipoleMoments1 == NULL) goto _fail;
#if (FMM_DIM==2)
		dipoleMoments_data = (_FLOAT_ (*)[2])((PyArrayObject*)dipoleMoments1)->data;
#else
		dipoleMoments_data = (_FLOAT_ (*)[3])((PyArrayObject*)dipoleMoments1)->data;
#endif
	}	
	
        dims[0] = NTargets;
	pot = (PyArrayObject *)PyArray_SimpleNew(1, dims, PyArray_MYTYPE);
	if (pot == NULL) goto _fail;

	grad_data = NULL;
	if (computeGradient) {
		dims[0] = NTargets;
#if (FMM_DIM==2)
		dims[1] = 2;
#else
		dims[1] = 3;
#endif
		grad = (PyArrayObject *)PyArray_SimpleNew(2, dims, PyArray_MYTYPE);
		if (grad == NULL) goto _fail;
#if (FMM_DIM==2)
		grad_data = (_FLOAT_ (*)[2])((PyArrayObject*)grad)->data;
#else
		grad_data = (_FLOAT_ (*)[3])((PyArrayObject*)grad)->data;
#endif
	}	

	if (extendedPrecision) {
                directEvalAccuracy = -1;
        }

#if (FMM_DIM==2)
#ifdef USE_SINGLE_PRECISION
       	fmmv2df_direct
#else	
       	fmmv2d_direct
#endif		
#else	
#ifdef USE_SINGLE_PRECISION
       	fmmv3df_direct
#else	
       	fmmv3d_direct
#endif		
#endif		
    	(	NParticles, 
#if (FMM_DIM==2)
		(_FLOAT_ (*)[2])((PyArrayObject*)particles1)->data,
#else
		(_FLOAT_ (*)[3])((PyArrayObject*)particles1)->data,
#endif		
		charges_data,  
		dipoleMoments_data,
		NTargets,		
		targets_data,
		(_FLOAT_ *)((PyArrayObject*)pot)->data,
		grad_data,
                directEvalAccuracy, 
/* #if (FMM_DIM==3) */
		beta,
/* #endif */
        	&time,
                &err
    	);
	

	if (err) {
		PyErr_Format(ErrorObject,"%s", err);
		goto _fail;
	}	


	Py_XDECREF(particles1);
	Py_XDECREF(charges1);
	Py_XDECREF(dipoleMoments1);
	Py_XDECREF(targets1);

	if (computeGradient) {
		return Py_BuildValue("(NN)", pot, grad); 
	}	
	else {
		return PyArray_Return(pot);
	}		
_fail:
	Py_XDECREF(particles1);
	Py_XDECREF(charges1);
	Py_XDECREF(dipoleMoments1);
	Py_XDECREF(targets1);
	return NULL;	
}








/* List of methods defined in the module */
#if (FMM_DIM==2)
#ifdef USE_SINGLE_PRECISION
#define PREFIX "fmmv2df"
#else
#define PREFIX "fmmv2d"
#endif
#else
#ifdef USE_SINGLE_PRECISION
#define PREFIX "fmmv3df"
#else
#define PREFIX "fmmv3d"
#endif
#endif

static struct PyMethodDef fmmv_methods[] = {
	{PREFIX , (PyCFunction)py_fmmv, METH_VARARGS | METH_KEYWORDS, py_fmmv__doc__},
	{PREFIX "_initialize", (PyCFunction)py_fmmv_initialize, METH_VARARGS | METH_KEYWORDS, py_fmmv_initialize__doc__},
	{PREFIX "_evaluate", (PyCFunction)py_fmmv_evaluate, METH_VARARGS | METH_KEYWORDS, py_fmmv_evaluate__doc__},
	{PREFIX "_finalize", (PyCFunction)py_fmmv_finalize, METH_VARARGS | METH_KEYWORDS, py_fmmv_finalize__doc__},
	{PREFIX "_direct", (PyCFunction)py_fmmv_direct, METH_VARARGS | METH_KEYWORDS, py_fmmv_direct__doc__},
	{NULL,		NULL}		/* sentinel */
};


/* Initialization function for the module (*must* be called initfmmv) */

static char fmmv_module_documentation[] = 
""
;

PyMODINIT_FUNC
#if (FMM_DIM==2)
#ifdef EXTRA_ICC_VERSION
#ifdef USE_SINGLE_PRECISION
initfmmv2df_icc(void)
#else
initfmmv2d_icc(void)
#endif	
#else	
#ifdef USE_SINGLE_PRECISION
initfmmv2df(void)
#else
initfmmv2d(void)
#endif	
#endif	
#else
#ifdef EXTRA_ICC_VERSION
#ifdef USE_SINGLE_PRECISION
initfmmv3df_icc(void)
#else
initfmmv3d_icc(void)
#endif	
#else	
#ifdef USE_SINGLE_PRECISION
initfmmv3df(void)
#else
initfmmv3d(void)
#endif	
#endif	

#endif	
{
	PyObject *m, *d;

	fmmv_HandleType.tp_new = PyType_GenericNew;
	if (PyType_Ready(&fmmv_HandleType) < 0)
        return;
	
	/* Create the module and add the functions */
	m = Py_InitModule4(
#if (FMM_DIM==2)
#ifdef EXTRA_ICC_VERSION
#ifdef USE_SINGLE_PRECISION
			"fmmv2df_icc", 
#else			
			"fmmv2d_icc", 
#endif			
#else			
#ifdef USE_SINGLE_PRECISION
			"fmmv2df", 
#else			
			"fmmv2d", 
#endif			
#endif			
#else			
#ifdef EXTRA_ICC_VERSION
#ifdef USE_SINGLE_PRECISION
			"fmmv3df_icc", 
#else			
			"fmmv3d_icc", 
#endif			
#else			
#ifdef USE_SINGLE_PRECISION
			"fmmv3df", 
#else			
			"fmmv3d", 
#endif			
#endif			
#endif			
			fmmv_methods,
			fmmv_module_documentation,
			(PyObject*)NULL,PYTHON_API_VERSION);
	
	/* Import the array object */
	import_array();
	
	/* Add some symbolic constants to the module */
	d = PyModule_GetDict(m);
#if (FMM_DIM==2)
#ifdef USE_SINGLE_PRECISION
	ErrorObject = PyString_FromString("fmmv2df.error");
#else	
	ErrorObject = PyString_FromString("fmmv2d.error");
#endif	
#else	
#ifdef USE_SINGLE_PRECISION
	ErrorObject = PyString_FromString("fmmv3df.error");
#else	
	ErrorObject = PyString_FromString("fmmv3d.error");
#endif	
#endif	
	PyDict_SetItemString(d, "error", ErrorObject);
	
	Py_INCREF(&fmmv_HandleType);
	PyModule_AddObject(m, "Handle", (PyObject *)&fmmv_HandleType);	
	
	/* Check for errors */
	if (PyErr_Occurred())
#if (FMM_DIM==2)
#ifdef USE_SINGLE_PRECISION
		Py_FatalError("can't initialize module fmmv2df");
#else	
		Py_FatalError("can't initialize module fmmv2d");
#endif		
#else	
#ifdef USE_SINGLE_PRECISION
		Py_FatalError("can't initialize module fmmv3df");
#else	
		Py_FatalError("can't initialize module fmmv3d");
#endif		
#endif	
  
}

