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


from os import environ
from numpy.distutils.core import setup
from numpy.distutils.extension import Extension

NUMPY_INCLUDE = environ.get('NUMPY_INCLUDE')
#if NUMPY_INCLUDE==None:
#    NUMPY_INCLUDE = "/usr/lib/python2.7/dist-packages/numpy/core/include/numpy"
#    #NUMPY_INCLUDE = '/usr/local/lib/python2.6/site-packages/numpy/core/include/numpy'
#    #NUMPY_INCLUDE = '/usr/local/lib/python2.6/dist-packages/numpy/core/include/numpy'
#    #NUMPY_INCLUDE = '/System/Library/Frameworks/Python.framework/Versions/2.6/Extras/lib/python/numpy/core/include/numpy'


from os.path import expanduser
HOME = expanduser('~')

module2d = Extension('fmmv2d', 
	define_macros=[
		('FMM_DIM', '2')
	],
	include_dirs = ['../include', NUMPY_INCLUDE],
	libraries = ['fmmv2d'],
	library_dirs =['../lib'],
	sources = ['fmmvmodule.c'])

module2df = Extension('fmmv2df', 
	define_macros=[('USE_SINGLE_PRECISION', '1'), 
		('FMM_DIM', '2')
	],
	include_dirs = ['../include', NUMPY_INCLUDE],
	libraries = ['fmmv2df'],
	library_dirs =['../lib'],
	sources = ['fmmvmodule.c'])

module3d = Extension('fmmv3d', 
	define_macros=[
		('FMM_DIM', '3')
	],
	include_dirs = ['../include', NUMPY_INCLUDE],
	libraries = ['fmmv3d'],
	library_dirs =['../lib'],
	sources = ['fmmvmodule.c'])

module3df = Extension('fmmv3df', 
	define_macros=[('USE_SINGLE_PRECISION', '1'), 
		('FMM_DIM', '3')
	],
	include_dirs = ['../include', NUMPY_INCLUDE],
	libraries = ['fmmv3df'],
	library_dirs =['../lib'],
	sources = ['fmmvmodule.c'])


setup (name = 'fmmv',
	version = '2.0',
	description = 'The Fastest Multipole Method of Vienna',
	author = 'Harald Hofstaetter',
	author_email = 'hofi@harald-hofstaetter.at',
        url = "http://www.harald-hofstaetter.at",
	ext_modules = [module2d, module2df, module3d, module3df])
