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

include ../Makefile.inc

all: setup.py fmmvmodule.c
	NUMPY_INCLUDE=$(NUMPY_INCLUDE) $(PYTHON) setup.py build --force   --debug 
	$(PYTHON) setup.py install  --install-platlib=.
	LD_LIBRARY_PATH=../lib:$(LD_LIBRARY_PATH) $(PYTHON) fmmv_demo.py &

#	python fmmv_demo.py &

#	LD_LIBRARY_PATH=../lib:$(LD_LIBRARY_PATH) python fmmv_demo.py &

# Delete the "--install-platlib=." option above to install 
# the FMMV Python modules into the standard location for 
# third-party Python modules. You will need to have root 
# permissions for this.

clean:
	rm -f *.so
	rm -r -f build
