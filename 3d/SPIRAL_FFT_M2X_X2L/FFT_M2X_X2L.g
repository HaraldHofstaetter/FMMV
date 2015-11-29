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

for M in [2, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64]
do 

	F_M2X := SPLNonTerminal("DFT",[M,+1]); # +1 !!!
	Implement(F_M2X, rec(
			  search:="DP", 
			  file:=ConcatenationString("FFT_M2X_",
						    StringInt(M))
		 ));

	F_X2L := SPLNonTerminal("DFT",[M,-1]); # -1 !!!
	
	Implement(F_X2L, rec(search:="DP", 
			  file:=ConcatenationString("FFT_X2L_",
						    StringInt(M))
		 ));
od;		 


	 
