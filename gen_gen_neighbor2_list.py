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

N2   = {
	 0.0 : "0", 
	 1.0: "1",
	-1.0 : "m1",
	 2.0: "2",
	-2.0 : "m2"}

def gen_neighbor2_list_2d():
	childs = ['SW', 'NW', 'SE', 'NE',]
	neighbors = [ 'SW', 'MW', 'NW', 'SM', 'NM', 'SE', 'ME', 'NE'] 

	Dir = { 'SW': (-1, -1),
        'MW': (-1,  0),
        'NW': (-1, +1),
        'SM': ( 0, -1),
        'MM': ( 0,  0),
        'NM': ( 0, +1),
        'SE': (+1, -1),
        'ME': (+1,  0),
        'NE': (+1, +1) }
	print "enum {"
	i = 0
	n2_num = 0;
	#direct neighbors first
	for neighbor in neighbors:
		(dx, dy) = Dir[neighbor]
		print '\tN2_%s_%s,' % (N2[dx], N2[dy]),
		n2_num += 1
		i += 1
		if i==4:
			print
			i=0
	for dx in [-2, -1, 0, 1, 2]:		
		for dy in [-2, -1, 0, 1, 2]:	
			if not(abs(dx)<=1 and abs(dy)<=1):
				print '\tN2_%s_%s,' % (N2[dx], N2[dy]),
				n2_num += 1
				i += 1
				if i==4:
					print
					i=0
	print "};"
	print
	print "void gen_neighbor2_list(Box *box, Box **list)"
	print "{"
	print "\tBox *parent = box->parent;"
	print "\tBox *nb;"
	print
	print "\tmemcpy(list, box->neighbor, %i*sizeof(Box*));" % (8);
	print "\tmemset(list+8, 0, %i*sizeof(Box*));" % (n2_num-8);
	print
	for child in childs:
		n = 0
		dir_child = Dir[child]
                print "\tif (box->whichChild==%s) {" % child
		for neighbor in neighbors:
			dir_neighbor = Dir[neighbor]
			print "\t\tnb = parent->neighbor[%s_];" % neighbor
			print "\t\tif (nb) {"
			for child1 in childs:
				dir_child1 = Dir[child1]
				x = -0.5*dir_child[0] + 2*dir_neighbor[0] + 0.5*dir_child1[0]
				y = -0.5*dir_child[1] + 2*dir_neighbor[1] + 0.5*dir_child1[1]
				if ((abs(x)<=2 and abs(y)<=2)
				 and not (abs(x)<=1 and abs(y)<=1)):
					print '\t\t\tlist[N2_%s_%s] = nb->child[%s];' % ( N2[x], N2[y], child1)
					n += 1		
			print "\t\t}"
		print "\t}  /* %i */" % n
	print "}"	


def gen_neighbor2_list_3d():
	childs = ['SWD', 'NWD', 'SED', 'NED', 'SWU', 'NWU', 'SEU', 'NEU']
	neighbors = [ 'SWD', 'MWD', 'NWD', 'SMD', 'MMD', 'NMD', 'SED', 'MED', 'NED',
		      'SWM', 'MWM', 'NWM', 'SMM', 'NMM', 'SEM', 'MEM', 'NEM', 'SWU',
		      'MWU', 'NWU', 'SMU', 'MMU', 'NMU', 'SEU', 'MEU', 'NEU']

	Dir = { 'SWD': (-1, -1, -1),
        'MWD': (-1,  0, -1),
        'NWD': (-1, +1, -1),
        'SMD': ( 0, -1, -1),
        'MMD': ( 0,  0, -1),
        'NMD': ( 0, +1, -1),
        'SED': (+1, -1, -1),
        'MED': (+1,  0, -1),
        'NED': (+1, +1, -1),
	'SWM': (-1, -1,  0),
        'MWM': (-1,  0,  0),
        'NWM': (-1, +1,  0),
        'SMM': ( 0, -1,  0),
        'NMM': ( 0, +1,  0),
        'SEM': (+1, -1,  0),
        'MEM': (+1,  0,  0),
        'NEM': (+1, +1,  0),
        'SWU': (-1, -1, +1),
        'MWU': (-1,  0, +1),
        'NWU': (-1, +1, +1),
        'SMU': ( 0, -1, +1),
        'MMU': ( 0,  0, +1),
        'NMU': ( 0, +1, +1),
        'SEU': (+1, -1, +1),
        'MEU': (+1,  0, +1),
        'NEU': (+1, +1, +1) }

	print "enum {"
	i = 0
	n2_num = 0;
	#direct neighbors first
	for neighbor in neighbors:
		(dx, dy, dz) = Dir[neighbor]
		print '\tN2_%s_%s_%s,' % (N2[dx], N2[dy], N2[dz]),
		n2_num += 1
		i += 1
		if i==4:
			print
			i=0
	for dx in [-2, -1, 0, 1, 2]:		
		for dy in [-2, -1, 0, 1, 2]:	
			for dz in [-2, -1, 0, 1, 2]:
				if not(abs(dx)<=1 and abs(dy)<=1 and abs(dz)<=1):
					print '\tN2_%s_%s_%s,' % (N2[dx], N2[dy], N2[dz]),
					n2_num += 1
					i += 1
					if i==4:
						print
						i=0
	print "};"
	print
	print "void gen_neighbor2_list(Box *box, Box **list)"
	print "{"
	print "\tBox *parent = box->parent;"
	print "\tBox *nb;"
	print
	print "\tmemcpy(list, box->neighbor, %i*sizeof(Box*));" % (26);
	print "\tmemset(list+26, 0, %i*sizeof(Box*));" % (n2_num-26);
	print
	for child in childs:
		n = 0
		dir_child = Dir[child]
                print "\tif (box->whichChild==%s) {" % child
		for neighbor in neighbors:
			dir_neighbor = Dir[neighbor]
			print "\t\tnb = parent->neighbor[%s_];" % neighbor
			print "\t\tif (nb) {"
			for child1 in childs:
				dir_child1 = Dir[child1]
				x = -0.5*dir_child[0] + 2*dir_neighbor[0] + 0.5*dir_child1[0]
				y = -0.5*dir_child[1] + 2*dir_neighbor[1] + 0.5*dir_child1[1]
				z = -0.5*dir_child[2] + 2*dir_neighbor[2] + 0.5*dir_child1[2]
				if ((abs(x)<=2 and abs(y)<=2 and abs(z)<=2)
				 and not (abs(x)<=1 and abs(y)<=1 and abs(z)<=1)):
					print '\t\t\tlist[N2_%s_%s_%s] = nb->child[%s];' % ( N2[x], N2[y], N2[z], child1)
					n += 1		
			print "\t\t}"
		print "\t}  /* %i */" % n
	print "}"	



print '''	
/* This file is automatically generated by gen_gen_neighbor_list2.py */
/* DO NOT EDIT! */

#include "_fmmv.h"
#include <string.h> /* memset */
'''
print "#if (FMM_DIM==2)"
print
gen_neighbor2_list_2d()
print
print "#elif (FMM_DIM==3)"
print
gen_neighbor2_list_3d()
print
print "#endif /* FMM_DIM==3 */"
