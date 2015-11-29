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


childs = ['SW', 'NW', 'SE', 'NE']
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
	
InvDir = dict([(val, key) for key, val in Dir.items()])
	
N   = {    0 : "0",
	  .5 : "1",
	 -.5 : "m1",
	 1.0 : "2",
	-1.0 : "m2",
	 1.5 : "3",
	-1.5 : "m3",
	 2.0 : "4",
	-2.0 : "m4",
	 2.5 : "5",
	-2.5 : "m5",
	 3.0 : "6",
	-3.0 : "m6",	
	 3.5 : "7",
	-3.5 : "m7",	
	 4.0 : "8",
	-4.0 : "m8",	
	 4.5 : "9",
	-4.5 : "m9"}	


def gen_M2L_enums():
	print '''
enum {XU_all, XU_D, XD_all, XD_D,
      XU_NE, XU_SE, XU_NW, XU_SW,
      XD_NE, XD_SE, XD_NW, XD_SW,
      Xh};
	'''      
        print "enum {"
        for neighbor in neighbors:
                dir_neighbor = Dir[neighbor]
                for child in childs:
                        dir_child = Dir[child]
                        x = 2*dir_neighbor[0] + 0.5*dir_child[0]
                        y = 2*dir_neighbor[1] + 0.5*dir_child[1]
                        print '\ti_%s_%s,' % (N[x], N[y]),
                print
        print "};"
	print
        n = 0;
        print "enum {"
        for x in [-0.5, 0.5]:
            for y in [0.5]:
                for z in [-0.5, +0.5]:
                    print '\td_%s_%s,' % (N[x], N[y]),
		    n += 1
        print	
        for x in [-1.5, -0.5, 0.5, 1.5]:
            for y in [0.5, 1.5]:
                for z in [1.5]:
                    print '\td_%s_%s,' % (N[x], N[y]),
		    n += 1
            print
        for x in [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]:
            for y in [0.5, 1.5, 2.5]:
                for z in [2.5]:
                    print '\td_%s_%s,' % (N[x], N[y]),
		    n += 1
	    print	
        print "\tN_D_X2X /* = %i */" % n
        print "};"
	print

        print "enum {"
        for neighbor in neighbors:
                dir_neighbor = Dir[neighbor]
                for child in childs:
                        dir_child = Dir[child]
                        x = 2*dir_neighbor[0] + 0.5*dir_child[0]
                        y = 2*dir_neighbor[1] + 0.5*dir_child[1]
                        print '\ti2_%s_%s,' % (N[x], N[y]),
                print
        for xx in [-2, -1, 0, 1, 2]:
             for yy in [-2, -1, 0, 1, 2]:
                        if not(abs(xx)<=1 and abs(yy)<=1 and abs(zz)<=1):
                            for child in childs:
                                dir_child = Dir[child]
                                x = 2*xx + 0.5*dir_child[0]
                                y = 2*yy + 0.5*dir_child[1]
                                print '\ti2_%s_%s,' % (N[x], N[y]),
                            print
        print "};"
        print
        n = 0;
        print "enum {"
        for x in [-0.5, +0.5]:
            for y in [0.5]:
                for z in [-0.5, +0.5]:
                    print '\td2_%s_%s,' % (N[x], N[y]),
		    n += 1
        print	
        for x in [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]:
            for y in [0.5, 1.5, 2.5]:
                for z in [2.5]:
                    print '\td2_%s_%s,' % (N[x], N[y]),
		    n += 1
	    print	
        for x in [-4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5]:
            for y in [0.5, 1.5, 2.5, 3.5, 4.5]:
                for z in [3.5, 4.5]:
                    print '\td2_%s_%s,' % (N[x], N[y]),
		    n += 1
	    print	
        print "\tN_D_X2X_ws2 /* = %i */" % n
        print "};"
        print


        print
        i = 0;
        print "enum {"
        for x in [-2, -1, 0, 1, 2]:
          for y in [-2, -1, 0, 1, 2]:
                if not(abs(x)<=1 and abs(y)<=1 and abs(z)<=1):
                        print '\tir_%s_%s,' % (N[2*x], N[2*y]),
                        i += 1
                        if i==4:
                                print
                                i = 0
        print "\n};"
        print
        n = 0;
        print "enum {"
        for x in [-0.5, +0.5]:
            for y in [0.5]:
                for z in [-0.5, +0.5]:
                    print '\tdr_%s_%s,' % (N[x], N[y]),
		    n += 1
        print	
        for x in [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]:
            for y in [0.5, 1.5, 2.5]:
                for z in [2.5]:
                    print '\tdr_%s_%s,' % (N[x], N[y]),
		    n += 1
	    print	
        for x in [-4, -2, 0, 2, 4]:
            for y in [0, 2, 4]:
                for z in [4]:
                    print '\tdr_%s_%s,' % (N[x], N[y]),
		    n += 1
	    print	
        print "\tN_D_X2X_ws2_reduced /* = %i */" % n
        print "\n};"
        print


def Id(x, y, z):
	return (x, y, z)
	
def Ry_pi2(x, y, z):
	return (z, y, -x) 

def Ry_pi2_Rz_pi_2(x, y, z):	
	return (y, z, x) 


def gen_M2L_tables0():
	for (dir_no, dir1, dir2, J1) in [(0, "U", "D", Id), (1, "N", "S", Ry_pi2_Rz_pi_2), (2, "E", "W", Ry_pi2) ]:
		def J2(x, y, z):  # apply Ry_pi
			return J1(-x, y, -z)

		print "static int T0_%s%s[][9] = {" % (dir1, dir2)
		n = 0
                for child in childs:
                        (x, y, z) = Dir[child]

                        (tx1, ty1, tz1) = (-0.5*x, -0.5*y, -0.5*z)
                        (tx2, ty2, tz2) = (0.5*x, -0.5*y, 0.5*z)
			conj1 = 0
			conj2 = 0
			if ty1<0:
				tx1 = -tx1
				ty1 = -ty1
				conj1 = 1
			if ty2<0:
				tx2 = -tx2
				ty2 = -ty2
				conj2 = 1
			n += 1	
                        if z>0:
				print "\t{%s, d_%s_%s_%s, %i, Xh, XU_all," % (InvDir[J1(x,y,z)], N[tx1], N[ty1], N[tz1], conj1),
                        else:
                                if dir1 == "U":
					print "\t{%s, d_%s_%s_%s, %i, Xh, XU_D," % (InvDir[J1(x,y,z)], N[tx1], N[ty1], N[tz1], conj1),
                                else:
					print "\t{%s, d_%s_%s_%s, %i, XU_%s, XU_D," % (InvDir[J1(x,y,z)], N[tx1], N[ty1], N[tz1], conj1, InvDir[(x,y,z)][:2]),
                        if z<0:
				print "d_%s_%s_%s, %i, Xh, XD_all}," % (N[tx2], N[ty2], N[tz2], conj2)
                        else:
                                if dir1 == "U":
					print "d_%s_%s_%s, %i, Xh, XD_D}," % (N[tx2], N[ty2], N[tz2], conj2)
                                else:
					print "d_%s_%s_%s, %i, XD_%s, XD_D}," % (N[tx2], N[ty2], N[tz2], conj2, InvDir[(-x,y,-z)][:2])
		print "\t{-1, -1, -1, -1, -1} /* # = %i */" % n
		print "};"

		
def gen_M2L_tables_reduced():
	#alreadyHandled = set([])	
	alreadyHandled = {}	
	for (dir_no, dir1, dir2, J1) in [(0, "U", "D", Id), (1, "N", "S", Ry_pi2_Rz_pi_2), (2, "E", "W", Ry_pi2) ]:
		def J2(x, y, z):  # apply Ry_pi
			return J1(-x, y, -z)
		
		for (dir, tdir, J)  in [("U", dir1, J1), ("D", dir2, J2)]:
			print "static int T1_%s_ws2_reduced[][3] = {" % tdir
			n = 0
			print "\t/* upper (reduced) */"
			coords = []
			for x in [-4, -2, 0, 2, 4]:
				for y in [0, 2, 4]:
					for z in [4]:
						coords.append((x,y,z))
			threshold = 4 
			for (dx, dy, dz) in coords:
				(tx, ty, tz) = J(dx, dy, dz)
				if not(dir1=="N" and abs(tz)>=threshold) and  not(dir1=="E" and (abs(tz)>=threshold or abs(ty)>=threshold)):
					#Note: case  dx==0 and dy==0 needs extra handling !!!
					if dy!=0 or dx>0:
						mdx = -dx
			 			mdy = -dy
						mdz = dz
						(mtx, mty, mtz) = J(mdx, mdy, mdz)
						print "\t{ir_%s_%s_%s, ir_%s_%s_%s, dr_%s_%s_%s}," % (N[tx], N[ty], N[tz], 
												  N[mtx], N[mty], N[mtz], 
												  N[dx], N[dy], N[dz])
						n += 1												  
			print "\t{-1, -1, -1} /* # = %i */"% n
			print "};"

def gen_M2L_tables1(ws=1):
	if ws==2:
		suf ="2"
	else:
		suf =""
	#alreadyHandled = set([])	
	alreadyHandled = {}
	for (dir_no, dir1, dir2, J1) in [(0, "U", "D", Id), (1, "N", "S", Ry_pi2_Rz_pi_2), (2, "E", "W", Ry_pi2) ]:
		def J2(x, y, z):  # apply Ry_pi
			return J1(-x, y, -z)
		
		for (dir, tdir, J)  in [("U", dir1, J1), ("D", dir2, J2)]:
			if ws==2:
				print "static int T1_%s_ws2[][3] = {" % tdir
			else:	
				print "static int T1_%s[][3] = {" % tdir
			n = 0
			print "\t/* lower A */"
			if ws==2:
				coords = []
				for x in [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]:
					for y in [0.5, 1.5, 2.5]:
						for z in [2.5]:
							coords.append((x,y,z))
				threshold = 2.5
			else:
				coords = []
				for x in [-1.5, -0.5, 0.5, 1.5]:
					for y in [0.5, 1.5]:
						for z in [1.5]:
							coords.append((x,y,z))
				threshold = 1.5
			for (dx, dy, dz) in coords:
				(tx, ty, tz) = J(dx, dy, dz)
				if not(dir1=="N" and abs(tz)>=threshold) and  not(dir1=="E" and (abs(tz)>=threshold or abs(ty)>=threshold)):
					mdx = -dx
				 	mdy = -dy
					mdz = dz 
					(mtx, mty, mtz) = J(mdx, mdy, mdz)
					for (sx,sy,sz) in [J(-1,-1,-1), J(-1,1,-1), J(1,-1,-1), J(1,1,-1)]:
						#alreadyHandled.add((sx, sy, sz, tx, ty, tz))
						#alreadyHandled.add((sx, sy, sz, mtx, mty, mtz))
						alreadyHandled[(sx, sy, sz, tx, ty, tz)] = True
						alreadyHandled[(sx, sy, sz, mtx, mty, mtz)] = True
	
					print "\t{i%s_%s_%s_%s, i%s_%s_%s_%s, d%s_%s_%s_%s}," % (suf, N[tx], N[ty], N[tz], 
												 suf, N[mtx], N[mty], N[mtz], 
												 suf, N[dx], N[dy], N[dz])
					n += 1
			print "\t/* # = %i */" % n		
			print "\t/* upper */"
			if ws==2:
				coords = []
				for x in [-4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5]:
					for y in [0.5, 1.5, 2.5, 3.5, 4.5]:
						for z in [3.5, 4.5]:
							coords.append((x,y,z))
				threshold = 3.5
			  	
			else:
				coords = []
				for x in [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]:
					for y in [0.5, 1.5, 2.5]:
						for z in [2.5]:
							coords.append((x,y,z))
				threshold = 2.5
			for (dx, dy, dz) in coords:
				(tx, ty, tz) = J(dx, dy, dz)
				if not(dir1=="N" and abs(tz)>=threshold) and  not(dir1=="E" and (abs(tz)>=threshold or abs(ty)>=threshold)):
					mdx = -dx
			 		mdy = -dy
					mdz = dz
					(mtx, mty, mtz) = J(mdx, mdy, mdz)
					print "\t{i%s_%s_%s_%s, i%s_%s_%s_%s, d%s_%s_%s_%s}," % (suf, N[tx], N[ty], N[tz], 
												 suf, N[mtx], N[mty], N[mtz], 
												 suf, N[dx], N[dy], N[dz])
					n += 1							 
			print "\t{-1, -1, -1} /* # = %i */"% n
			print "};"
	return alreadyHandled

			
def gen_M2L_tables2(alreadyHandled, ws=1):
	if ws==2:
		suf ="2"
	else:
		suf =""
	#for (dir_no, dir1, dir2, J1) in [(0, "U", "D", Id), (1, "N", "S", Ry_pi2_Rz_pi_2), (2, "E", "W", Ry_pi2) ]: # Note: T2_U, T2_D, T2_U_ws2, T2_D_ws2 empty
	for (dir_no, dir1, dir2, J1) in [(1, "N", "S", Ry_pi2_Rz_pi_2), (2, "E", "W", Ry_pi2) ]:
		def J2(x, y, z):  # apply Ry_pi
			return J1(-x, y, -z)
			
		for (dir, tdir, J)  in [("U", dir1, J1), ("D", dir2, J2)]:
			print "\t/* lower B */"
			if ws==2:
				print "static int T2_%s_ws2[][4] = {" % tdir
			else:
				print "static int T2_%s[][4] = {" % tdir
			n = 0
			if ws==2 or ws=="2r":
				coords = [-2.5, -1.5, -0.5, 0.5, 1.5, 2.5]
				tz0 = 2.5
			else:
				coords = [-1.5, -0.5, 0.5, 1.5]
				tz0 = 1.5
			for tx0 in coords:
			    for ty0 in coords:
				(tx, ty, tz) = J(tx0, ty0, tz0)
				if ty0<0:
					tx1 = -tx0
					ty1 = -ty0
					tz1 = tz0
					conj = 1
				else:					
					tx1 = tx0
					ty1 = ty0
					tz1 = tz0
					conj = 0
				for (sx0, sy0, sz0) in [(-1,-1,-1), (1,-1,-1), (-1,1,-1), (1,1,-1)]:
					(sx, sy, sz) = J(sx0, sy0, sz0)
					if not (sx, sy, sz, tx, ty, tz) in alreadyHandled:
						#alreadyHandled.add((sx, sy, sz, tx, ty, tz))
						alreadyHandled[(sx, sy, sz, tx, ty, tz)] = True
						source = InvDir[(sx0, sy0, sz0)][:2]
						print "\t{i%s_%s_%s_%s, d%s_%s_%s_%s, %i, X%s_%s}," % (suf, N[tx], N[ty], N[tz], suf, N[tx1], N[ty1], N[tz1], conj, dir, source)
						n +=1
			print "\t{-1, -1, -1, -1} /* # = %i */" % n
			print "};"
			print 
	
	print '/*'
	print "# = %i" % len(alreadyHandled)
	print '*/'

import sys

print "/* This file is automatically generated by gen_M2L_h_2d.py */"
print "/* DO NOT EDIT! */"
print

gen_M2L_enums()
gen_M2L_tables0()
alreadyHandled = gen_M2L_tables1()
gen_M2L_tables2(alreadyHandled)
alreadyHandled = gen_M2L_tables1(ws=2)
gen_M2L_tables2(alreadyHandled, ws=2)
gen_M2L_tables_reduced()
		
