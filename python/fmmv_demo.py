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

extra_icc_lib = False
dim = 3

from Tkinter import * 
import tkMessageBox

if dim==2:
        import fmmv2d
        import fmmv2df
else:
        import fmmv3d
        import fmmv3df

if extra_icc_lib:
        if dim==2:
	        import fmmv_icc 
	        import fmmvf_icc 
        else:
	        import fmmv_icc 
	        import fmmvf_icc 
	gcc_version ="GNU gcc 4.0.3"
	icc_version ="Intel icc 9.0"


from math import *
from random import randrange
#import webbrowser
import time

from numpy import *
from numpy.random import rand, seed
#from numarray import *
#from numarray.random_array import random, seed

try:
	import Pmw
except:
        tkMessageBox.showerror(
            "Error",
            "FMMV needs the Pmw package.\nSee: pmw.sourceforge.net\n"
        )
	raise ImportError, "could not find Pmw package"



help = {
"random":
"""Select different seeds for the random number generator to generate
different problem data.
Choosing the same seed allows to compare several runs of FMMV on the 
same problem data (while varying parameters of the FMMV algorithms).""",
"sources":
"""Sources Distribution
...""",
"charges":
"""FMMV supports monopole and/or dipole sources.""",

"targets": 
"""FMMV supports independent target locations.  If no targets are specified, 
the sources are taken as target locations.  In this case, the additional 
symmetry is fully exploited.""",

"boundaryConditions":
"""FMMV supports periodic boundary conditions. 
They are implemented as described in
  M.Challlacombe, C.White, M. Head-Gordon: Periodic boundary conditions 
  and the fast multipole method, J.Chem.Phys. 107(23), 10131.
To get results equivalent to those obtained by Ewald summation, 
"extrinsic correction" has to be activated (cf. this ref.).  """,

"gradient":
"""Additionally to the Coulomb (Yukawa) potential, FMMV optionally computes the 
corresponding gradients, i.e., the Coulomb (Yukawa) forces between the particles.""",
"yukawa":
"""FMMV implements a FMM for Coulomb interactions ( 1/r, lambda=0 ) 
and a FMM for Yukawa interactions ( exp(-lambda*r)/r, lambda!=0 ).""",
"WS":
"""FMMV is the first and only implementation of the FMM based on exponential
expansions which supports the well separateness parameter ws=2. 
For ws=2 FMMV optionally supports a reduced scheme where exponential expansions 
are transformed to local expansions already at the parent level whenever possible, 
such that the number of X2X translations is significantly reduced.""",
"adaptive":
"""FMMV implements a non-adaptive FMM (for uniformly distributed particles)
and also an adaptive FMM (for highly non-uniform particle distributions).
Currently, the adaptive FMM is only implemented for free space boundary
conditions and well separateness parameter ws=1.""",
"sorting":
"""While building the oct-tree data structure, all particles are sorted either 
according to Molton order or along a Hilbert space filling curve, such that 
all particles within a certain box are contiguously stored in memory.
The more complex Hilbert order may lead to better performance due to 
improved cache utilization.""",
"scale":
"""FMMV supports the specification of an additional scaling factor which allows
fine control of the average number of particles per leaf box, see
  C.White, M.Head-Gordon: Fractional tiers in fast multipole calculations,
  Chem.Phys.Letters 257 (1996) 647-650.""",
"threads":
"""In FMMV the evaluations of the nearfield and the farfield can be
performed in parallel. This only makes sense if more than one
processor (or a dual core processor) is available.""",
"precision":
"""FMMV comes in a single as well as in a double precision version. On
many architectures (e.g., Intel/SSE, PowerPC/AltiVec) FMMV is fully
vectorized  using  4-way (single precision) or 2-way (double precision) 
SIMD vectorization techniques.""",
"iter":
"""FMMV allows a  specification (typically dependent on the desired accuracy) 
of the number of Newton iterations used for the evaluation of 1/sqrt(x) during 
direct evaluation of the nearfield. The accuracy is respecively 3-4, 7-8, 
or 15-16 digits for 0, 1, or 2 iterations, such that full accuracy means one 
iteration for single precision and two iterations for double precision.""",
"p":
"""Generally, the lengths of the multipole expansions will be the same as those 
of the local expansions. But in some cases  it may be advantageous to choose 
different expansion lengths (e.g., p_local > p_multipole can compensate a 
possible accuracy reduction of the reduced scheme for ws=2).""",
"s":
"""FMMV implements an FMM with fast translation operators based on 
exponential (plane wave) expansions as described in
  H. Cheng, L.Greengard, V.Rokhlin: A fast adaptive multipole algorithm in 
  three dimensions, J.Comp.Phys. 155 (1999), 468-498.
Here a certain integral is approximated by a special generalized Gaussian 
quadrature formula.  Altough in this and other references nodes and weights 
for this quadrature are published only for a few values of the number s of
nodes (and only for well separateness parameter ws=1), FMMV supports all 
values of s within a certain range (and also for ws=2).""",
"aurora":
"""This project was supported by the 
Special Research Program SFB F011 'AURORA' 
(Project 5: Advanced Scientific Computing)
of the Austrian Science Fund FWF.""",
"libVersion":
"""libVersion...""",
}

if dim==2:
        p_max_single = 30
        p_max_double = 40
        s_max_single = 30
        s_max_double = 40
else:
        p_max_single = 20
        p_max_double = 30
        s_max_single = 20
        s_max_double = 30

periodic_dx = 0.78901234567
periodic_dy = 0.45678901234
periodic_dz = 0.12345678901

def get_dist_square(n): 
	return rand(n, 2)
		
def get_dist_circle(n):
	phi = 2.0*pi*rand(n,1)
	return concatenate((
	  0.50 + 0.49999*cos(phi),
	  0.50 + 0.499999*sin(phi)
	  ), 1)
	  # Note: factor .5 instead of .49999 could cause that
	  # some coordinate is exactly equal 1, whereas FMMV
	  # requires all coordinates to be >= 0 and strictly <1.


def get_dist_cube(n): 
	return rand(n, 3)
		
def get_dist_sphere(n):
	phi = 2.0*pi*rand(n,1)
	theta = pi*rand(n,1)
	return concatenate((
	  0.50 + 0.49999*sin(theta)*cos(phi),
	  0.50 + 0.499999*sin(theta)*sin(phi),
	  0.50 + 0.499999*cos(theta)
	  ), 1)
	  # Note: factor .5 instead of .49999 could cause that
	  # some coordinate is exactly equal 1, whereas FMMV
	  # requires all coordinates to be >= 0 and strictly <1.
	
def get_dist_cylinder(n):
	phi = 2.0*pi*rand(n,1)
	return concatenate((
	  rand(n,1),
	  0.5 + 0.05*cos(phi),
	  0.5 + 0.05*sin(phi),
	  ), 1)


if dim==2:
        particle_distributions = [
	        ("uniformly on square" , get_dist_square),
	        ("circle", get_dist_circle)
	        ] 
else:
        particle_distributions = [
	        ("uniformly on cube" , get_dist_cube),
	        ("sphere", get_dist_sphere),
	        ("cylinder h=1, r=0.05", get_dist_cylinder)
	        ] 


class ResultWindow(Frame):
	def __init__(self, master, app, 
		seed0,
		sourcesDist, numberOfParticles, chargesMonopoles, chargesKind, chargesDipoles,
		targetsDist, numberOfTargets,
		boundaryConditions, extrinsicCorrection,
		gradient,
		WS, reducedScheme,
		adaptive, levels, splitThreshold, splitTargetThreshold,
		sorting, scale, beta, threads,
		precision, iter, fullAccuracy, 
		pM0, pL0, s0, libVersion):
		
		self.app = app
		self.seed = seed0
		self.sourcesDist = sourcesDist
		self.numberOfParticles = numberOfParticles
		self.chargesMonopoles = chargesMonopoles
		self.chargesKind = chargesKind
		self.chargesDipoles = chargesDipoles
		self.targetsDist = targetsDist
		self.numberOfTargets = numberOfTargets
		self.boundaryConditions  = boundaryConditions 
		self.extrinsicCorrection = extrinsicCorrection
		self.gradient = gradient
		self.WS = WS
		self.reducedScheme = reducedScheme
		self.adaptive = adaptive
		self.levels = levels
		self.splitThreshold = splitThreshold
		self.splitTargetThreshold = splitTargetThreshold
		self.sorting = sorting
		self.scale = scale
		self.beta = beta
		self.threads = threads
		self.precision = precision
		self.iter = iter
		self.fullAccuracy = fullAccuracy
		self.pM0 = pM0
		self.pL0 = pL0
		self.s0 = s0
		self.libVersion = libVersion
			
		if adaptive==0:
                        splitThreshold = 0
                        splitTargetThreshold = 0
			levels = levels
		else:	
                        splitThreshold = splitThreshold
                        splitTargetThreshold = splitTargetThreshold
			levels = -1
		ws = WS+1	
		pM = int(pM0)
		if pL0=="p_multipole":
			pL = pM 
		else:	
			pL = int(pL0)
		s = int(s0)
		
		n_particles = numberOfParticles
		if targetsDist!="targets=sources":
			n_targets = numberOfTargets
		else:	
			n_targets = n_particles

		call_txt=""
		if  gradient:
			call_txt =  "(pot, grad) = " 
		else:	
			call_txt  = "pot = "

		if libVersion==0:	
		        if precision==0:
			        call_txt += "fmmv%idf.fmmv%idf(\n\tparticles=particles,\n" % (dim, dim)
		        else:	
			        call_txt += "fmmv%id.fmmv%idf(\n\tparticles=particles,\n" % (dim, dim)
		else:	
		        if precision==0:
			        call_txt += "fmmv%idf_icc.fmmv%idf(\n\tparticles=particles,\n" % (dim, dim)
		        else:
			        call_txt += "fmmv%id_icc.fmmv%idf(\n\tparticles=particles,\n" % (dim, dim)
		if chargesMonopoles==1:	
                        call_txt += "\tcharges=charges,\n"
		#else:	
                #        call_txt += "\tcharges=None,\n"
		if chargesDipoles==1:	
                        call_txt += "\tdipoleMoments=dipoleMoments,\n"
		#else:	
                #        print "\tcharges=None,\n"
		if targetsDist!="targets=sources":
			call_txt += "\ttargets=targets,\n"
		#else:	
                #        print "\ttargets=None,\n"
		if  gradient:
                        call_txt += "\tcomputeGradient=True,\n"
                call_txt += "\tsplitThreshold=%i,\n" % splitThreshold
		if splitThreshold==0:
                        call_txt += "\tlevels=%i,\n" % levels
		if targetsDist!="targets=sources" and splitThreshold!=0:
                	call_txt += "\tsplitTargetThreshold=%i,\n" % splitTargetThreshold
		periodicBoundaryConditions = boundaryConditions==1
		if periodicBoundaryConditions:	
                        call_txt += "\tperiodicBoundaryConditions=True,\n"
			if extrinsicCorrection:	
        	                call_txt += "\textrinsicCorrection=True,\n"
		useHilbertOrder = (sorting==1)
		if useHilbertOrder:		
			call_txt += "\tuseHilbertOrder=True,\n"
	        useFarfieldNearfieldThreads = (threads==1)
	        if useFarfieldNearfieldThreads:
			call_txt += "\tuseFarfieldNearfieldThreads=True,\n"
		if scale!=1.0:	
			call_txt += "\tscale=%.2f,\n" % scale
		if scale!=0.0:	
			call_txt += "\tbeta=%.2f,\n" % beta
		if fullAccuracy==1:	
			if precision==0:
				directEvalAccuracy=1
			else:	
				directEvalAccuracy=2
		else:
			directEvalAccuracy = iter
			call_txt += "\tdirectEvalAccuracy=%i,\n" % directEvalAccuracy
		if ws==2:
			call_txt += "\tws=%i,\n" % ws
			if reducedScheme:
				call_txt += "\treducedScheme=True,\n"
		call_txt += "\tpM=%i, pL=%i, s=%i)\n" % ( pM, pL, s)

		Frame.__init__(self, master) 
		self.grid() 
		self.master.title("FMMV Results") 

		row = 0
		text1 = Pmw.ScrolledText(self,
		#text1 = Text(self,
			borderframe=5, # a bit space around the text...
			vscrollmode='dynamic', 
			hscrollmode='dynamic',
			#labelpos='n', label_text='LABEL',
			text_width=40, text_height=30,
			text_wrap='none', # do not break too long lines
			)
		text1.grid(row=row, column=0) 
		text1.insert('end', call_txt)
		text1.configure(text_state=DISABLED)
		row +=1
		self.buttonsFrame = Frame(self)
		self.restoreButton = Button (self.buttonsFrame, text="Restore Parameters", command=self.restoreParameters) 
		self.restoreButton.pack(side=LEFT) 
		self.closeButton = Button (self.buttonsFrame, text="Close", command=self.close)
		self.closeButton.pack(side=LEFT) 
		self.buttonsFrame.grid(row=row, column=0)
		self.master.protocol('WM_DELETE_WINDOW', self.close)
		
		
		ss=int(app.seedVar.get())
		seed([int(ss % 2**31), int(ss / 2**31)])
		
		dist_dict = dict(particle_distributions)

		particles = dist_dict[sourcesDist](n_particles) 
                if precision==0:
                        particles = array(particles, dtype=float32)

		if chargesMonopoles==1:	
			if chargesKind=="general":
				charges = rand(n_particles)
				for i in xrange(0,n_particles/2):
				        charges[2*i] = -charges[2*i+1]
			elif chargesKind=="positive":
				charges = rand(n_particles)
			elif chargesKind=="all +1":
				charges = ones(n_particles)
			elif chargesKind=="+/-1":
				charges = ones(n_particles)
				for i in xrange(0,n_particles/2):
				        charges[2*i] = -charges[2*i+1]
                        if precision==0:
                                charges = array(charges, dtype=float32)
                        
		else:
			charges = None
		
		if chargesDipoles==1:	
			dipoleMoments = rand(n_particles,dim)-0.5
                        if precision==0:
                                dipoleMoments = array(dipoleMoments, dtype=float32)
		else:
			dipoleMoments = None
		
		if targetsDist!="targets=sources":
			targets = dist_dict[targetsDist](n_targets) 
                        if precision==0:
                                targets = array(targets, dtype=float32)
		else:	
			targets = None

		n_direct = min(100, n_targets)
                if dim==2:
		    if libVersion==0:
		        if precision==0:
			        coulomb=fmmv2df.fmmv2df
			        coulomb_direct=fmmv2df.fmmv2df_direct
		        else:	
			        coulomb=fmmv2d.fmmv2d
			        coulomb_direct=fmmv2d.fmmv2d_direct
		    else:		
		        if precision==0:
			        coulomb=fmmv2df_icc.fmmv2df
			        coulomb_direct=fmmv2df_icc.fmmv2df_direct
		        else:	
			        coulomb=fmmv2d_icc.coulomb
			        coulomb_direct=fmmv2d_icc.fmmv2d_direct
                else:
		    if libVersion==0:
		        if precision==0:
			        coulomb=fmmv3df.fmmv3df
			        coulomb_direct=fmmv3df.fmmv3df_direct
		        else:	
			        coulomb=fmmv3d.fmmv3d
			        coulomb_direct=fmmv3d.fmmv3d_direct
		    else:		
		        if precision==0:
			        coulomb=fmmv3df_icc.fmmv3df
			        coulomb_direct=fmmv3df_icc.fmmv3df_direct
		        else:	
			        coulomb=fmmv3d_icc.coulomb
			        coulomb_direct=fmmv3d_icc.fmmv3d_direct

		if periodicBoundaryConditions:
		    if extrinsicCorrection:
                        if dim==2:
			    particles2 = particles + array([[periodic_dx, periodic_dy]])
                        else:
			    particles2 = particles + array([[periodic_dx, periodic_dy, periodic_dz]])
			particles2 = particles2-(particles2>=1)
                        if precision==0:
                                particles2 = array(particles2, dtype=float32)
			if targets==None:
				#targets2 = particles2[:n_direct,:]
				targets2 = None 
			else:
                                if dim==2:
				    targets2 = targets[:n_direct,:] + array([[periodic_dx, periodic_dy]])
                                else:
				    targets2 = targets[:n_direct,:] + array([[periodic_dx, periodic_dy, periodic_dz]])
				targets2 = targets2-(targets2>=1)
                                if precision==0:
                                        targets2 = array(targets2, dtype=float32)

			if gradient:
			   t0 = time.clock(); te0 = time.time()
		    	   (pot_ex, grad_ex) = coulomb(particles=particles2,
                        	charges=charges,
                        	dipoleMoments=dipoleMoments,
				targets=targets2,
                        	splitThreshold=splitThreshold,
                        	splitTargetThreshold=splitTargetThreshold,
                        	levels=levels,
                        	getStatistics=False,
                        	periodicBoundaryConditions=periodicBoundaryConditions,
                        	extrinsicCorrection=extrinsicCorrection,
				useHilbertOrder=useHilbertOrder,
				useFarfieldNearfieldThreads=useFarfieldNearfieldThreads,
				scale=scale,
				beta=beta,
				directEvalAccuracy=directEvalAccuracy,
                        	computeGradient=gradient,
                        	ws=ws,
				reducedScheme=reducedScheme,
                        	pM=pM, pL=pL, s=s )
			   t1 = time.clock(); te1 = time.time()
			   grad_ex = grad_ex[:n_direct,:]
			else:	
			   t0 = time.clock(); te0 = time.time()
		    	   pot_ex = coulomb(particles=particles2,
                        	charges=charges,
                        	dipoleMoments=dipoleMoments,
				targets=targets2,
                        	splitThreshold=splitThreshold,
                        	splitTargetThreshold=splitTargetThreshold,
                        	levels=levels,
                        	getStatistics=False,
                        	periodicBoundaryConditions=periodicBoundaryConditions,
                        	extrinsicCorrection=extrinsicCorrection,
				useHilbertOrder=useHilbertOrder,
				useFarfieldNearfieldThreads=useFarfieldNearfieldThreads,
				scale=scale,
				beta=beta,
				directEvalAccuracy=directEvalAccuracy,
                        	computeGradient=gradient,
                        	ws=ws,
				reducedScheme=reducedScheme,
                        	pM=pM, pL=pL, s=s)
			   t1 = time.clock(); te1 = time.time()
			pot_ex = pot_ex[:n_direct]

				
		else:
			if gradient:
			    td0 = time.clock(); tde0 = time.time()
		    	    (pot_direct, grad_direct) = coulomb_direct(particles=particles,
                        	charges=charges,
                        	dipoleMoments=dipoleMoments,
				targets=targets,
				NTargets=n_direct, 
				beta=beta,
                        	computeGradient=gradient,
				extendedPrecision=False)
			    td1 = time.clock(); tde1 = time.time()
		    	    (pot_ex, grad_ex) = coulomb_direct(particles=particles,
                        	charges=charges,
                        	dipoleMoments=dipoleMoments,
				targets=targets,
				NTargets=n_direct, 
				beta=beta,
                        	computeGradient=gradient,
				extendedPrecision=True)
			else:
			    td0 = time.clock(); tde0 = time.time()
		    	    pot_direct = coulomb_direct(particles=particles,
                        	charges=charges,
                        	dipoleMoments=dipoleMoments,
				targets=targets,
				NTargets=n_direct,
				beta=beta,
                        	computeGradient=gradient,
				extendedPrecision=False)
			    td1 = time.clock(); tde1 = time.time()
		    	    pot_ex = coulomb_direct(particles=particles,
                        	charges=charges,
                        	dipoleMoments=dipoleMoments,
				targets=targets,
				NTargets=n_direct,
				beta=beta,
                        	computeGradient=gradient,
				extendedPrecision=True)

		if gradient:
		    t0 = time.clock(); te0 = time.time()
		    (pot, grad, stat) = coulomb(particles=particles,
                        charges=charges,
                        dipoleMoments=dipoleMoments,
			targets=targets,
                        splitThreshold=splitThreshold,
                        splitTargetThreshold=splitTargetThreshold,
                        levels=levels,
                        getStatistics=True,
                        periodicBoundaryConditions=periodicBoundaryConditions,
                        extrinsicCorrection=extrinsicCorrection,
			useHilbertOrder=useHilbertOrder,
			useFarfieldNearfieldThreads=useFarfieldNearfieldThreads,
			scale=scale,
			beta=beta,
			directEvalAccuracy=directEvalAccuracy,
                        computeGradient=gradient,
                        ws=ws,
			reducedScheme=reducedScheme,
                        pM=pM, pL=pL, s=s )
		    t1 = time.clock(); te1 = time.time()
		else:	
		    t0 = time.clock(); te0 = time.time()
		    (pot, stat) = coulomb(particles=particles,
                        charges=charges,
                        dipoleMoments=dipoleMoments,
			targets=targets,
                        splitThreshold=splitThreshold,
                        splitTargetThreshold=splitTargetThreshold,
                        levels=levels,
                        getStatistics=True,
                        periodicBoundaryConditions=periodicBoundaryConditions,
                        extrinsicCorrection=extrinsicCorrection,
			useHilbertOrder=useHilbertOrder,
			useFarfieldNearfieldThreads=useFarfieldNearfieldThreads,
			scale=scale,
			beta=beta,
			directEvalAccuracy=directEvalAccuracy,
                        computeGradient=gradient,
                        ws=ws,
			reducedScheme=reducedScheme,
                        pM=pM, pL=pL, s=s)
		    t1 = time.clock(); te1 = time.time()

		err_txt=""
                if periodicBoundaryConditions and extrinsicCorrection:
			err_txt += "Test for periodicity:\n"
                        #print "pot-pot_ex", pot-pot_ex
		if not periodicBoundaryConditions or extrinsicCorrection:	
			err_abs = abs(pot[:n_direct]-pot_ex)	
			err_rel = abs((pot[:n_direct]-pot_ex)/pot_ex)	
			err_rel_L2 = sqrt(sum(err_abs**2)/sum(pot_ex**2))
			if not periodicBoundaryConditions:
				err_abs_direct = abs(pot_direct[:n_direct]-pot_ex)	
				err_rel_direct = abs((pot_direct[:n_direct]-pot_ex)/pot_ex)	
				err_rel_L2_direct = sqrt(sum(err_abs_direct**2)/sum(pot_ex**2))
			#for i in range(min(n_direct, n_targets)):
			#	err_txt += "%i  %+10.5e  %+10.5e  %.3e  %.3e\n" \
			#		% (i, pot[i], pot_ex[i], err_abs[i], err_rel[i])
			#err_txt += "\n"
			err_txt += "err_L2(pot) = %.3e (FMM)\n" % err_rel_L2
			if not periodicBoundaryConditions:
				err_txt += "err_L2(pot) = %.3e (direct method)\n" % err_rel_L2_direct
			if gradient:
				err_grad_rel_L2 = sqrt(sum(sum((grad[:n_direct,:]-grad_ex)**2))/sum(sum(grad_ex**2)))
				err_txt += "err_L2(grad) = %.3e (FMM)\n" % err_grad_rel_L2
				if not periodicBoundaryConditions:
					err_grad_rel_L2_direct = sqrt(sum(sum((grad_direct[:n_direct,:]-grad_ex)**2))/sum(sum(grad_ex**2)))
					err_txt += "err_L2(grad) = %.3e (direct method)\n" % err_grad_rel_L2_direct
	
		stat_txt = ""
		for (index, name) in [
			("noOfParticles", "noOfParticles"),
			("noOfSourceLeafBoxes", "noOfSourceLeafBoxes"),
			("averageNoOfParticlesPerLeafBox", "averageNoOfParticlesPerLeafBox"),
			("noOfTargets", "noOfTargets"),
			("noOfTargetLeafBoxes", "noOfTargetLeafBoxes"),
			("averageNoOfTargetsPerLeafBox", "averageNoOfTargetsPerLeafBox"),
			("maxAllocatedMemory", "maxAllocatedMemory"),
			("noOfDirectInteractions", " noOfDirectInteractions")]:
			if index in stat and stat[index]>=0:
				stat_txt += "%s = %i\n" % (name, stat[index])
		stat_txt +="\n"
		stat_txt +="_____Timing:_____user_____real__\n"
			
		
		for (index, name) in [
			("buildTree", "buildTree"),
			("genM", "genM"),
			("M2M", "M2M"),
			("M2L", "M2L"),
			("L2L", "L2L"),
			("evalL", "evalL"),
			("list1", "evalDirect"),
			#("initialize", "initialize"),
			#("evaluate", "evaluate"),
			("nearfield", "nearfield"),
			("farfield", "farfield"),
			#("total", "total"),  # measure total time from Python, see below 
			]:
			if index in stat["time"]:
				stat_txt +="%11s: %7.2fs %7.2fs\n" % (name, stat["time"][index], stat["etime"][index])
		stat_txt += "%11s: %7.2fs %7.2fs\n" % ('total', t1-t0, te1-te0)	
		if not periodicBoundaryConditions:
			td_est = (td1-td0)*(float(n_targets)/float(n_direct))
			tde_est = (tde1-tde0)*(float(n_targets)/float(n_direct))
        		stat_txt += "%11s: %7.2fs %7.2fs (est.)\n" % ('direct meth', td_est, tde_est)	

		#print call_txt
		#print stat_txt
		#print err_txt
		#print stat
		tot_energy_txt = "total energy = %.4e" % (0.5*(sum(pot)))
		
		#####################################################

		text1.configure(text_state=NORMAL)
		text1.insert('end', "\n" + tot_energy_txt + "\n"+stat_txt + "\n"+err_txt)
		text1.configure(text_state=DISABLED)

	def restoreParameters(self):	
		self.app.seedVar.set(self.seed)
		self.app.sourcesDistVar.set(self.sourcesDist)
		self.app.numberOfParticlesVar.set(self.numberOfParticles)
		self.app.chargesMonopolesVar.set( self.chargesMonopoles )
		self.app.chargesKindVar.set(self.chargesKind)
		self.app.chargesDipolesVar.set(self.chargesDipoles)
		self.app.targetsDistVar.set(self.targetsDist)
		self.app.numberOfTargetsVar.set(self.numberOfTargets)
		self.app.boundaryConditionsVar.set(self.boundaryConditions)
		self.app.extrinsicCorrectionVar.set(self.extrinsicCorrection)
		self.app.gradientVar.set(self.gradient)
		self.app.WSVar.set(self.WS)
		self.app.reducedSchemeVar.set(self.reducedScheme)
		self.app.adaptiveVar.set(self.adaptive)
		self.app.levelsVar.set(self.levels)
		self.app.splitThresholdVar.set(self.splitThreshold)
		self.app.splitTargetThresholdVar.set(self.splitTargetThreshold)
		self.app.sortingVar.set(self.sorting)
		self.app.scaleVar.set(self.scale)
		self.app.yukawaVar.set(self.beta)
		self.app.threadsVar.set(self.threads)
		self.app.precisionVar.set(self.precision)
		self.app.iterVar.set(self.iter)
		self.app.fullAccuracyVar.set(self.fullAccuracy)
		self.app.pMVar.set(self.pM0 )
		self.app.pLVar.set(self.pL0)
		self.app.sVar.set(self.s0)
		self.app.libVersionVar.set(self.libVersion)

	def close(self):
		i = 0
		for win in self.app.resultWindows:
			if win==self:
				del self.app.resultWindows[i]
				break
			i += 1
		self.master.destroy()
		



class Application(Frame): 
	def __init__(self, master=None, ):
		Frame.__init__(self, master) 
		self.grid() 
		self.master.title("FMMV Demo Application") 
		self.createWidgets()
		self.resetWidgets()
		self.resultWindows = []

	def showAboutInfo(self):
		tkMessageBox.showinfo(
			"About FMMV",
#A Fully Vektorized Implementation
#of the Fast Multipole Method
#a.k.a.
"""FMMV
The Fastest Multipole Method 
of Vienna

written by
Harald Hofstaetter
http://harald-hofstaetter.at

FMMV is 
(c) 2006-2015 Harald Hofstaetter,
released under the GNU 
General Public License (GPL)
"""
		)	


	def run_fmm(self):
		if not self.seedEntry.valid():
	        	tkMessageBox.showerror(
		            "Error",
        		    "Invalid seed\n"
	        	)
			return

		if not self.numberOfParticlesEntry.valid():
	        	tkMessageBox.showerror(
		            "Error",
        		    "Invalid number of sources\n"
	        	)
			return
		
		if self.targetsDistVar.get()!="targets=sources" and not self.numberOfTargetsEntry.valid():
	        	tkMessageBox.showerror(
		            "Error",
        		    "Invalid number of targets\n"
	        	)
			return

		if self.adaptiveVar.get()==1 and not self.splitThresholdEntry.valid():
	        	tkMessageBox.showerror(
		            "Error",
        		    "Invalid split threshold\n"
	        	)
			return

		if self.adaptiveVar.get()==1 and self.targetsDistVar.get()!="targets=sources" and not self.splitTargetThresholdEntry.valid():
	        	tkMessageBox.showerror(
		            "Error",
        		    "Invalid split threshold (targets)\n"
	        	)
			return

		self.resultWindows.append(ResultWindow(Tk(), self,
		seed0 = self.seedVar.get(),
		sourcesDist = self.sourcesDistVar.get(),
		numberOfParticles = self.numberOfParticlesVar.get(),
		chargesMonopoles = self.chargesMonopolesVar.get(),
		chargesKind = self.chargesKindVar.get(),
		chargesDipoles = self.chargesDipolesVar.get(),
		beta = self.yukawaVar.get(),
		
		targetsDist = self.targetsDistVar.get(),
		numberOfTargets = self.numberOfTargetsVar.get(),

		boundaryConditions = self.boundaryConditionsVar.get(),
		extrinsicCorrection = self.extrinsicCorrectionVar.get(),
		gradient = self.gradientVar.get(),

		WS = self.WSVar.get(),
		reducedScheme = self.reducedSchemeVar.get(),

		adaptive = self.adaptiveVar.get(),
		levels = self.levelsVar.get(),
                splitThreshold = self.splitThresholdVar.get(),
                splitTargetThreshold = self.splitTargetThresholdVar.get(),

		sorting = self.sortingVar.get(),
		scale = self.scaleVar.get(),
		threads = self.threadsVar.get(),
		
		precision = self.precisionVar.get(),
		iter = self.iterVar.get(),
		fullAccuracy = self.fullAccuracyVar.get(),

		pM0 = self.pMVar.get(),
		pL0 = self.pLVar.get(),
		s0 = self.sVar.get(),
		libVersion = self.libVersionVar.get()
		))

	def closeResultWindows(self):
		for win in self.resultWindows:
			win.master.destroy()
		self.resultWindows = []	


	def resetWidgets(self):
		self.seedVar.set(12345678901234L)
		self.seedEntry.configure(validate={"validator": "integer", "min": 1, "max": 99999999999999})
		self.sourcesDistVar.set(particle_distributions[0][0])
		self.numberOfParticlesVar.set(10000)
		self.numberOfParticlesEntry.configure(entryfield_validate={"validator": "integer", "min": 1, "max": 10000000})
		self.chargesMonopolesVar.set(1)
		self.chargesDipolesVar.set(0)
		self.chargesKindVar.set("general") 
		self.yukawaVar.set(0.0)
                if dim==2: # Yukawa not implemented for dim==2
		    self.yukawaEntry.configure(entry_state=DISABLED)
		self.targetsDistVar.set("targets=sources")
		self.numberOfTargetsVar.set(10000)
		self.numberOfTargetsEntry.configure(entryfield_validate={"validator": "integer", "min": 1, "max": 10000000})
		self.chargesMonopolesVar.set(1)
		self.numberOfTargetsEntry.configure(entry_textvariable=self.numberOfParticlesVar)
		self.numberOfTargetsEntry.configure(label_state=DISABLED)
		self.numberOfTargetsEntry.configure(entry_state=DISABLED)
		self.numberOfTargetsEntry.configure(listbox_state=DISABLED)
		self.gradientVar.set(0)
		self.boundaryConditionsVar.set(0)
		self.extrinsicCorrectionVar.set(1)
		self.extrinsicCorrectionCheckbutton.configure(state=DISABLED)
		self.adaptiveSavedValue=0
		self.adaptiveVar.set(0)
		self.adaptiveRadiobutton.configure(state=NORMAL)
		self.levelsVar.set(3)
		self.levelsOptionMenu.configure(label_state=NORMAL)
		self.levelsOptionMenu.configure(menubutton_state=NORMAL)
		self.splitThresholdVar.set(100)
		self.splitTargetThresholdVar.set(100)
		self.splitThresholdEntry.configure(label_state=DISABLED)
		self.splitThresholdEntry.configure(entry_state=DISABLED)
		self.splitThresholdEntry.configure(listbox_state=DISABLED)
		self.splitTargetThresholdEntry.configure(label_state=DISABLED)
		self.splitTargetThresholdEntry.configure(entry_state=DISABLED)
		self.splitTargetThresholdEntry.configure(listbox_state=DISABLED)
		self.splitThresholdEntry.configure(entryfield_validate={"validator": "integer", "min": 2, "max": 10000})
		self.splitTargetThresholdEntry.configure(entryfield_validate={"validator": "integer", "min": 2, "max": 10000})
		self.sortingVar.set(0)
		self.scaleVar.set(1.0)
		self.scaleScale.configure(state=NORMAL)
		self.threadsVar.set(0)
		self.precisionVar.set(0)
		self.iterVar.set(1)
		self.iter0Radiobutton.configure(state=DISABLED)
		self.iter1Radiobutton.configure(state=DISABLED)
		self.iter2Radiobutton.configure(state=DISABLED)
		self.fullAccuracyVar.set(1)
		self.WSVar.set(0)
                if dim==2: # ws=2 not implemented for 2D
		    self.WS1RB.configure(state=DISABLED)
		    self.WS2RB.configure(state=DISABLED)
		self.reducedSchemeVar.set(0)
		self.reducedSchemeCheckbutton.configure(state=DISABLED)
		self.pMVar.set(6)
		self.pLVar.set("p_multipole")
		self.sVar.set(8)
		self.libVersionVar.set(0)

	
	def boundaryConditionsAction(self):
		if self.boundaryConditionsVar.get() == 1: #periodic
			self.extrinsicCorrectionCheckbutton.configure(state=NORMAL)
			self.scaleSavedValue = self.scaleVar.get()
			self.scaleVar.set(1.0)
			self.scaleScale.configure(state=DISABLED)
			self.adaptiveRadiobutton.configure(state=DISABLED)
			if self.WSVar.get() != 1: #ws!=2
				self.adaptiveSavedValue = self.adaptiveVar.get()
			self.adaptiveVar.set(0)
			self.adaptiveAction()
		else: #free space
			self.extrinsicCorrectionCheckbutton.configure(state=DISABLED)
			self.scaleVar.set(self.scaleSavedValue)
			self.scaleScale.configure(state=NORMAL)
			if self.WSVar.get() != 1: #ws!=2
				self.adaptiveRadiobutton.configure(state=NORMAL)
				self.adaptiveVar.set(self.adaptiveSavedValue)
				self.adaptiveAction()

	def adaptiveAction(self):
		if self.adaptiveVar.get() == 1: #adaptive
			self.levelsOptionMenu.configure(label_state=DISABLED)
			self.levelsOptionMenu.configure(menubutton_state=DISABLED)
			self.splitThresholdEntry.configure(label_state=NORMAL)
			self.splitThresholdEntry.configure(entry_state=NORMAL)
			self.splitThresholdEntry.configure(listbox_state=NORMAL)
			if self.targetsDistVar.get()!="targets=sources":
				self.splitTargetThresholdEntry.configure(label_state=NORMAL)
				self.splitTargetThresholdEntry.configure(entry_state=NORMAL)
				self.splitTargetThresholdEntry.configure(listbox_state=NORMAL)
		else: #non adaptive
			self.levelsOptionMenu.configure(label_state=NORMAL)
			self.levelsOptionMenu.configure(menubutton_state=NORMAL)
			self.splitThresholdEntry.configure(label_state=DISABLED)
			self.splitThresholdEntry.configure(entry_state=DISABLED)
			self.splitThresholdEntry.configure(listbox_state=DISABLED)
			self.splitTargetThresholdEntry.configure(label_state=DISABLED)
			self.splitTargetThresholdEntry.configure(entry_state=DISABLED)
			self.splitTargetThresholdEntry.configure(listbox_state=DISABLED)

	def WSAction(self):
		if self.WSVar.get() == 1: #ws==2
			self.reducedSchemeCheckbutton.configure(state=NORMAL)
			self.adaptiveRadiobutton.configure(state=DISABLED)
			if self.boundaryConditionsVar.get() != 1: #!= periodic
				self.adaptiveSavedValue = self.adaptiveVar.get()
			self.adaptiveVar.set(0)
			self.adaptiveAction()
		else: #ws==1
			self.reducedSchemeCheckbutton.configure(state=DISABLED)
			if self.boundaryConditionsVar.get() != 1: #!= periodic
				self.adaptiveRadiobutton.configure(state=NORMAL)
				self.adaptiveVar.set(self.adaptiveSavedValue)
				self.adaptiveAction()

	def chargesMonopolesAction(self):
		if self.chargesMonopolesVar.get()==1:
			self.chargesKindOptionMenu.configure(label_state=NORMAL)
			self.chargesKindOptionMenu.configure(menubutton_state=NORMAL)
		else:
			self.chargesKindOptionMenu.configure(label_state=DISABLED)
			self.chargesKindOptionMenu.configure(menubutton_state=DISABLED)
			if self.chargesDipolesVar.get()==0:
				self.chargesDipolesVar.set(1)

	def chargesDipolesAction(self):
		if self.chargesDipolesVar.get()==0:
			if self.chargesMonopolesVar.get()==0:
				self.chargesMonopolesVar.set(1)
				self.chargesMonopolesAction()

	def fullAccuracyAction(self):
		if self.fullAccuracyVar.get()==1:
			if self.precisionVar.get()==0: # single
				self.iterVar.set(1)
			else:	
				self.iterVar.set(2)
			self.iter0Radiobutton.configure(state=DISABLED)
			self.iter1Radiobutton.configure(state=DISABLED)
			self.iter2Radiobutton.configure(state=DISABLED)
		else:	
			self.iter0Radiobutton.configure(state=NORMAL)
			self.iter1Radiobutton.configure(state=NORMAL)
			self.iter2Radiobutton.configure(state=NORMAL)
		if self.precisionVar.get()==0: # single
			if int(self.pMVar.get())>p_max_single:
				self.pMVar.set(str(p_max_single)) 
			if self.pLVar.get()!="p_multipole" and int(self.pLVar.get())>p_max_single:
				self.pLVar.set(str(p_max_single)) 
			if int(self.sVar.get())>s_max_single:
				self.sVar.set(str(s_max_single)) 

	def pMAction(self, dummy):
		if self.precisionVar.get()==0 and int(self.pMVar.get())>p_max_single:
			self.pMVar.set(str(p_max_single)) 
	def pLAction(self, dummy):
		if self.precisionVar.get()==0 and self.pLVar.get()!="p_multipole" and int(self.pLVar.get())>p_max_single:
			self.pLVar.set(str(p_max_single)) 
	def sAction(self, dummy):
		if self.precisionVar.get()==0 and int(self.sVar.get())>s_max_single:
				self.sVar.set(str(s_max_single)) 

	def targetsAction(self, dummy):
		if self.targetsDistVar.get()=="targets=sources":
			self.numberOfTargetsEntry.configure(label_state=DISABLED)
			self.numberOfTargetsEntry.configure(entry_state=DISABLED)
			self.numberOfTargetsEntry.configure(listbox_state=DISABLED)
			self.numberOfTargetsEntry.configure(entry_textvariable=self.numberOfParticlesVar)
			self.splitTargetThresholdEntry.configure(label_state=DISABLED)
			self.splitTargetThresholdEntry.configure(entry_state=DISABLED)
			self.splitTargetThresholdEntry.configure(listbox_state=DISABLED)
		else:
			self.numberOfTargetsEntry.configure(label_state=NORMAL)
			self.numberOfTargetsEntry.configure(entry_state=NORMAL)
			self.numberOfTargetsEntry.configure(listbox_state=NORMAL)
			self.numberOfTargetsEntry.configure(entry_textvariable=self.numberOfTargetsVar)
			if self.adaptiveVar.get() == 1: #adaptive
				self.splitTargetThresholdEntry.configure(label_state=NORMAL)
				self.splitTargetThresholdEntry.configure(entry_state=NORMAL)
				self.splitTargetThresholdEntry.configure(listbox_state=NORMAL)

	def newSeedAction(self):
		self.seedVar.set(randrange(27814431486575L))
		self.seedEntry.configure(validate={"validator": "integer", "min": 1, "max": 99999999999999})

	def createWidgets(self):
		self.balloon = Pmw.Balloon(self.master)
		self.titleFrame = Frame(self, borderwidth=2, relief='groove')
		
		#self.titleLabel = Label(self.titleFrame, text="FMMV - A Fully Vectorized Implementation of the Fast Multipole Method", font="times 18 bold")
		self.titleLabel = Label(self.titleFrame, text="FMMV - The Fastest Multipole Method of Vienna", font="times 18 bold")
		self.titleLabel.pack(side='left', ipady=3)

		self.problemFrame = Frame(self, borderwidth=2, relief='groove')
		row = 0
		self.problemFrame.columnconfigure(0, minsize=260)
		self.problemFrame.columnconfigure(1, minsize=370)
		self.problemLabel = Label(self.problemFrame, text="Problem Specification", font="times 18 bold")
		self.problemLabel.grid(row=row, column=0, sticky=W)

		row+=1
		self.seedVar = IntVar(self)
		self.seedVar.set(12345678901234)
		self.randomLabel = Label(self.problemFrame, text="Random Data: ")
		self.randomLabel.grid(row=row, column=0, sticky=E)
		self.balloon.bind(self.randomLabel, help["random"])
		self.randomFrame = Frame(self.problemFrame)
		self.seedLabel = Label(self.randomFrame, text="seed = ")
		self.seedLabel.pack(side=LEFT)
		self.seedEntry = Pmw.EntryField(self.randomFrame, entry_textvariable=self.seedVar, entry_width=14, validate={"validator": "integer", "min": 1, "max": 99999999999999})
		self.seedEntry.pack(side=LEFT)
		self.newSeedButton = Button(self.randomFrame, text="new seed", command=self.newSeedAction)
		self.newSeedButton.pack(side=LEFT)
		self.randomFrame.grid(row=row, column=1, sticky=W)
		
		row+=1
		self.sourcesLabel = Label(self.problemFrame, text="Source Distribution: ")
		self.sourcesLabel.grid(row=row, column=0, sticky=E)
		#self.balloon.bind(self.sourcesLabel, help["sources"])
		self.sourcesFrame = Frame(self.problemFrame)
		self.sourcesDistVar = StringVar(self)
		self.sourcesDistOptionMenu = Pmw.OptionMenu(self.sourcesFrame, menubutton_textvariable=self.sourcesDistVar, items=[a for (a,b) in particle_distributions])
		self.sourcesDistOptionMenu.pack(side=LEFT)
		self.numberOfParticlesVar = IntVar(self)
		self.numberOfParticlesVar.set(10000)
		self.numberOfParticlesEntry = Pmw.ComboBox(self.sourcesFrame, labelpos="w", label_text="number of sources:", entry_textvariable=self.numberOfParticlesVar, entry_width=8, entryfield_validate={"validator": "integer", "min": 1, "max": 10000000}, scrolledlist_items=[10, 100, 1000, 10000, 100000, 1000000], listheight=110, dropdown=True)
		self.numberOfParticlesEntry.pack(side=LEFT)
		self.sourcesFrame.grid(row=row, column=1, sticky=W)

		
		
		row += 1
		self.chargesMonopolesVar = IntVar(self)
		self.chargesDipolesVar = IntVar(self)
		self.chargesKindVar = StringVar(self)
		self.chargesLabel = Label(self.problemFrame, text="Charges: ")
		self.chargesLabel.grid(row=row, column=0, sticky=E)
		self.balloon.bind(self.chargesLabel, help["charges"])
		self.chargesFrame = Frame(self.problemFrame)
		self.chargesMonopolesCheckbutton = Checkbutton(self.chargesFrame, text="monopoles", variable=self.chargesMonopolesVar, command=self.chargesMonopolesAction)
		self.chargesDipolesCheckbutton = Checkbutton(self.chargesFrame, text="dipoles", variable=self.chargesDipolesVar, command=self.chargesDipolesAction)
		self.chargesMonopolesCheckbutton.pack(side=LEFT)
		self.chargesKindOptionMenu = Pmw.OptionMenu(self.chargesFrame, labelpos="w", label_text="/ kind:", menubutton_textvariable=self.chargesKindVar, items=["general", "positive", "all +1", "+/-1"] )
		self.chargesKindOptionMenu.pack(side=LEFT)

		self.chargesDipolesCheckbutton.pack(side=LEFT)
		self.chargesFrame.grid(row=row, column=1, sticky=W)

                row +=1
		self.yukawaVar = DoubleVar(self)
		self.yukawaVar.set(0)
		self.yukawaLabel = Label(self.problemFrame, text="Yukawa Potential Coefficient: ")
		self.yukawaLabel.grid(row=row, column=0, sticky=E)
		self.balloon.bind(self.yukawaLabel, help["yukawa"])
		self.yukawaFrame = Frame(self.problemFrame)
		self.yukawaLabel = Label(self.yukawaFrame, text="lambda = ")
		self.yukawaLabel.pack(side=LEFT)
		self.yukawaEntry = Pmw.EntryField(self.yukawaFrame, entry_textvariable=self.yukawaVar, entry_width=14, validate={"validator": "real", "min": -10, "max": 10})
		self.yukawaEntry.pack(side=LEFT)
		self.yukawaFrame.grid(row=row, column=1, sticky=W)



		row +=1
		self.targetsLabel = Label(self.problemFrame, text="Target Distribution: ")
		self.targetsLabel.grid(row=row, column=0, sticky=E)
		self.balloon.bind(self.targetsLabel, help["targets"])
		self.targetsFrame = Frame(self.problemFrame)
		self.targetsDistVar = StringVar(self)
		self.targetsDistOptionMenu = Pmw.OptionMenu(self.targetsFrame, menubutton_textvariable=self.targetsDistVar, items=["targets=sources"] + [a for (a,b) in particle_distributions], command=self.targetsAction)
		self.targetsDistOptionMenu.pack(side=LEFT)
		self.numberOfTargetsVar = IntVar(self)
		self.numberOfTargetsVar.set(10000)
		self.numberOfTargetsEntry = Pmw.ComboBox(self.targetsFrame, labelpos="w", label_text="number of targets:", entry_textvariable=self.numberOfParticlesVar, entry_width=8, entryfield_validate={"validator": "integer", "min": 1, "max": 10000000}, scrolledlist_items=[10, 100, 1000, 10000, 100000, 1000000], listheight=110, dropdown=True)
		self.numberOfTargetsEntry.pack(side=LEFT)
		self.targetsFrame.grid(row=row, column=1, sticky=W)

		row+=1
		self.boundaryConditionsVar = IntVar()
		self.extrinsicCorrectionVar = IntVar()
		self.boundaryConditionsLabel = Label(self.problemFrame, text="Boundary Conditions: ")
		self.boundaryConditionsLabel.grid(row=row, column=0, sticky=E)
		self.balloon.bind(self.boundaryConditionsLabel, help["boundaryConditions"])
		self.boundaryConditionsFrame = Frame(self.problemFrame)
		Radiobutton(self.boundaryConditionsFrame, text="free space", variable=self.boundaryConditionsVar, value=0, command=self.boundaryConditionsAction).pack(side=LEFT)
		Radiobutton(self.boundaryConditionsFrame, text="periodic", variable=self.boundaryConditionsVar, value=1, command=self.boundaryConditionsAction).pack(side=LEFT)
		self.extrinsicCorrectionCheckbutton = Checkbutton(self.boundaryConditionsFrame, text="extrinsic correction", variable=self.extrinsicCorrectionVar)
		self.extrinsicCorrectionCheckbutton.pack(side=LEFT)
		self.boundaryConditionsFrame.grid(row=row, column=1, sticky=W)
		
		
		row +=1
		self.gradientVar = IntVar(self)
		self.gradientLabel = Label(self.problemFrame, text="Gradients: ")
		self.gradientLabel.grid(row=row, column=0, sticky=E)
		self.balloon.bind(self.gradientLabel, help["gradient"])
		self.gradientFrame = Frame(self.problemFrame)
		self.gradientCheckbutton = Checkbutton(self.gradientFrame, text="compute gradients", variable=self.gradientVar)
		self.gradientCheckbutton.pack(side=LEFT)
		self.gradientFrame.grid(row=row, column=1, sticky=W)

		self.algorithmFrame = Frame(self, borderwidth=2, relief='groove')
		self.algorithmFrame.columnconfigure(0, minsize=260)
		self.algorithmFrame.columnconfigure(1, minsize=370)
		row=0
		self.algorithmLabel = Label(self.algorithmFrame, text="FMM Parameter Specification", font="times 18 bold")
		self.algorithmLabel.grid(row=row, column=0, sticky=W)
		

		row+=1
		self.WSVar = IntVar()
		self.reducedSchemeVar = IntVar()
		self.WSLabel = Label(self.algorithmFrame, text="Well Separateness Parameter: ")
		self.WSLabel.grid(row=row, column=0, sticky=E)
		self.balloon.bind(self.WSLabel, help["WS"])
		self.WSFrame = Frame(self.algorithmFrame)
		self.WS1RB = Radiobutton(self.WSFrame, text="ws=1", variable=self.WSVar, value=0, command=self.WSAction)
                self.WS1RB.pack(side=LEFT)
		self.WS2RB = Radiobutton(self.WSFrame, text="ws=2", variable=self.WSVar, value=1, command=self.WSAction)
                self.WS2RB.pack(side=LEFT)
		self.reducedSchemeCheckbutton = Checkbutton(self.WSFrame, text="reduced scheme", variable=self.reducedSchemeVar)
		self.reducedSchemeCheckbutton.pack(side=LEFT, ipadx=15)
		self.WSFrame.grid(row=row, column=1, sticky=W)

		row+=1
		self.adaptiveVar = IntVar()
		self.levelsVar = IntVar()
		self.splitThresholdVar = IntVar()
		self.splitTargetThresholdVar = IntVar()
		self.adaptiveLabel = Label(self.algorithmFrame, text="Adaptive Scheme: ")
		self.adaptiveLabel.grid(row=row, column=0, sticky=E)
		self.balloon.bind(self.adaptiveLabel, help["adaptive"])
		self.adaptiveFrame = Frame(self.algorithmFrame)
		self.adaptiveRadiobutton = Radiobutton(self.adaptiveFrame, text="adaptive", variable=self.adaptiveVar, value=1, command=self.adaptiveAction)
		self.adaptiveRadiobutton.pack(side=LEFT)
		Radiobutton(self.adaptiveFrame, text="non adaptive", variable=self.adaptiveVar, value=0, command=self.adaptiveAction).pack(side=LEFT)
		self.levelsOptionMenu = Pmw.OptionMenu(self.adaptiveFrame, labelpos="w", label_text="/ levels:", menubutton_textvariable=self.levelsVar, items=range(0,8))
		self.levelsOptionMenu.pack(side=LEFT)
		self.adaptiveFrame.grid(row=row, column=1, sticky=W)

		row+=1
		self.adaptiveFrame1 = Frame(self.algorithmFrame)
		self.splitThresholdEntry = Pmw.ComboBox(self.adaptiveFrame1, labelpos="w", label_text="split threshold :", entry_textvariable=self.splitThresholdVar, entry_width=4, entryfield_validate={"validator": "integer", "min": 2, "max": 10000}, scrolledlist_items=range(10,260,10), listheight=110, dropdown=True)
		self.splitThresholdEntry.pack(side=LEFT)
		self.splitThresholdEntry.configure(entryfield_validate={"validator": "integer", "min": 2, "max": 10000})
		self.adaptiveFrame1.grid(row=row, column=1, sticky=W)
		self.splitTargetThresholdEntry = Pmw.ComboBox(self.adaptiveFrame1, labelpos="w", label_text="split threshold (targets):", entry_textvariable=self.splitTargetThresholdVar, entry_width=4, entryfield_validate={"validator": "integer", "min": 2, "max": 10000}, scrolledlist_items=range(10,260,10), listheight=110, dropdown=True)
		self.splitTargetThresholdEntry.configure(entryfield_validate={"validator": "integer", "min": 2, "max": 10000})
		self.numberOfParticlesEntry.pack(side=LEFT)
		self.splitTargetThresholdEntry.pack(side=LEFT)
		self.adaptiveFrame1.grid(row=row, column=1, sticky=W)

		row+=1
		self.sortingVar = IntVar()
		self.sortingLabel = Label(self.algorithmFrame, text="Particle Sorting: ")
		self.sortingLabel.grid(row=row, column=0, sticky=E)
		self.balloon.bind(self.sortingLabel, help["sorting"])
		self.sortingFrame = Frame(self.algorithmFrame)
		Radiobutton(self.sortingFrame, text="Molton order", variable=self.sortingVar, value=0).pack(side=LEFT)
		Radiobutton(self.sortingFrame, text="Hilbert Order", variable=self.sortingVar, value=1).pack(side=LEFT)
		self.sortingFrame.grid(row=row, column=1, sticky=W)

		row+=1
		self.scaleLabel = Label(self.algorithmFrame, text="Fractional Tiers Scaling: ")
		self.scaleLabel.grid(row=row, column=0, sticky=E)
		self.balloon.bind(self.scaleLabel, help["scale"])
		self.scaleFrame = Frame(self.algorithmFrame)
		self.scaleVar = DoubleVar()
		self.scaleScale = Scale(self.scaleFrame, from_=0.5, to=1, resolution=0.01, length=200, orient=HORIZONTAL, variable=self.scaleVar)
		self.scaleScale.pack(side=LEFT)
		self.scaleFrame.grid(row=row, column=1, sticky=W)

		row+=1
		self.threadsVar = IntVar()
		self.threadsLabel = Label(self.algorithmFrame, text="Concurrent Threads: ")
		self.threadsLabel.grid(row=row, column=0, sticky=E)
		self.balloon.bind(self.threadsLabel, help["threads"])
		self.threadsFrame = Frame(self.algorithmFrame)
		Radiobutton(self.threadsFrame, text="one thread only", variable=self.threadsVar, value=0).pack(side=LEFT)
		Radiobutton(self.threadsFrame, text="farfield/nearfield threads", variable=self.threadsVar, value=1).pack(side=LEFT)
		self.threadsFrame.grid(row=row, column=1, sticky=W)

		row+=1
		self.precisionVar = IntVar()
		self.precisionFrame = Frame(self.algorithmFrame)
		self.precisionLabel = Label(self.algorithmFrame, text="Floating Point Precision: ")
		self.precisionLabel.grid(row=row, column=0, sticky=E)
		self.balloon.bind(self.precisionLabel, help["precision"])
		self.singlePrecisionRadiobutton = Radiobutton(self.precisionFrame, text="single", variable=self.precisionVar, value=0, command=self.fullAccuracyAction)
		self.singlePrecisionRadiobutton.pack(side=LEFT) 
		self.doublePrecisionRadiobutton = Radiobutton(self.precisionFrame, text="double", variable=self.precisionVar, value=1, command=self.fullAccuracyAction)
		self.doublePrecisionRadiobutton.pack(side=LEFT) 
		self.precisionFrame.grid(row=row, column=1, sticky=W)

		row+=1
		self.iterVar = IntVar()
		self.fullAccuracyVar = IntVar()
		self.iterLabel = Label(self.algorithmFrame, text="Accuracy of direct evaluations: ")
		self.iterLabel.grid(row=row, column=0, sticky=E)
		self.balloon.bind(self.iterLabel, help["iter"])
		self.iterFrame = Frame(self.algorithmFrame)
		self.iter0Radiobutton = Radiobutton(self.iterFrame, text="0", variable=self.iterVar, value=0)
		self.iter0Radiobutton.pack(side=LEFT)
		self.iter1Radiobutton = Radiobutton(self.iterFrame, text="1", variable=self.iterVar, value=1)
		self.iter1Radiobutton.pack(side=LEFT)
		self.iter2Radiobutton = Radiobutton(self.iterFrame, text="2", variable=self.iterVar, value=2)
		self.iter2Radiobutton.pack(side=LEFT)
		self.fullAccuracyCheckbutton = Checkbutton(self.iterFrame, text="full accuracy", variable=self.fullAccuracyVar, command=self.fullAccuracyAction)
		self.fullAccuracyCheckbutton.pack(side=LEFT, ipadx=15)
		self.iterFrame.grid(row=row, column=1, sticky=W)

		row+=1
		self.pMVar = IntVar()
		self.pLVar = StringVar()
		self.pFrame = Frame(self.algorithmFrame)
		self.pLabel = Label(self.algorithmFrame, text="Lengths of Multipole/Local Expansions: ")
		self.pLabel.grid(row=row, column=0, sticky=E)
		self.balloon.bind(self.pLabel, help["p"])
		self.pMOptionMenu = Pmw.OptionMenu(self.pFrame, labelpos="w", label_text="p_multipole =", menubutton_textvariable=self.pMVar, items=range(0,p_max_double+1), command=self.pMAction)
		self.pMOptionMenu.pack(side=LEFT)
		self.pLOptionMenu = Pmw.OptionMenu(self.pFrame, labelpos="w", label_text=" p_local =", menubutton_textvariable=self.pLVar, items=["p_multipole"]+range(0,p_max_double+1), command=self.pLAction)
		self.pLOptionMenu.pack(side=LEFT)
		self.pFrame.grid(row=row, column=1, sticky=W)

		row+=1
		self.sVar = IntVar()
		self.sFrame = Frame(self.algorithmFrame)
		self.sLabel = Label(self.algorithmFrame, text="Order of (generalized) Gauss Quadrature: ")
		self.sLabel.grid(row=row, column=0, sticky=E)
		self.balloon.bind(self.sLabel, help["s"])
		self.sOptionMenu = Pmw.OptionMenu(self.sFrame, labelpos="w", label_text="s =", menubutton_textvariable=self.sVar, items=range(2,s_max_double+1), command=self.sAction)
		self.sOptionMenu.pack(side=LEFT)

		self.sFrame.grid(row=row, column=1, sticky=W)

		self.libVersionVar = IntVar()
		if extra_icc_lib:
			self.techFrame = Frame(self, borderwidth=2, relief='groove')
			self.techFrame.columnconfigure(0, minsize=260)
			self.techFrame.columnconfigure(1, minsize=370)

			row = 0
			self.libVersionFrame = Frame(self.techFrame)
			self.libVersionLabel = Label(self.techFrame, text="FMMV library compiled with: ")
			self.libVersionLabel.grid(row=row, column=0, sticky=E)
			self.balloon.bind(self.libVersionLabel, help["libVersion"])
			self.gccRadiobutton = Radiobutton(self.libVersionFrame, text=gcc_version, variable=self.libVersionVar, value=0, command=self.fullAccuracyAction)
			self.gccRadiobutton.pack(side=LEFT) 
			self.iccRadiobutton = Radiobutton(self.libVersionFrame, text=icc_version, variable=self.libVersionVar, value=1, command=self.fullAccuracyAction)
			self.iccRadiobutton.pack(side=LEFT) 
			self.libVersionFrame.grid(row=row, column=1, sticky=W)

		self.buttonsFrame = Frame(self)

		self.runButton = Button ( self.buttonsFrame, text="Run", 
			command=self.run_fmm ) 
		self.runButton.pack(side=LEFT) 
		self.resetButton = Button ( self.buttonsFrame, text="Reset Parameters", 
			command=self.resetWidgets ) 
		self.resetButton.pack(side=LEFT) 
		self.closeButton = Button ( self.buttonsFrame, text="Close Result Windows", 
			command=self.closeResultWindows ) 
		self.closeButton.pack(side=LEFT) 

		#self.helpButton = Button ( self.buttonsFrame, text="Help", 
		#	command=webbrowser.open_new("google.at")) 
		#self.helpButton.pack(side=LEFT) 

		self.aboutButton = Button ( self.buttonsFrame, text="About", 
			command=self.showAboutInfo ) 
		self.aboutButton.pack(side=LEFT) 
		
		self.quitButton = Button ( self.buttonsFrame, text="Quit", 
			command=self.quit ) 
		self.quitButton.pack(side=LEFT)
		try:
			self.aurorapic=PhotoImage(file="aurora.gif")
			self.auroraLabel = Label(self.buttonsFrame, image = self.aurorapic)
			self.auroraLabel.pack(side=LEFT)
			self.balloon.bind(self.auroraLabel, help["aurora"])
		except:
			None
                

		self.titleFrame["relief"] = "groove"
		self.titleFrame["borderwidth"] = 2
		self.problemFrame["relief"] = "groove"
		self.problemFrame["borderwidth"] = 2
		self.algorithmFrame["relief"] = "groove"
		self.algorithmFrame["borderwidth"] = 2

		self.titleFrame.pack(fill='x', pady=1)
		self.problemFrame.pack(fill='x', pady=1)
		self.algorithmFrame.pack(fill='x', pady=1)
		if extra_icc_lib:
			self.techFrame["relief"] = "groove"
			self.techFrame["borderwidth"] = 2
			self.techFrame.pack(fill='x', pady=1)
		self.buttonsFrame.pack()

app = Application() 
app.mainloop()
