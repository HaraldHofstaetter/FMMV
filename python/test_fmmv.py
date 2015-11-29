import sys
import time
from numpy import *
from math import *

def set_from_command_line(varname, var):
    if type(var)==bool:
        if varname in sys.argv:
            print "option '%s'" % varname
            return True
    else:
        lhs = [i.split('=')[0] for i in sys.argv]
        if varname in lhs:
            var_new = sys.argv[lhs.index(varname)].split('=')[1]
            print "changed %s from %s to %s " % (varname, str(var), var_new)
            if type(var)==int:
                return int(var_new)
            elif type(var)==float:
                return float(var_new)
            else:    
                return var_new
    return var


def test_it():
    dim = set_from_command_line("dim", 3)
    use_single_precision = set_from_command_line("single", False)
    with_dipoleMoments =  set_from_command_line("dipoleMoments", False)
    with_gradients = set_from_command_line("gradients", False)
    printResult = set_from_command_line("printResult", False)
    splitCall = set_from_command_line("splitCall", False)
    directEvalExtendedPrecision = set_from_command_line("directEvalExtendedPrecision", False)

    NParticles = set_from_command_line("NParticles", 10000)
    NTargets = set_from_command_line("NTargets", 0)
    NDirect = set_from_command_line("NDirect", 100)

    pM = set_from_command_line("pM", 6)
    pL = set_from_command_line("pL", 6)
    s = set_from_command_line("s", 8)
    
    ws = set_from_command_line("ws", 1)
    reducedScheme = set_from_command_line("reducedScheme", False)
   
    scale = set_from_command_line("scale", 1.0)
    splitThreshold = set_from_command_line("splitThreshold", 190)
    splitTargetThreshold =  set_from_command_line("splitTargetThreshold", 190)
    levels =  set_from_command_line("levels", 51)
    directEvalThreshold = set_from_command_line("directEvalThreshold", -1)
    periodicBoundaryConditions = set_from_command_line("periodicBoundaryConditions", False)
    extrinsicCorrection = set_from_command_line("extrinsicCorrection", False)
    useHilbertOrder = set_from_command_line("useHilbertOrder", False)
    if use_single_precision:
        if dim==2:
            import fmmv2df
        else:
            import fmmv3df
        my_float = float32
        directEvalAccuracy = set_from_command_line("directEvalAccuracy", 1)
        exAcc = 1
    else:    
        if dim==2:
            import fmmv2d
        else:
            import fmmv3d
        my_float = float64
        directEvalAccuracy = set_from_command_line("directEvalAccuracy", 2)
        exAcc = 2
    useFarfieldNearfieldThreads = set_from_command_line("useFarfieldNearfieldThreads", False)

    random.seed([123,1234])
    particles = array(random.rand(NParticles,dim), dtype=my_float)
    charges = ones(NParticles, dtype=my_float)

    with_targets = (NTargets>0)
    if with_targets:
        targets = array(random.rand(NTargets,dim), dtype=my_float)
    else:
        targets = None
        NTargets = NParticles

    if with_dipoleMoments:
        dipoleMoments = array(random.rand(NParticles,dim), dtype=my_float)
    else:
        dipoleMoments = None

    if dim==2:
	        if use_single_precision:
		        coulomb = fmmv2df.fmmv2df
		        coulomb_initialize = fmmv2df.fmmv2df_initialize
		        coulomb_evaluate = fmmv2df.fmmv2df_evaluate
		        coulomb_finalize = fmmv2df.fmmv2df_finalize
		        coulomb_direct = fmmv2df.fmmv2df_direct
	        else:	
		        coulomb = fmmv2d.fmmv2d
		        coulomb_initialize = fmmv2d.fmmv2d_initialize
		        coulomb_evaluate = fmmv2d.fmmv2d_evaluate
		        coulomb_finalize = fmmv2d.fmmv2d_finalize
		        coulomb_direct = fmmv2d.fmmv2d_direct
    else:
	        if use_single_precision:
		        coulomb = fmmv3df.fmmv3df
		        coulomb_initialize = fmmv3df.fmmv3df_initialize
		        coulomb_evaluate = fmmv3df.fmmv3df_evaluate
		        coulomb_finalize = fmmv3df.fmmv3df_finalize
		        coulomb_direct = fmmv3df.fmmv3df_direct
	        else:	
		        coulomb = fmmv3d.fmmv3d
		        coulomb_initialize = fmmv3d.fmmv3d_initialize
		        coulomb_evaluate = fmmv3d.fmmv3d_evaluate
		        coulomb_finalize = fmmv3d.fmmv3d_finalize
		        coulomb_direct = fmmv3d.fmmv3d_direct

        

    if splitCall:
        (handle, stat) = coulomb_initialize(particles=particles,
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
                directEvalAccuracy=directEvalAccuracy,
                computeGradient=with_gradients,
                dipoleSources=with_dipoleMoments,
                ws=ws,
                reducedScheme=reducedScheme,
                pM=pM, pL=pL, s=s )
        if with_gradients:
            (pot, grad, stat) = coulomb_evaluate(handle=handle,
                charges=charges,
                dipoleMoments=dipoleMoments,
                getStatistics=True)
        else:
            (pot, stat) = coulomb_evaluate(handle=handle,
                charges=charges,
                dipoleMoments=dipoleMoments,
                getStatistics=True)
        stat =  coulomb_finalize(handle=handle,
                getStatistics=True)
    
    else:
        if with_gradients:
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
                directEvalAccuracy=directEvalAccuracy,
                computeGradient=True,
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
                directEvalAccuracy=directEvalAccuracy,
                computeGradient=False,
                ws=ws,
                reducedScheme=reducedScheme,
                pM=pM, pL=pL, s=s)
            t1 = time.clock(); te1 = time.time()
    if with_gradients:
            (exPot, exGrad) = coulomb_direct(particles=particles,
                charges=charges,
                dipoleMoments=dipoleMoments,
                targets=targets,
                directEvalAccuracy=exAcc,
                extendedPrecision=directEvalExtendedPrecision,
                computeGradient=True
                )
    else:	
            exPot = coulomb_direct(particles=particles,
                charges=charges,
                dipoleMoments=dipoleMoments,
                targets=targets,
                directEvalAccuracy=exAcc,
                extendedPrecision=directEvalExtendedPrecision,
                computeGradient=False
                )
    
    errL2 = sqrt(sum((pot-exPot)**2)/sum(exPot**2))
    print "err_L2(pot) = %12.4e" % errL2,
    if with_gradients:
        errL2grad = sqrt(sum((grad-exGrad)**2)/sum(exGrad**2))
        print "err_L2(grad) = %12.4e" % errL2grad,
    print    
test_it()    
