from fmmv_setup import 	name,additional_options, additional_statistics, float_types, dim, eval_direct_number_of_accuracies, drivers
from fmmv_setup import charge_kinds as charge_kinds0
from fmmv_setup import pot_kinds as pot_kinds0 
from fmmv_setup import charge_pot_combs as charge_pot_combs0
from fmmv_setup import charge_combs as charge_combs0
from fmmv_setup import pot_combs as pot_combs0


def Aa(name):
	return name[0].upper()+name[1:]

papi = True

def get_subdictionaryies(subset):
        charge_pot_combs =  {x: charge_pot_combs0[x] for x in subset}
        s = [x[0] for x in charge_pot_combs.values()]
        charge_combs = {x: charge_combs0[x] for x in s}
        s = sum([x for x in charge_combs.values()], [])
        charge_kinds = {x: charge_kinds0[x] for x in s}
        s = [x[1] for x in charge_pot_combs.values()]
        pot_combs = {x: pot_combs0[x] for x in s}
        s = sum([x for x in pot_combs.values()], [])
        pot_kinds = {x: pot_kinds0[x] for x in s}
        return (charge_kinds, pot_kinds, charge_combs, pot_combs, charge_pot_combs)

def print_charge_args(charge_kinds, float_type):
	charge_kinds_sorted = range(len(charge_kinds))
	for i in charge_kinds:
		charge_kinds_sorted[charge_kinds[i][2]] = i
	for i in charge_kinds_sorted:
		n = charge_kinds[i][0]
                if float_type==None:
		    print "\t\t%s," % n
                else:
		    d = len(charge_kinds[i][1])
		    if d==1:
			    print "\t%s %s[]," % (float_type, n)
		    else:
			    print "\t%s %s[][%i]," % (float_type, n, d)

def print_pot_args(pot_kinds, float_type):
	pot_kinds_sorted = range(len(pot_kinds)) 
	for i in pot_kinds:
		pot_kinds_sorted[pot_kinds[i][2]] = i
	for i in pot_kinds_sorted:
		n = pot_kinds[i][0]
                if float_type==None:
		    print "\t\t%s," % n
                else:
		    d = len(pot_kinds[i][1])
		    if d==1:
			    print "\t%s %s[]," % (float_type, n)
		    else:
			    print "\t%s %s[][%i]," % (float_type, n, d)

def gen_header():
	print "#ifndef _FMMV_%s_H_" % name["generic"].upper()
	print "#define _FMMV_%s_H_" % name["generic"].upper()
	print 
	print "#define STAT_MAX 20"
	if papi: 
		print "#define MAX_NUM_PAPI_EVENTS 10"
	
	print """
enum StatStep {
        STAT_TOTAL,
        STAT_BUILD_TREE, 
	STAT_GEN_M,
        STAT_M2M, 
	STAT_M2L, 
	STAT_L2L,
        STAT_EVAL_L,
        STAT_LIST1, 
	STAT_LIST3, 
	STAT_LIST4, 
	STAT_LIST34,
        STAT_FARFIELD, 
	STAT_NEARFIELD,
	STAT_INITIALIZE,
	STAT_EVALUATE,
	STAT_FINALIZE,
	_STAT_LAST_ /* mark last entry, do not remove! */	
};"""
	
	print """
struct FmmvStatistics {
	double time[STAT_MAX];
	double etime[STAT_MAX];

	int PAPIeventSet;
	long long int PAPIvalues[STAT_MAX][MAX_NUM_PAPI_EVENTS];

	int pM;
	int pL;
	int s_eps;
	int s_exp;
		
	int noOfParticles;
	int noOfTargets;
	int noOfSourceLevels;
	int noOfTargetLevels;
	
	int maxNoOfStoredXin;
	int maxAllocatedMemory;
	long long int noOfDirectInteractions;

	int noOfSourceBoxes;
	int noOfTargetBoxes;
	int noOfSourceLeafBoxes;
	int noOfTargetLeafBoxes;
	float averageNoOfParticlesPerLeafBox;
	float averageNoOfTargetsPerLeafBox;

	int noOfParticlesInLevel[52];
	int noOfTargetsInLevel[52];
	int noOfSourceBoxesInLevel[52];
	int noOfTargetBoxesInLevel[52];
	int noOfSourceLeafBoxesInLevel[52];
	int noOfTargetLeafBoxesInLevel[52];
	float averageNoOfParticlesPerLeafBoxInLevel[52];
	float averageNoOfTargetsPerLeafBoxInLevel[52];
	"""
	for data in additional_statistics:
		print "\t%s;" % data
	print "};"

	print """
void printFmmvStatistics(struct FmmvStatistics *stat);"""

	print """	
struct FmmvOptions {
	double beta; /* default: 0 */

	int pM; /* default:6 */
	int pL; /* default:6 */
	int s; 

	int ws; /* default: 1 */
	int reducedScheme;

	double scale; /* default: 1 */
	int splitThreshold;  
	int splitTargetThreshold; 
	int levels;
	int directEvalThreshold;  
	int periodicBoundaryConditions; /* default: 0 */
	int extrinsicCorrection; /* default: 0 */
	int useHilbertOrder;  	/* default: 0, i.e. use Molton order */
	int directEvalAccuracy;	/* default: 2 (double) resp. 1 (single) */
	int useFarfieldNearfieldThreads;	/* default: 0 */
	"""
	if papi:
		print "\tint PAPIeventSet;"
		print
	for data in additional_options:
		print "\t%s;" % data
	print "};"	

	print """
struct FmmvOptions fmmvGetDefaultOptions(void);
	"""

        for (driver_suf, subset) in drivers:
		(charge_kinds, pot_kinds, charge_combs, pot_combs, charge_pot_combs) = \
                  	get_subdictionaryies(subset)  


		for float_type in float_types:	
			print "void fmmv%s%s(" % (driver_suf, Aa(name[float_type]))
			print "\tunsigned int NParticles, "
			print "\t%s particles[][%i]," % (float_type, dim)
			print_charge_args(charge_kinds, float_type)
			print "\tunsigned int NTargets, "
			print "\t%s targets[][%i]," % (float_type, dim)
			print_pot_args(pot_kinds, float_type)
			print "\tstruct FmmvOptions *options, /* if NULL use default options */"
			print "\tstruct FmmvStatistics *statistics, /* can be NULL */"
                	print "\tchar **errorMessage);"
			print
	
		#if len(charge_kinds)>1:
		if len(charge_kinds)>0:
			j = 0
			for i in charge_kinds:
				n = charge_kinds[i][0].upper()
				print "#define %s %i" % (n, 2**j)
				j += 1
			print
	
		#if len(pot_kinds)>1:
		if len(pot_kinds)>0:
			j = 0
			for i in pot_kinds:
				n = pot_kinds[i][0].upper()
				print "#define %s %i" % (n, 2**j)
				j += 1
			print	

		for float_type in float_types:	
			print "void fmmv%s%s_initialize(" % (driver_suf, Aa(name[float_type]))
			print "\tvoid** fmmvHandle,"
			print "\tunsigned int NParticles,"
			print "\t%s particles[][%i]," % (float_type, dim)
			if len(charge_kinds)>1:
				print "\tint typeSources,"
			print "\tunsigned int NTargets, "
			print "\t%s targets[][%i]," % (float_type, dim)
			if len(pot_kinds)>1:
				print "\tint typeTargets,"
			print "\tstruct FmmvOptions *options, /* if NULL use default options */"
			print "\tstruct FmmvStatistics *statistics, /* can be NULL */"
                	print "\tchar **errorMessage);"
			print


		for float_type in float_types:	
			print "void fmmv%s%s_evaluate(" % (driver_suf, Aa(name[float_type]))
			print "\tvoid* fmmvHandle,"
			print_charge_args(charge_kinds, float_type)
			print_pot_args(pot_kinds, float_type)
			print "\tstruct FmmvStatistics *statistics, /* can be NULL */"
                	print "\tchar **errorMessage);"
			print

		for float_type in float_types:
			print "void fmmv%s%s_finalize(" % (driver_suf, Aa(name[float_type]))
			print "\tvoid* fmmvHandle,"
			print "\tstruct FmmvStatistics *statistics, /* can be NULL */"
                	print "\tchar **errorMessage);"
			print

		#print
		#print "/* Note: fmmv%s_direct is optimized for speed," % Aa(name["generic"])
		#print "*        it uses the same kernel routines as fmmv%s" % Aa(name["generic"])
		#print "*        for direct evaluations."
		#print "*"
		#print "*        fmmv%s_direct_xp is optimized for accuracy, " % Aa(name["generic"])
		#print "*        it uses extended precision ('double-double') "
		#print "*        for accumulating potentials."
		#print "*/"
		print

		#for suf in ["", "_xp"]:
		for suf in [""]:
			for float_type in float_types:	
				print "void fmmv%s%s_direct%s(" % (driver_suf, Aa(name[float_type]), suf)
				print "\tunsigned int NParticles, "
				print "\t%s particles[][%i]," % (float_type, dim)
				print_charge_args(charge_kinds, float_type)
				print "\tunsigned int NTargets, "
				print "\t%s targets[][%i]," % (float_type, dim)
				print_pot_args(pot_kinds, float_type)
				if suf=="":
					print "\tint accuracy,"
				print "\tdouble beta,"
				print "\tdouble *time,"
                		print "\tchar **errorMessage);"
				print

	print "#endif /* _FMMV_%s_H_ */" % name["generic"].upper()


def gen_fmmv_all_in_one():
	for (driver_suf, subset) in drivers:
		(charge_kinds, pot_kinds, charge_combs, pot_combs, charge_pot_combs) = \
			get_subdictionaryies(subset)  

		print "void"
		i = 0 ; print "#if (FMM_PRECISION==%i)" % i	
		for float_type in float_types:
			if i>0: print "#elif (FMM_PRECISION==%i)" % i
			i += 1
			print "fmmv%s%s" % (driver_suf, Aa(name[float_type]))
		print "#endif"
		print "(\tunsigned int NParticles, "
		print "\t%s particles[][%i]," % ("_FLOAT_", dim)
		print_charge_args(charge_kinds, "_FLOAT_")
		print "\tunsigned int NTargets, "
		print "\t%s targets[][%i]," % ("_FLOAT_", dim)
		print_pot_args(pot_kinds, "_FLOAT_")
		print "\tstruct FmmvOptions *options, /* if NULL use default options */"
		print "\tstruct FmmvStatistics *statistics, /* can be NULL */"
        	print "\tchar **errorMessage)"
		print "{"
		print "\tFmmvHandle *FMMV;"
		print "\tchar *err;"
		if len(charge_kinds)>1:
	    		print "\tint typeSources = 0;"
		if len(pot_kinds)>1:
	    		print "\tint typeTargets = 0;"
		print
		if len(charge_kinds)>1:
			for i in charge_kinds:
				n = charge_kinds[i][0]
				print "\tif (%s) {" % n
				print "\t\ttypeSources |= %s;" % n.upper() 
				print "\t}"
	    	print
		if len(pot_kinds)>1:
	                for i in pot_kinds:
				n = pot_kinds[i][0]
				print "\tif (%s) {" % n
				print "\t\ttypeTargets |= %s;" % n.upper() 
				print "\t}"
		print
	
		i = 0 ; print "#if (FMM_PRECISION==%i)" % i
		for float_type in float_types:
			if i>0: print "#elif (FMM_PRECISION==%i)" % i
			i += 1
			print "\tfmmv%s%s_initialize" % (driver_suf, Aa(name[float_type]))
		print "#endif"
		print "\t(\t(void**) &FMMV,"
		print "\t\tNParticles, particles,", 
		if len(charge_kinds)>1:
			print "typeSources," 
		print "\t\tNTargets, targets,",
		if len(pot_kinds)>1:
			print "typeTargets," 
		print "\t\toptions, statistics, &err);"
		print "\tif (err) goto _err;"
		print 
		i = 0 ; print "#if (FMM_PRECISION==%i)" % i
		for float_type in float_types:
			if i>0: print "#elif (FMM_PRECISION==%i)" % i
			i += 1
			print "\tfmmv%s%s_evaluate" % (driver_suf, Aa(name[float_type]))
		print "#endif"
		print "\t(\t(void*) FMMV,"
		print "\t\t",
		charge_kinds_sorted = range(len(charge_kinds))
		for i in charge_kinds:
			charge_kinds_sorted[charge_kinds[i][2]] = i
		for i in charge_kinds_sorted:
			n = charge_kinds[i][0]
			print "%s," % n,
		print "\n\t\t",
		pot_kinds_sorted = range(len(pot_kinds))
		for i in pot_kinds:
			pot_kinds_sorted[pot_kinds[i][2]] = i
		for i in pot_kinds_sorted:
			n = pot_kinds[i][0]
			print "%s," % n,
		print "\n\t\t statistics, &err);"	
		print "\tif (err) goto _err;"
		print 
		i = 0 ; print "#if (FMM_PRECISION==%i)" % i
		for float_type in float_types:
			if i>0: print "#elif (FMM_PRECISION==%i)" % i
			i += 1
			print "\tfmmv%s%s_finalize((void*) FMMV, statistics, &err);" % (driver_suf, Aa(name[float_type]))
		print "#endif"
		print "\tif (err) goto _err;"
		print
		print "\t*errorMessage = 0;"
		print "_err:"
		print "\t*errorMessage = err;"
		print "}"
		print


def gen_fmmv_initialize():
	for (driver_suf, subset) in drivers:
		(charge_kinds, pot_kinds, charge_combs, pot_combs, charge_pot_combs) = \
			get_subdictionaryies(subset)  
		print "void"
		i = 0 ; print "#if (FMM_PRECISION==%i)" % i
		for float_type in float_types:
			if i>0: print "#elif (FMM_PRECISION==%i)" % i
			i += 1
			print "fmmv%s%s_initialize" % (driver_suf, Aa(name[float_type]))
		print "#endif"
		print "(\tvoid** fh,"
		print "\tunsigned int NParticles,"
		print "\t%s particles[][%i]," % ("_FLOAT_", dim)
		if len(charge_kinds)>1:
			print "\tint typeSources,"
		print "\tunsigned int NTargets, "
		print "\t%s targets[][%i]," % ("_FLOAT_", dim)
		if len(pot_kinds)>1:
			print "\tint typeTargets,"
		print "\tstruct FmmvOptions *options, /* if NULL use default options */"
		print "\tstruct FmmvStatistics *statistics, /* can be NULL */"
		print "\tchar **errorMessage)"
		print "{"
		print "\tFmmvHandle *FMMV;"
		print "\tchar *err;"
		print 
		# TODO: test magic number in evaluate, etc.
		print """
		FMMV = (FmmvHandle*) malloc(sizeof(FmmvHandle));
		FMMV->magicNumber = 1079252849; /* TODO: fmmv kind specific... */
		FMMV->NParticles = NParticles;
		FMMV->particles = particles;
		FMMV->targets = targets;
		if (targets) {	
			FMMV->NTargets = NTargets;
		}
		else {
			FMMV->NTargets = NParticles;
		}	
		"""
		if len(charge_kinds)==1:
			n = charge_kinds[charge_kinds.keys()[0]][0]
			print "\tint typeSources = %s;" % n.upper() 
		if len(pot_kinds)==1:
			n = pot_kinds[pot_kinds.keys()[0]][0]
			print "\tint typeTargets = %s;" % n.upper() 
		print 
		charge_pot_combs_sorted = range(len(charge_pot_combs))
		for i in charge_pot_combs:
			charge_pot_combs_sorted[charge_pot_combs[i][2]] = i
	        for st in ["_st", ""]:                        
			if st=="_st":
				print "\tif (FMMV->targets) {"
			else:
				print "\telse{"	
			first = True
			for i in charge_pot_combs_sorted:
				charge_comb = charge_pot_combs[i][0]
				pot_comb = charge_pot_combs[i][1]
				if first:
					print "\t\tif(",
					first = False
				else:
					print "\t\telse if(",
				for j in charge_combs[charge_comb]:
					n = charge_kinds[j][0].upper()
					print "(typeSources&%s) &&" % n,
				for j in pot_combs[pot_comb]:
					n = pot_kinds[j][0].upper()
					print "(typeTargets&%s) &&" % n,
				print "1) {"
				print "\t\t\tFMMV->dataKind = FMM%s;" % (st+i).upper()
				print "\t\t}"
			print "\t\telse {"
			print "\t\t\tassert(0);" # TODO error handling, fall backs
			print "\t\t}"
			print "\t}"	
		print
        	print "\tFMMV->beta = options->beta; "
		print "\terr = fmmv_initialize(&FMMV, options, statistics);"
		print "\tif (err) goto _err;"
		print
        	print "\tswitch(FMMV->dataKind) {"
		for st in ["", "_st"]:
			for i in charge_pot_combs: 	
	       			print "\tcase FMM%s:" % (st+i).upper()
				print "\t\tFMMV->gen_M = gen_M%s;" % (st.upper()+i)
				print "\t\tFMMV->eval_L = eval_L%s;" % (st.upper()+i)
				if st=="_st":
					print "\t\tFMMV->eval_M = eval_M%s;" % (st.upper()+i)
					print "\t\tFMMV->gen_L = gen_L%s;" % (st.upper()+i)
					print "\t\tFMMV->gen_L_eval_M = 0;" 
				else:
					print "\t\tFMMV->eval_M = 0;"
					print "\t\tFMMV->gen_L = 0;"
					print "\t\tFMMV->gen_L_eval_M = gen_L_eval_M%s;" % (st.upper()+i)
				print "\t\tFMMV->extrinsic_correction = extrinsic_correction%s;" % (st.upper()+i)
				print "\t\tswitch(FMMV->directEvalAccuracy) {"
				for j in range(eval_direct_number_of_accuracies["generic"]):
					print "\t\tcase %i:" %j
					print "\t\t\tFMMV->eval_direct = eval_direct%s_acc%i;" % (st.upper()+i, j)
					print "\t\t\tFMMV->eval_direct_periodic = eval_direct_periodic%s_acc%i;" % (st.upper()+i, j)
					print "\t\t\tbreak;"
				print "\t\t}"	
				print "\t\tbreak;"
        	print "\t}"		
		print 
		print "\t*fh = (void*) FMMV;"
		print "\t*errorMessage = 0;"
		print "_err:"
		print "\t*errorMessage = err;"
		print "}"
		print
	
def gen_fmmv_evaluate():
        for (driver_suf, subset) in drivers:
		(charge_kinds, pot_kinds, charge_combs, pot_combs, charge_pot_combs) = \
                  	get_subdictionaryies(subset)  
		print "void"
		i = 0 ; print "#if (FMM_PRECISION==%i)" % i
		for float_type in float_types:
			if i>0: print "#elif (FMM_PRECISION==%i)" % i
			i += 1
			print "fmmv%s%s_evaluate" % (driver_suf, Aa(name[float_type])) 
		print "#endif"
		print "(\tvoid* fh,"
		print_charge_args(charge_kinds, "_FLOAT_")
		print_pot_args(pot_kinds, "_FLOAT_")
		print "\tstruct FmmvStatistics *statistics, /* can be NULL */"
		print "\tchar **errorMessage)"
		print "{"
		print "\tchar *err;"
		print "\tFmmvHandle* FMMV = (FmmvHandle*) fh;"
		print
		for i in charge_kinds:
			n = charge_kinds[i][0]
			print "\tFMMV->%s = %s;" % (n, n)
		for i in pot_kinds:
			n = pot_kinds[i][0]
			print "\tFMMV->%s = %s;" % (n, n)
		print
		print "\terr = fmmv_evaluate(FMMV, statistics);";
		print "\tif (err) goto _err;"
		print 
		print "\t*errorMessage = 0;"
		print "_err:"
		print "\t*errorMessage = err;"
		print "}"
		print

def gen_fmmv_finalize():
        for (driver_suf, subset) in drivers:
		(charge_kinds, pot_kinds, charge_combs, pot_combs, charge_pot_combs) = \
                  	get_subdictionaryies(subset)  
		print "void"
		i = 0 ; print "#if (FMM_PRECISION==%i)" % i
		for float_type in float_types:
			if i>0: print "#elif (FMM_PRECISION==%i)" % i
			i += 1
			print "fmmv%s%s_finalize" % (driver_suf, Aa(name[float_type]))
		print "#endif"
		print "(\tvoid* fh,"
		print "\tstruct FmmvStatistics *statistics, /* can be NULL */"
		print "\tchar **errorMessage)"
		print "{"
		print "\tFmmvHandle* FMMV = (FmmvHandle*) fh;"
		print "\tchar *err;"
		print
		print "\terr = fmmv_finalize(FMMV, statistics);";
		print "\tif (err) goto _err;"
		print 
		print "\t*errorMessage = 0;"
		print "_err:"
		print "\t*errorMessage = err;"
		print "}"
		print

def gen_direct_method():
        for (driver_suf, subset) in drivers:
		(charge_kinds, pot_kinds, charge_combs, pot_combs, charge_pot_combs) = \
                  	get_subdictionaryies(subset)  
		print "void"
		i = 0 ; print "#if (FMM_PRECISION==%i)" % i
		for float_type in float_types:
			if i>0: print "#elif (FMM_PRECISION==%i)" % i
			i += 1
			print "fmmv%s%s_direct" % (driver_suf, Aa(name[float_type]))
		print "#endif"
		print "(\tunsigned int NParticles, "
		print "\t%s particles[][%i]," % ("_FLOAT_", dim)
		print_charge_args(charge_kinds, "_FLOAT_")
		print "\tunsigned int NTargets, "
		print "\t%s targets[][%i]," % ("_FLOAT_", dim)
		print_pot_args(pot_kinds, "_FLOAT_")
		print "\tint accuracy,"
		print "\tdouble beta,"
		print "\tdouble *time,"
		print "\tchar **errorMessage)"
		print "{"
		print """
		FmmvHandle _FMMV;
		FmmvHandle *FMMV = &_FMMV;
		Box _box;
		Box *box = &_box;
		int i;
		int err;
		int typeSources = 0;
		int typeTargets = 0;
        	"""
        	print "\tif (accuracy<0) {"
		print "\t\tdirect_method%s_xp(NParticles, particles," % driver_suf
		print_charge_args(charge_kinds, None)
		print "\t\tNTargets, targets,"
		print_pot_args(pot_kinds, None)
		print "\t\tbeta,"
		print "\t\ttime, errorMessage);"
        	print "\t\treturn;"
        	print "\t}"
        	print """
        	FMMV->NParticles = NParticles;
        	FMMV->particles = particles;
		FMMV->targets = targets;
		if (targets||((!targets)&&(NTargets>0)&&(NTargets<=NParticles))) {	
			FMMV->NTargets = NTargets;
		}
		else {
			FMMV->NTargets = NParticles;
		}	
		FMMV->directEvalAccuracy = accuracy;
		"""
		for i in charge_kinds:
			n = charge_kinds[i][0]
			print "\tFMMV->%s = %s;" % (n, n)
		for i in pot_kinds:
			n = pot_kinds[i][0]
			print "\tFMMV->%s = %s;" % (n, n)
		print	
	
		for i in charge_kinds:
			n = charge_kinds[i][0]
			print "\tif (%s) {" % n
			print "\t\ttypeSources |= %s;" % n.upper() 
			print "\t}"
		print
		for i in pot_kinds:
			n = pot_kinds[i][0]
			print "\tif (%s) {" % n
			print "\t\ttypeTargets |= %s;" % n.upper() 
			print "\t}"
		print
		charge_pot_combs_sorted = range(len(charge_pot_combs))
		for i in charge_pot_combs:
			charge_pot_combs_sorted[charge_pot_combs[i][2]] = i
		for st in ["_st", ""]:
			if st=="_st":
				print "\tif (FMMV->targets) {"
			else:
				print "\telse{"	
			first = True
			for i in charge_pot_combs_sorted:
				charge_comb = charge_pot_combs[i][0]
				pot_comb = charge_pot_combs[i][1]
				if first:
					print "\t\tif(",
					first = False
				else:
					print "\t\telse if(",
				for j in charge_combs[charge_comb]:
					n = charge_kinds[j][0].upper()
					print "(typeSources&%s) &&" % n,
				for j in pot_combs[pot_comb]:
					n = pot_kinds[j][0].upper()
					print "(typeTargets&%s) &&" % n,
				print "1) {"
				print "\t\t\tFMMV->dataKind = FMM%s;" % (st+i).upper()
				print "\t\t}"
			print "\t\telse {"
			print "\t\t\tassert(0);" # TODO error handling, fall backs
			print "\t\t}"
			print "\t}"	
		print
		print """
		FMMV->perm = (int *) malloc(NParticles*sizeof(int));
		if (FMMV->perm==0) goto _err;
		
		for (i=0; i<FMMV->NParticles; i++) {
       			FMMV->perm[i] = i;
		}
		if (targets) {
                	FMMV->permTargets = (int *) malloc(NTargets*sizeof(int));
                	if (FMMV->permTargets==0) goto _err;
			for (i=0; i<FMMV->NTargets; i++) {
        			FMMV->permTargets[i] = i;
			}
		}
		else  {
              		FMMV->permTargets = FMMV->perm;
		}
	
        	err = ida_allocate(FMMV);
        	if (err) goto _err;
	
        	copy_particles(FMMV);
        	copy_charges(FMMV);
        	zero_pot(FMMV);
		
		box->firstParticle = 0;
		box->noOfParticles = FMMV->NParticles;
		box->firstTarget = 0;
		box->noOfTargets = FMMV->NTargets;
		"""	
        	print "\tFMMV->beta = beta;"
        	print "\tswitch(FMMV->dataKind) {"
		for st in ["", "_st"]:
			for i in charge_pot_combs: 	
				print "\tcase FMM%s:" % (st+i).upper()
				print "\t\tswitch(FMMV->directEvalAccuracy) {"
				for j in range(eval_direct_number_of_accuracies["generic"]):
					print "\t\tcase %i:" % j
					print "\t\t\teval_direct%s_acc%s(FMMV, box, box);" % (st.upper()+i, j)
					print "\t\t\tbreak;"
				print "\t\t}"	
				print "\t\tbreak;"
   		print "\t}"		
		print
        	print "\tbackcopy_pot(FMMV);"
        	print "\tida_free(FMMV);"
		print
		print "\t*errorMessage = 0;"
		print "_err:"
		print "\t*errorMessage = err;"
		print "}"
		print



	
import sys

if "header" in sys.argv[1]:
	print "/* This file is automatically generated by gen_frontend.py */"
	print "/* DO NOT EDIT! */"
        print
	gen_header()
elif "source" in sys.argv[1]:
	print "/* This file is automatically generated by gen_frontend.py */"
	print "/* DO NOT EDIT! */"
        print
	print '#include"_fmmv.h"'
	print '#include<assert.h>'
	print
	gen_fmmv_all_in_one()
	gen_fmmv_initialize()
	gen_fmmv_evaluate()
	gen_fmmv_finalize()
	gen_direct_method()
