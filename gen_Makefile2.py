from fmmv_setup import charge_pot_combs, float_types, eval_direct_number_of_accuracies, dim

print "# This file is automatically generated by gen_Makefile2.py */"
print "# DO NOT EDIT! */"
print

data = [
("gen_eval_access",      False),
("eval_direct",          True),
("extrinsic_correction", False), 
]

simds = ["", "_simd2", "_simd4"]

print "NAMES_ACCESS = \\"
for (name, acc) in data:
	for i in charge_pot_combs:
	    for st in ["", "_st"]:
		if  acc:
			for j in range(eval_direct_number_of_accuracies["generic"]):	
				print "\t"+name+st+i+("_acc%i" % j)+"\\"
				print "\t"+name+"_periodic"+st+i+("_acc%i" % j)+"\\"
		else:		
			print "\t"+name+st+i+"\\"
print			
print

float_type_nr = 0
for float_type in float_types:
	for (name, acc) in data:
	     for simd in simds:
		for i in charge_pot_combs:
		    for st in ["", "_st"]:
			if  acc:
				for j in range(eval_direct_number_of_accuracies["generic"]):	
					dep = name+simd+".c"
					target = "$(ODIR)/"+name+st+i+("_acc%i" % j)+simd+"_"+float_type+".o"
					print "%s: %s" % (target, dep)
					print "\t$(LIBTOOL_COMPILE) -c -D%s -DACCURACY=%i -DFMM_PRECISION=%i $(CFLAGS) -I.. %s -o %s" % (
						(st+i)[1:].upper(), j, float_type_nr, dep, target)
					print	
					target = "$(ODIR)/"+name+"_periodic"+st+i+("_acc%i" % j)+simd+"_"+float_type+".o"
					print "%s: %s" % (target, dep)
					print "\t$(LIBTOOL_COMPILE) -c -DPERIODIC -D%s -DACCURACY=%i -DFMM_PRECISION=%i $(CFLAGS) -I.. %s -o %s" % (
						(st+i)[1:].upper(), j, float_type_nr, dep, target)
					print	
			else:
				target = "$(ODIR)/"+name+st+i+simd+"_"+float_type+".o"
				dep = name+simd+".c"
				print "%s: %s" % (target, dep)
				print "\t$(LIBTOOL_COMPILE) -c -D%s -DFMM_PRECISION=%i $(CFLAGS) -I.. %s -o %s" % (
					(st+i)[1:].upper(), float_type_nr, dep, target)
				print	
	float_type_nr += 1					

