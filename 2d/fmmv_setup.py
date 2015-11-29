dim = 2

dir = "2d"

float_types= ["float" , "double"]
type_name = {
        "float": "float",
        "double": "double",
        }
name = {
	"generic":   "2d",
	"float":  "2df",
	"double": "2d",
	}
simd_vector_length = {
	"float":  1,   # for the time being...
	"double": 1,
	}
p_max = { 
	"float":  30, 
	"double": 40,
	}
s_eps_max = {	
	"float":  30, 
	"double": 40, 
	}
s_exp_max = {	
	"float":  30, # equal to s_eps_max 
	"double": 40, 
	}

eval_direct_number_of_accuracies = {
	"generic": 3, # for the time being 
	"float":  2, 
	"double": 3, 
	}

#####	

charge_kinds = { 
	"_charge": ( "charges",       ["q"],        0 ),
	"_dipole": ( "dipoleMoments", ["mx", "my"], 1 ),
	}

pot_kinds = {
	"_pot":    ("potentials",     ["pot"],            0 ),
	"_grad":   ("gradients",      ["gradx", "grady"], 1 ),
	}

charge_combs = {        
	"_charge": ["_charge"],
	"_dipole": ["_charge", "_dipole"],
	}
	
pot_combs = {        
	"_pot":    ["_pot"],
	"_grad":   ["_pot", "_grad"],
	}

charge_pot_combs = { #  (charge_comb,  pot_comb, look for match )
	"_standard":    ("_charge",    "_pot",  3 ),
	"_dipole":      ("_dipole",    "_pot",  2 ),
	"_grad":        ("_charge",    "_grad", 1 ),
	"_dipole_grad": ("_dipole",    "_grad", 0 ),
	}


########

ida= { 
    "float": {
	"_standard":    [['x', 'y', 'q', 'pot']],
	"_dipole":	[['x', 'y', 'q', 'mx', 'my', 'pot']],
	"_grad":	[['x', 'y', 'q', 'pot', 'gradx', 'grady']],
	"_dipole_grad":	[['x', 'y', 'q', 'mx', 'my', 'pot', 'gradx', 'grady']],
    },
    "double": {
	
	"_standard":    [['x', 'y', 'q', 'pot']],
	"_dipole":	[['x', 'y', 'q', 'mx', 'my', 'pot']],
	"_grad":	[['x', 'y', 'q', 'pot', 'gradx', 'grady']],
	"_dipole_grad":	[['x', 'y', 'q', 'mx', 'my', 'pot', 'gradx', 'grady']],
    },
}
      
ida_st = {	
    "float": {
	"_standard":    [['x', 'y', 'q'], ['tx', 'ty', 'pot']],
	"_dipole":	[['x', 'y', 'q', 'mx', 'my'],  ['tx', 'ty', 'pot']],
	"_grad":	[['x', 'y', 'q'], ['tx', 'ty', 'pot', 'gradx', 'grady']],
	"_dipole_grad": [['x', 'y', 'q', 'mx', 'my'],  ['tx', 'ty', 'pot', 'gradx', 'grady']],
    },
    "double": {
	"_standard":    [['x', 'y', 'q'], ['tx', 'ty', 'pot']],
	"_dipole":	[['x', 'y', 'q', 'mx', 'my'],  ['tx', 'ty', 'pot']],
	"_grad":	[['x', 'y', 'q'], ['tx', 'ty', 'pot', 'gradx', 'grady']],
	"_dipole_grad": [['x', 'y', 'q', 'mx', 'my'],  ['tx', 'ty', 'pot', 'gradx', 'grady']],
    },
}



######

additional_options = [
]

additional_statistics = []

additional_data_in_box = [] # no additional data

additional_global_data = [
	"_FLOAT_ beta",
        "_FLOAT_ *recip",
	"_FLOAT_ *lambda",
	"_FLOAT_ *w",
	"_FLOAT_ *A[FMM_P_MAX+2]",
	"_FLOAT_ *Ac[FMM_P_MAX+2]",
	"_FLOAT_ *vander_over_fact[FMM_S_EXP_MAX]",
	"_FLOAT_ C",
       	
	"_FLOAT_ k0_correction",
	"double *FF",
	"double *FF_INV",
]






     
