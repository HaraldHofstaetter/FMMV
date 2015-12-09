dim = 3

dir = "3d"

float_types = ["float" , "double"]

type_name = {
        "float": "float",
        "double": "double",
        }
        
name = {
	"generic":   "3d",
	"float":  "3df",
	"double": "3d",
	}
        
simd_vector_length = {
	"float":  4, # for sse
	"double": 2, # for sse2
	# "double": 4, # for avx
	}
        
p_max = { 
	"float":  20, 
	"double": 30,
	}

s_eps_max = {	
	"float":  20, 
	"double": 30, 
	}

s_exp_max = {	
	"float":  1200, 
	"double": 2400, 
	}

eval_direct_number_of_accuracies = {
	"generic": 3, 
	"float":  2, 
	"double": 3, 
	}

#####	

charge_kinds = { 
	"_charge": ( "charges",       ["q"],              0 ),
	"_dipole": ( "dipoleMoments", ["mx", "my", "mz"], 1 ),
	}

pot_kinds = {
	"_pot":    ("potentials",     ["pot"],                     0 ),
	"_grad":   ("gradients",      ["gradx", "grady", "gradz"], 1 ),
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

drivers = [
        ("",  ["_standard", "_dipole", "_grad", "_dipole_grad"]),
        ]

########

ida= { 
    "float": {
	"_standard":    [['x', 'y', 'z', 'q', 'pot']],
	"_dipole":	[['x', 'y', 'z', 'q', 'mx', 'my', 'mz', 'pot']],
	"_grad":	[['x', 'y', 'z', 'q', 'pot', 'gradx', 'grady', 'gradz']],
	"_dipole_grad":	[['x', 'y', 'z', 'q', 'mx', 'my', 'mz', 'pot', 'gradx', 'grady', 'gradz']],
    },
    "double": {
	
	"_standard":    [['x', 'y', 'z', 'q', 'pot']],
	"_dipole":	[['x', 'y', 'z', 'q', 'mx', 'my', 'mz', 'pot']],
	"_grad":	[['x', 'y', 'z', 'q', 'pot', 'gradx', 'grady', 'gradz']],
	"_dipole_grad":	[['x', 'y', 'z', 'q', 'mx', 'my', 'mz', 'pot', 'gradx', 'grady', 'gradz']],
    },
}
      
ida_st = {	
    "float": {
	"_standard":    [['x', 'y', 'z', 'q'], ['tx', 'ty', 'tz', 'pot']],
	"_dipole":	[['x', 'y', 'z', 'q', 'mx', 'my', 'mz'], ['tx', 'ty', 'tz', 'pot']],
	"_grad":	[['x', 'y', 'z', 'q'], ['tx', 'ty', 'tz', 'pot', 'gradx', 'grady', 'gradz']],
	"_dipole_grad": [['x', 'y', 'z', 'q', 'mx', 'my', 'mz'], ['tx', 'ty', 'tz', 'pot', 'gradx', 'grady', 'gradz']],
    },
    "double": {
	"_standard":    [['x', 'y', 'z', 'q'], ['tx', 'ty', 'tz', 'pot']],
	"_dipole":	[['x', 'y', 'z', 'q', 'mx', 'my', 'mz'], ['tx', 'ty', 'tz', 'pot']],
	"_grad":	[['x', 'y', 'z', 'q'], ['tx', 'ty', 'tz', 'pot', 'gradx', 'grady', 'gradz']],
	"_dipole_grad": [['x', 'y', 'z', 'q', 'mx', 'my', 'mz'], ['tx', 'ty', 'tz', 'pot', 'gradx', 'grady', 'gradz']],
    },
}

ida_nosimd = ida
ida_st_nosimd = ida_st

######

additional_options = [
	# user defined quadrature coeffs: 
	"double *x",
	"double *w",
	"int *M",
]

additional_statistics = []

additional_data_in_box = [] # no additional data

additional_global_data = [
	"_FLOAT_ beta",

	"_FLOAT_ *lambda",
	"_FLOAT_ *w",
	"int *M",
	"int *FFT_rep",
	"int len_F",

	"_FLOAT_ **Tz_M2M",
	"_FLOAT_ **Tz_L2L",
	"_FLOAT_ *Tz_M2M_coeffs",
	"_FLOAT_ *Tz_L2L_coeffs",

	"_FLOAT_ *CMX_coeffs",
	"_FLOAT_ *CXL_coeffs",
	"_FLOAT_ *CXL_reduced_coeffs",
	
	"_FLOAT_ **Ry_theta",
	"_FLOAT_ **Ry_minus_theta",
	"_FLOAT_ **Ry_pi_minus_theta",
	"_FLOAT_ **Ry_minus_pi_minus_theta",
	"_FLOAT_ **Ry_pi2",
	"_FLOAT_ **Ry_minus_pi2",
	"_FLOAT_ *Ry_theta_coeffs",
	"_FLOAT_ *Ry_minus_theta_coeffs",
	"_FLOAT_ *Ry_pi_minus_theta_coeffs",
	"_FLOAT_ *Ry_minus_pi_minus_theta_coeffs",
	"_FLOAT_ *Ry_pi2_coeffs",
	"_FLOAT_ *Ry_minus_pi2_coeffs",

	"int *P_MRT",
	"int *P_MRT_plus",
	"int *P_MRT_minus",
	"int *P_Mriri2rrii",
	"int *P_LRT",
	"int *P_LRT_plus",
	"int *P_LRT_minus",
	"int *P_Lriri2rrii",
	"int *P_X_riri2rrii",
	"int *P_VF",
	"int *neg_F",
	"int len_neg_F",
	"int *P_FV",
	"int *neg_V",
	"int len_neg_V",
]






     
