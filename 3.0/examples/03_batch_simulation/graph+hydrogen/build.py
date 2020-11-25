from slabCalc.simulation import *
from slabCalc import *
import qe.pw as pw
import os

PREFIX = "graph+hydrogen"

# BUILDING ROUTINE GLOBAL VARIABLES
BUILD = True
SURFACES_DIR  = os.path.join("..","..","surfaces","")
MOLECULES_DIR = os.path.join("..","..","molecules","")
STRUC_DIR = os.path.join("structures","")


# CALCULATION GLOBAL VARIABLES
CALCULATE = True
PSEUDO_DIR  = os.path.join("..","pseudopotentials","")
CALC_FROM_INPUT = False
INPUT_MODEL = "input_model"
CMD = "mpirun -np 8 pw.x"
SAVEOUT = True
OUTFILE = "calc.out"
SAVECOORDS = True
COORDSFILE = "final_coords.xyz"

a = 2.46 # surface cell parameter in angstroms
n,m = 2,1 # surface repetition numbers

def create(fname):
	"""
	creates a mol3D structure using the file name passed in fname
	"""
	structure = mol3D()
	structure.readfromxyz(fname)
	return structure
	
def make_dirs():
	dirs = ["structures"]
	
	for i in dirs:
		if i not in os.listdir("."):
			os.mkdir(i)

def set_calc():
	calc = pw.calc()
	
	if CALC_FROM_INPUT:
		calc.read_input(INPUT_MODEL)
	else:
		""" SET THE CALCULATION PARAMETERS HERE """
		# &CONTROL
		calc.control["calculation"]   = "\"relax\""
		calc.control["restart_mode"]  = "\"from_scratch\""
		calc.control["pseudo_dir"]    = "\""+PSEUDO_DIR+"\""
		calc.control["prefix"]        = "\"simple_example\""
		calc.control["outdir"]        = "\"outputs\""
		calc.control["etot_conv_thr"] = "1.0E-5"
		calc.control["forc_conv_thr"] = "1.0D-4"
		calc.control["nstep"]         = 200
		
		# &SYSTEM
		calc.system["ibrav"]       = 0
		calc.system["celldm(1)"]   = a*1.889725989 # cell parameter in bohr
		#calc.system["nat"]         =
		#calc.system["ntyp"]        =
		calc.system["nspin"]       = 1
		calc.system["occupations"] = "\"smearing\""
		calc.system["degauss"]     = 0.02
		calc.system["ecutwfc"]     = 32
		calc.system["ecutrho"]     = 320.0
		calc.system["smearing"]    = "\"gaussian\""
		calc.system["input_dft"]   = "\"vdW-DF\""
		
		# &ELECTRONS
		calc.electrons["mixing_beta"]      = 0.5
		calc.electrons["electron_maxstep"] = 100
		calc.electrons["conv_thr"]         = "1.0D-6"
		calc.electrons["diagonalization"]  = "\"david\""
		calc.electrons["startingpot"]      = "\"atomic\""
		calc.electrons["startingwfc"]      = "\"atomic+random\""
		
		# &IONS
		calc.ions["ion_dynamics"] = "\"bfgs\""
		
		# CELL_PARAMETERS
		calc.cell_parameters_units = "alat"
		calc.v1 = [n,0,0]
		calc.v2 = [0,m*sqrt(3),0]
		calc.v3 = [0,0,10]
		
		# ATOMIC_SPECIES
		calc.pseudopotential = ["C.pbe-rrkjus.UPF","H.pbe-rrkjus_psl.1.0.0.UPF"]
		
		# ATOMIC_POSITIONS
		calc.atomic_positions_units = "angstrom"
		
		# K_POINTS
		calc.k_points_type = "automatic"
		calc.nk1 = 4
		calc.nk2 = 4
		calc.nk3 = 1
		calc.sk1 = 0
		calc.sk2 = 0
		calc.sk3 = 0
					
	return calc
	
def build_structures():
	"""SETTING BUILDING OPTIONS"""
	
	# Selecting a surface and molecules
	surface_file = "graphene_2x1.xyz"
	molecule_file = "hydrogen.xyz"
	surface = create(SURFACES_DIR + surface_file)
	molecule = create(MOLECULES_DIR + molecule_file)
	
	# Molecule's alocation parameters
	params = list()
	
	# Defining sites
	site = surface.get_coord_between(2,5) #hexagon center site
	
	# Defining align distances
	dists = [3,4,5]
	
	# Filling the params list
	for d in dists:
		iparam = {    "surface" : surface,
					"molecules" : [molecule],
						 "site" : site,
			   	  "align_point" : molecule.centersym(), 
					     "dist" : d,
					  "saveinp" : True,
					  "inpfile" : "slab%.0f.in"%d
				 }
		
		params.append(iparam)
		      
	"""BUILDING"""
	
	# Creating Sim() object
	sim = Sim()
	
	# Building
	sim.create_slabs(params)
	
	sim.build_slabs(silent=True)
	
	#Writing resulting structure to xyz file
	for slab in sim.slabs:
		slab.writexyz(STRUC_DIR+"%.0f.xyz"%slab.molecules[0].dist)
		
	sim.save()
	
	return sim
	
if __name__=="__main__":

	make_dirs()
	
	#Building Routine
	if BUILD:
		sim = build_structures()
	else:
		sim = Sim()
		sim.load()
	
	#Calculations
	if CALCULATE:
		calc = set_calc()
		sim.set_qe(calc)
		sim.run_qe(cmd=CMD,silent=False)
		
	sim.save()
	
	print("Done!")
