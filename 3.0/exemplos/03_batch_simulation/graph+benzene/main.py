from slabCalc.simulation import *
from slabCalc import *
import qe.pw as pw
import os

PREFIX = "graph+benzene"

# BUILDING ROUTINE GLOBAL VARIABLES
BUILD = True
SURFACES_DIR  = os.path.join("..","..","surfaces","")
MOLECULES_DIR = os.path.join("..","..","molecules","")

HC_DIR = os.path.join("hexagon_center","")
AA_DIR = os.path.join("above_atom","")
AB_DIR = os.path.join("above_bond","")

# CALCULATION GLOBAL VARIABLES
NP = 128 # number of processors
CALCULATE = True
PSEUDO_DIR  = os.path.join("..","..","pseudopotentials","")
CALC_FROM_INPUT = False
INPUT_MODEL = "input_model"
CMD = "mpirun -np %d pw.x"%NP
SAVEOUT = True
OUTFILE = "calc.out"
SAVECOORDS = True
COORDSFILE = "final_coords.xyz"

a = 2.46 # surface cell parameter in angstroms
n,m = 3,2 # surface repetition numbers

def create(fname):
	"""
	creates a mol3D structure using the file name passed in fname
	"""
	structure = mol3D()
	structure.readfromxyz(fname)
	return structure
	
def make_dirs():
	dirs = ["hexagon_center", "above_atom", "above_bond"]
	
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
		calc.control["prefix"]        = "\"batch_%s\""%PREFIX
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
	surface_file = "graphene_%dx%d.xyz"%(n,m)
	molecule_file = "benzene.xyz"
	surface = create(SURFACES_DIR + surface_file)
	molecule = create(MOLECULES_DIR + molecule_file)
	
	# Molecule's alocation parameters
	params = list()
	
	# Defining sites
	hc_site = surface.get_coord_between(4,11) #hexagon center site
	aa_site = surface.get_coord_between(11) #above atom site
	ab_site = surface.get_coord_between(11,12) #above bond site
	sites = [hc_site,aa_site,ab_site]
	
	# Defining align distances
	angles = [0, 15, 30, 45]
	
	# Filling the params list
	for s in sites:
		for a in angles:
			# Label slabs according to site
			if s==hc_site:
				label = "hexagon_center"
			elif s==aa_site:
				label = "above_atom"
			elif s==ab_site:
				label = "above_bond"
		
			iparam = {    "surface" : surface,
						"molecules" : [molecule],
							 "site" : s,
							"label" : label,
				   	  "align_point" : molecule.centersym(), 
						     "dist" : 3,
						 "rotation" : True,
						    "angle" : a
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
		slabpath = "%d.xyz"%(slab.molecules[0].angle)
		if slab.molecules[0].site==hc_site:
			slabpath = HC_DIR + slabpath
		elif slab.molecules[0].site==aa_site:
			slabpath = AA_DIR + slabpath
		elif slab.molecules[0].site==ab_site:
			slabpath = AB_DIR + slabpath
		slab.writexyz(slabpath)
		
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
		sim.run_qe(cmd=CMD)
		
	sim.save()
	
	print("Done!")
