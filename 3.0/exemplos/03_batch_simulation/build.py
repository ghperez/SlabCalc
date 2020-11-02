from slabCalc.simulation import *
from slabCalc import *
import os

SURFACES_DIR  = os.path.join("..","surfaces","")
MOLECULES_DIR = os.path.join("..","molecules","")

STRUC_DIR = os.path.join("structures","")

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

if __name__=="__main__":
	make_dirs()

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
					     "dist" : d
				 }
		
		params.append(iparam)
		      
	"""BUILDING"""
	
	# Creating Sim() object
	sim = Sim()
	
	# Building
	sim.create_slabs(params)
	
	sim.build_slabs(silent=True)
	
	# Writing resulting structure to xyz file
	for slab in sim.slabs:
		slab.writexyz(STRUC_DIR+"%.0f.xyz"%slab.molecules[0].dist)
	
	print("Done!")
