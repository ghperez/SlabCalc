from slabCalc import *
import os

SURFACES_DIR  = os.path.join("..","surfaces","")
MOLECULES_DIR = os.path.join("..","molecules","")
SURFACE_NAME = "graphene_2x1"
MOLECULE_NAME = "hydrogen"

def create(fname):
	"""
	creates a mol3D structure using the file name passed in fname
	"""
	structure = mol3D()
	structure.readfromxyz(fname)
	return structure

if __name__=="__main__":

	"""SETTING BUILDING OPTIONS"""
	
	# Selecting a surface and molecules
	surface = create(SURFACES_DIR + SURFACE_NAME + ".xyz")
	molecule = create(MOLECULES_DIR + MOLECULE_NAME + ".xyz")
	
	# Molecule's alocation parameters
	params = {       "site" : surface.get_coord_between(5,2), #get coordination between atoms indexes
		   	  "align_point" : molecule.centersym(), 
		             "dist" : 3, #align distance between "site" in surface and "align_point" in molecule
		         "rotation" : True,
		            "angle" : 45
		     } 
	molecule.set_molecule_parameters(params)
		      
	"""BUILDING"""
	
	# Creating slab object
	slab = Slab(surface,molecule)
	
	# Building
	slab.build()

	# Writing resulting structure to xyz file
	slab.writexyz("graphene_hydrogen.xyz")
	
	print("Done!")
