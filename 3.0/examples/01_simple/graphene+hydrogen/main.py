from slabCalc import *
import os

SURFACES_DIR  = os.path.join("..","..","surfaces","")
MOLECULES_DIR = os.path.join("..","..","molecules","")
PREFIX = "graph+hydrogen"

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
	surface_file = "graphene_2x1.xyz"
	molecule_file = "hydrogen.xyz"
	surface = create(SURFACES_DIR + surface_file)
	molecule = create(MOLECULES_DIR + molecule_file)
	
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
	slab.writexyz("%s.xyz"%PREFIX)
	
	print("Done!")
