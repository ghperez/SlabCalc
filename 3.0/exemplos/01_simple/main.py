from slabCalc import *
from math import *
import os

SURFACES_DIR  = os.path.join("..","surfaces","")
MOLECULES_DIR = os.path.join("..","molecules","")

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
	surface_file = "graphene.xyz"
	molecule_file = "benzene.xyz"
	surface = create(SURFACES_DIR + surface_file)
	molecule = create(MOLECULES_DIR + molecule_file)
	
	# Molecule's alocation parameters
	inds = (14,21)
	params = {       "site" : surface.get_coord_between(14,21), #get coordination between atoms indexes
		   	  "align_point" : molecule.centersym(), 
		             "dist" : 3, #align distance between "site" in surface and "align_point" in molecule
		         "rotation" : True,
		            "angle" : 30
		     } 
	molecule.set_molecule_parameters(params)
		      
	"""BUILDING"""
	
	# Creating slab object
	slab = Slab(surface,molecule)
	
	# Building
	slab.build()
		
	# Writing resulting structure to xyz file
	slab.writexyz("graphene_benzene.xyz")
	
	print("Done!")
