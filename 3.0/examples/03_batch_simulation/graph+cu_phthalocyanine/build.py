from slabCalc.simulation import *
from slabCalc import *
import os

SURFACES_DIR  = os.path.join("..","surfaces","")
MOLECULES_DIR = os.path.join("..","molecules","")

HC_DIR = os.path.join("hexagon_center","")
AA_DIR = os.path.join("above_atom","")
AB_DIR = os.path.join("above_bond","")

def create(fname):
	"""
	creates a mol3D structure using the file name passed in fname
	"""
	structure = mol3D()
	structure.readfromxyz(fname)
	return structure
	
def make_dirs():
	dirs = ["hexagon_center","above_atom","above_bond"]
	
	for i in dirs:
		if i not in os.listdir("."):
			os.mkdir(i)

if __name__=="__main__":
	make_dirs()

	"""SETTING BUILDING OPTIONS"""
	
	# Selecting a surface and molecules
	surface_file = "graphene_7x4.xyz"
	molecule_file = "cu_phthalocyanine.xyz"
	surface = create(SURFACES_DIR + surface_file)
	molecule = create(MOLECULES_DIR + molecule_file)
	
	# Molecule's alocation parameters
	
	params = list()
	
	# Defining sites
	hc_site = surface.get_coord_between(40,55) #hexagon center site
	aa_site = surface.get_coord_between(55) #above atom site
	ab_site = surface.get_coord_between(55,56) #above bond site
	sites = [hc_site,aa_site,ab_site]
	
	# Defining rotational parameters
	angles = [0, 15, 30, 45]
	
	# Filling the params list
	for s in sites:
		for a in angles:
			iparam = {    "surface" : surface,
						"molecules" : molecule,
							 "site" : s,
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
	
	# Writing resulting structure to xyz file
	for slab in sim.slabs:
		slabpath = "%.0f.xyz"%(slab.molecules[0].angle)
		if slab.molecules[0].site==hc_site:
			slabpath = HC_DIR + slabpath
		elif slab.molecules[0].site==aa_site:
			slabpath = AA_DIR + slabpath
		elif slab.molecules[0].site==ab_site:
			slabpath = AB_DIR + slabpath
		slab.writexyz(slabpath)
	
	print("Done!")
