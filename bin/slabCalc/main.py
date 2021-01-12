"""
    SlabCalc
    --------
	Set of objects for studying 2D atomic surfaces with python

author:       Gabriel
institution:  Universidade Federal do ABC (UFABC)
"""

from molSimplify.Scripts.cellbuilder_tools import*
from molSimplify.Scripts.cellbuilder import*
import molSimplify.Classes.mol3D as ms
from copy import deepcopy
import qe.pw as pw
import os

class mol3D(ms.mol3D):
	"""
	Modification of Molsimplify's mol3D object with new implemented and superposed methods
	"""
	def __init__(self):
		super(mol3D, self).__init__()
		self.site = list()
		self.align_point = list()
		self.dist = None
		self.rotation = None
		self.angle = None
		self.axis = list()
		self.rpoint = list()
		
	def set_molecule_parameters(self,params):
		"""
		Set the alocation parameters

		Arguments
        ---------
        (list) site         : coordinations of surface's site
        (list) align_point  : coordinations of the molecule's align point
        (float) dist        : align distance
		(bool) rotation     : If true then rotate the molecule
			
		Rotational args:
		---------------
		(float) angle        
        (list) axis         : z axis is the default
        (list) rpoint       : intercept point of the axis
		"""
		self.site   = params["site"].copy()
		self.dist   = params["dist"]
		self.align_point = params["align_point"].copy()
		
		if "rotation" in params.keys() and params["rotation"]==True:
			self.rotation = True
			self.angle    = params["angle"]
			if "axis" in params.keys():
				self.axis = params["axis"].copy()
			else:
				self.axis = [0,0,1]
			if "rpoint" in params.keys():
				self.rpoint = params["rpoint"].copy()
			else:
				self.rpoint = self.centersym()
		
	def get_coord_between(self,*inds):
		if type(inds[0])==list or type(inds[0])==tuple:
			inds = inds[0]
		return center_of_sym([self.getAtom(i).coords() for i in inds])
		
	def copymol3D(self,other,copymols=False):
		super(mol3D, self).copymol3D(other)
		if copymols:
			self.copy_molecule_parameters(other)
		
	def copy_molecule_parameters(self,other):
		params = {       "site" : other.site.copy(),
		   		  "align_point" : other.align_point.copy(),
		   		  		 "dist" : other.dist,
		   		     "rotation" : other.rotation,
		   		        "angle" : other.angle,
		   		         "axis" : other.axis.copy(),
		   		       "rpoint" : other.rpoint.copy()
				 }
				 
		self.set_molecule_parameters(params)

class Slab(mol3D,object):
	"""
	Surface's class

	mol3D's child class, for more information please read Molsimplify's documentation
	"""
	def __init__(self,surface=None,*molecules,cell_vector=None,extents=None):
		"""
		Slab's properties:
		-----------------
			(object)  surface       :  mol3D's surface instance
		    (list)  molecules       :  list of mol3D's molecules instances
			(list)  cell_vector     :  unity cell vector used for establish surface's extention.
									   The user may pass this parameter or the 'extents' one.
			(list)  extents         :  vector of surface's boarders
		"""
		super(Slab, self).__init__()
		
		self.molecules = list(molecules)
		if len(self.molecules)!=0:
			while type(self.molecules[0])!=mol3D:
				self.molecules = self.molecules[0]
		
		if surface==None:
			self.surface = mol3D()
		else:
			self.surface = surface

		if cell_vector!=None:
			self.extents = find_extents_cv(cell_vector)
		elif extents!=None:
			self.extents = extents
		else:
			self.extents = [0,0,0]

	def build(self,silent=False):
		"""
		Insert the molecule above the surface
            
            (list) msite : align point in molecule
           (float) adist : align distance
		"""
		# INITIALISING THE COMBINED SLAB
		cslab = mol3D()
		if self.surface!=None:
			cslab.copymol3D(self.surface)
		else:
			print("\n!!! Please set the slab's surface !!!\n")
			
		for i in range(len(self.molecules)):
			# INITIALISING THE PAYLOAD MOLECULE
			payload = mol3D()
			payload.copymol3D(self.molecules[i],copymols=True)
			# ROTATION
			if payload.rotation==True and payload.angle!=0:
				try:
					rotate_around_axis(payload,payload.rpoint,payload.axis,payload.angle)
					if not silent :
						print("Molecule Rotated")
				except:
					print("\n!!! Error during molecule%i rotation !!!\nPlease, check the molecule rotational parameters\n"%(i+1))
			# TRANSLATION
			try:
				trans_vec = vecdiff(payload.site,payload.align_point)
				payload.translate(trans_vec)
				payload.translate([0, 0, self.extents[2]+payload.dist])
				# COMBINING MOLECULE AND SURFACE INTO SLAB
				cslab.combine(payload)
			except:
				print("\n!!! Error while inserting molecule %d !!!\nPlease, check the molecule parameters.\n"%(i+1))
			# CHECKING SANITY
			self.sanity = self.sanitycheck(silence=True)
			if self.sanity[0] and not silent:
				print("!!! Molecule %i might be overlapping!\nThe minimum distance is %f.\n"%(i+1,self.sanity[1]))
		
		self.copymol3D(cslab)
		if not silent:
			print("Slab built")

	def set_qe_calculation(self,in_calc):
		"""
		Returns a Quantum Espresso's input as a string according to a calc object pre-set
        
        (object) calc      : QE's library object for setting Quantum Espresso's calculations
		"""
		out_calc = deepcopy(in_calc)
		
		for atom in self.atoms:
			out_calc.x.append(atom.coords()[0])
			out_calc.y.append(atom.coords()[1])
			out_calc.z.append(atom.coords()[2])
			out_calc.atom_type.append(atom.symbol())
			if atom.symbol() not in out_calc.atomic_species:
				out_calc.atomic_mass.append(atom.mass)
				out_calc.atomic_species.append(atom.symbol())

		out_calc.system["nat"] = len(self.atoms)
		out_calc.system["ntyp"] = len(out_calc.atomic_species)
		
		return out_calc
    
	def copy(self):
		copyslab = Slab()
		copyslab.surface.copymol3D(self.surface)
		copyslab.molecules = self.copy_molecules()
		copyslab.extents = self.extents.copy()
		copyslab.copymol3D(self)
		
		return copyslab
		
	def copy_molecules(self):
		"""
		Return a copy of molecules list
		"""
		molecules_copy = list()
		if len(self.molecules)!=0:
			for m in self.molecules:
				copymol = mol3D()
				copymol.copymol3D(m,copymols=True)
				molecules_copy.append(copymol)
				
		return molecules_copy

	def clear(self):
		"""
		Resets the slab to the pre-built form
		"""
		self.__init__(self.surface,self.molecules,extents=self.extents)
		print("Slab reset!")
