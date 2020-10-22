#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    SlabCalc
    --------
    Conjunto de objetos para o estudo de superfícies 2D utilizando métodos
    computacionais.

autor:       Gabriel
instituição: Universidade Federal do ABC (UFABC)
"""

from molSimplify.Scripts.cellbuilder_tools import*
from molSimplify.Scripts.cellbuilder import*
import molSimplify.Classes.mol3D as ms
import qe.pw as pw
import os

class mol3D(ms.mol3D):
	"""
	Redefinicao da classe mol3D para uso no slabCalc
	"""
	def __init__(self):
		super(mol3D, self).__init__()
		self.site = None
		self.align_point = None
		self.dist = None
		self.rotation = None
		self.angle = None
		self.axis = None
		self.rpoint = None
		
	def set_molecule_parameters(self,params):
		"""
		Define os parametros para a alocacao da molecula sobre a superficie
		
		Parametros
        ----------
        	(list) site         : sitio sobre a superficie
        	(list) align_point  : ponto de alinhamento na molecula
           (float) dist         : distancia de alinhamento
            (bool) rotation     : se Verdadeiro então a molecula vai ser rotacionada
           (float) angle        : angulo
            (list) axis         : eixo (z axis is the default)
            (list) rpoint       : ponto por onde passa o eixo de rotacao
		"""
		self.site   = params["site"].copy()
		self.dist   = params["dist"]
		self.align_point = params["align_point"].copy()
		
		if "rotation" in params.keys() and params["rotation"]==True:
			try:
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
			except:
				pass
		else:
			self.rotation = False 
			self.angle    = None
			self.axis     = None
			self.rpoint   = None
		
	def get_coord_between(self,*inds):
		if type(inds[0])==list or type(inds[0])==tuple:
			inds = inds[0]
		return center_of_sym([self.getAtom(i).coords() for i in inds])
		
	def copymol3D(self,other):
		super(mol3D, self).copymol3D(other)
		try:
			self.copy_molecule_parameters(other)
		except:
			pass
		
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
	Classe da superfície

	subclasse da classe mol3D, para mais informações leia a documentação do molSimplify
	"""
	def __init__(self,surface=None,*molecules,cell_vector=None,extents=None):
		"""
		Método Construtor

		Propriedades
			(object)  surface       :  objeto mol3D da superfície
              
			(object)  molecules     :  objetos mol3D de moléculas
              
			(list)  cell_vector  :  vetores de célula, usados para determinar
                                      a extensão da superfície. Você pode tanto
                                      passar este parâmetro ou o parâmetro extents
                                      
			(list)  extents       :  um vetor de coordenadas dos limites da superfície
		"""
		super(Slab, self).__init__()
		
		if type(molecules)==tuple:
			self.molecules = list(molecules)
		
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
		Insere a molécula sob a superfície

            
            (list) msite : ponto que vai ser alinhado na molécula
           (float) adist : distância de alinhamento em angstrons

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
			payload.copymol3D(self.molecules[i])
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
				print("\n!!! Error while building molecule %d !!!\nPlease, check the molecule parameters.\n"%(i+1))
			# CHECKING SANITY
			self.sanity = self.sanitycheck(silence=True)
			if self.sanity[0] and not silent:
				print("!!! Molecule %i might be overlapping!\nThe minimum distance is %f.\n"%(i+1,self.sanity[1]))
		
		self.copymol3D(cslab)
		if not silent:
			print("Slab built")

	def set_qe_calculation(self,calc):
		"""
        Retorna um input do Quantum Espresso em formato de string confome um
        um objeto calc pré-configurado
        
        (object) calc      : objeto da biblioteca QE
		"""
		for atom in self.atoms:
			calc.x.append(atom.coords()[0])
			calc.y.append(atom.coords()[1])
			calc.z.append(atom.coords()[2])
			calc.atom_type.append(atom.symbol())
			if atom.symbol() not in calc.atomic_species:
				calc.atomic_mass.append(atom.mass)
				calc.atomic_species.append(atom.symbol())

		calc.system["nat"] = len(self.atoms)
		calc.system["ntyp"] = len(calc.atomic_species)
        
		return calc
    
	def copy(self):
		"""
        Método para fazer cópias das superfícies
		"""
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
				copymol.copymol3D(m)
				molecules_copy.append(copymol)
		return molecules_copy

	def clear(self):
		"""
        Reinicia a superfície
		"""
		self.__init__(self.surface,self.molecules,extents=self.extents)
		print("Slab reset!")

"""
    Funções úteis para obter os pontos de alocação:

        - center_of_sym([molecule.getAtom(i).coords() for i in inds])
            # retorna as coordenadas do ponto entre os átomos de índicies
            especificados em "inds"
        - mol3D_object.centersym()
            # retorna o centro de simetria de um objeto mol3D
            
	Leitura e escrita de arquivos xyz:
		mol3D.readfromxyz(fname)
        mol3D.writexyz(fname)    
"""
