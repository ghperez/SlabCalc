#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 30 13:33:26 2020

    SlabCalc simulation
    -------------------
    Objects for simulations using structures built with SlabCalc

@author: gabriel
"""
from slabCalc.utilities import *
from slabCalc import *
import os
import pickle

class Sim(object):
	"""
    Classe usada para rodar os cálculos com as superfícies
	"""
	def __init__(self):
		"""
		"""
		self.slabs      = list()
        
	def create_slabs(self,*params):
		"""
		(list) params : lista de dicionarios de parametros
		"""
		if (type(params[0])==list or type(params[0])==tuple):
			params = params[0]
			
		for i in range(len(params)):
			iparams = params[i].copy()
			surface = iparams["surface"]#.pop("surface")
			molecule = iparams["molecule"]
			cell_vector = None
			extents = None
			if "cell_vector" in iparams.keys():
				cell_vector = iparams["cell_vector"] 
			elif "extents" in iparams.keys():
				extents = iparams["extents"] 
					
			molecule.set_molecule_parameters(iparams)
			
			islab = Slab(surface,molecule,cell_vector=cell_vector,extents=extents)
			islab.set_qe_parameters(iparams)
		
			self.slabs.append(islab.copy())
			
	def build_slabs(self,silent=False):
		"""
		"""
		for slab in self.slabs:
			slab.build(silent=silent)
		if not silent:
			print("All slabs built!")

	def run_qe_calculations(self,calc,calculate=False,cmd=None,silent=False):
		"""
		"""
		for i in range(len(self.slabs)):
		
			saveinp = self.slabs[i].saveinp
			inpfile = self.slabs[i].inpfile
			
			if saveinp and inpfile==None:
				inpfile = "slab%i.in"%(i+1)
			
			self.slabs[i].set_qe_calculation(calc)
			input_string = calc.build_input(saveinp,inpfile)
			
			if calculate:
				if self.slabs[i].status=="calculated":
					if not silent:
						print("Slab %i already calculated!"%(i+1))
				else:
					self.save()
					self.calculate_slab(i,cmd,calc,input_string)
				 
	def calculate_slab(self,i,cmd,calc,input_string):
		"""
		"""
		saveout    = self.slabs[i].saveout
		outfile    = self.slabs[i].outfile
		savecoords = self.slabs[i].savecoords
		coordsfile = self.slabs[i].coordsfile
		
		if saveout and outfile==None:
			outfile = "slab%i.out"%(i+1)
		if savecoords and coordsfile==None:
			coordsfile = "slab%i_final_coords.xyz"%(i+1)
			
		out  = calc.run(CMD, input_string, saveout, outfile, savecoords, coordsfile)
		
		if out.jobdone:
			try:
				self.slabs[i].out = out
				self.slabs[i].status = "calculated"
			except:
				print("!!! Error while calculating slab %i energy"%(i+1))
				
	def save(self,fname="sim.pickle"):
		with open(fname,"wb") as pickle_out:
			pickle.dump(self.slabs,pickle_out)
		
	def load(self,fname="sim.pickle"):
		with open(fname,"rb") as pickle_in:
			pickle.load(pickle_in)
