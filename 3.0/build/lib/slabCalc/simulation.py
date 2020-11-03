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
		self.slabs = list()
		self.dat   = list()
        
	def create_slabs(self,*params):
		"""
		(list) params : lista de dicionarios de parametros
		"""
		if (type(params[0])==list or type(params[0])==tuple):
			params = params[0]
			
		for i in range(len(params)):
			iparams = params[i].copy()
			surface = iparams.pop("surface")
			
			cell_vector = None
			extents = None
			if "cell_vector" in iparams.keys():
				cell_vector = iparams["cell_vector"] 
			elif "extents" in iparams.keys():
				extents = iparams["extents"]
			
			molecules = list()	 
			for molecule in iparams.pop("molecules"):
				molecule.set_molecule_parameters(iparams)
				molecules.append(molecule)
				
			self.dat.append(iparams)
			self._set_qe_params(i)
				
			islab = Slab(surface,molecules,cell_vector=cell_vector,extents=extents)
			
			self.slabs.append(islab.copy())
			
	def build_slabs(self,silent=False):
		"""
		"""
		for i in range(len(self.slabs)):
			self.slabs[i].build(silent=silent)
			self.dat[i]["status"] = "not calculated"
		if not silent:
			print("All slabs built!")

	def set_qe(self,calc):
		"""
		"""
		for i in range(len(self.slabs)):
			self.dat[i]["calc"] = self.slabs[i].set_qe_calculation(calc)
			
	def run_qe(self,cmd,silent=False,save_when_done=True,savefile=None):
		"""
		"""
		for i in range(len(self.slabs)):
			saveinp = self.dat[i]["saveinp"]
			inpfile = self.dat[i]["inpfile"]
			
			if saveinp and inpfile==None:
				inpfile = "slab%i.in"%(i+1)
		
			calc = self.dat[i]["calc"]
			if self.dat[i]["status"]=="calculated":
				if not silent:
					print("Slab %i already calculated!"%(i+1))
				continue
			elif self.dat[i]["status"]=="paused":
				calc.system["restart_mode"] = "\"restart\""
				
			input_string = calc.build_input(saveinp,inpfile)
			
			self.calculate_slab(i,cmd,calc,input_string,save_steps=True,savefile="temp.pickle")
			
		os.remove("temp.pickle")
		if save_when_done:
			self.save(savefile)
				 
	def calculate_slab(self,i,cmd,calc,input_string,save_steps=False,savefile=None):
		"""
		"""
		saveout    = self.dat[i]["saveout"]
		outfile    = self.dat[i]["outfile"]
		savecoords = self.dat[i]["savecoords"]
		coordsfile = self.dat[i]["coordsfile"]
		
		if saveout and outfile==None:
			outfile = "slab%i.out"%(i+1)
		if savecoords and coordsfile==None:
			coordsfile = "slab%i_final_coords.xyz"%(i+1)
			
		self.dat[i]["status"]="paused"
		out  = calc.run(cmd, input_string, saveout, outfile, savecoords, coordsfile)
		
		if out.jobdone:
			try:
				self.dat[i]["output"] = out
				self.dat[i]["status"] = "calculated"
			except:
				print("!!! Error while calculating slab %i energy"%(i+1))
		if save_steps:
			self.save(savefile)
				
	def save(self,fname="sim.dat"):
		"""
		"""
		with open(fname,"wb") as f:
			pickle.dump(self.dat, f)
			
	def load(self,fname="sim.dat"):
		"""
		"""
		with open(fname,"rb") as f:
			self.dat = pickle.load(f)
			
	def _set_qe_params(self,i):
		"""
		"""
		params = ["saveinp","inpfile","saveout",
				  "outfile","savecoords","coordsfile"]
		for p in params:
			if p not in self.dat[i].keys():
				self.dat[i][p] = None
