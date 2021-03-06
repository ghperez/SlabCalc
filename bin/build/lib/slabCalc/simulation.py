"""
Created on Sun Aug 30 13:33:26 2020

    SlabCalc simulation
    -------------------
    Objects for simulations using structures built with SlabCalc

author: Gabriel
institution:  Universidade Federal do ABC (UFABC)
"""
from slabCalc import *
import os
import pickle

class Sim(object):
	"""
	Class for running QE calculations with the slabs
	"""
	def __init__(self):
		"""
		"""
		self.slabs = list()
		self.dat   = list()
        
	def create_slabs(self,*params):
		"""
		(list) params : list of dictionaries containing slab's parameters (see set_molecule_parameters documentation)
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
		Build all slabs in the self.slabs list
		"""
		for i in range(len(self.slabs)):
			self.slabs[i].build(silent=silent)
			self.dat[i]["status"] = "not calculated"
		if not silent:
			print("All slabs built!")

	def set_qe(self,calc):
		"""
		Set QE calculations for all slabs in self.slabs list
		"""
		for i in range(len(self.slabs)):
			self.dat[i]["calc"] = self.slabs[i].set_qe_calculation(calc)
			
	def run_qe(self,cmd,silent=False,save_when_done=True,savefile="sim.dat"):
		"""
		Runs Quantum Espresso calculation for ALL slabs in self.slabs list
		"""
		for i in range(len(self.dat)):
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
			
		if "temp.pickle" in os.listdir("."):
			os.remove("temp.pickle")
		if save_when_done:
			self.save(savefile)
				 
	def calculate_slab(self,i,cmd,calc,input_string,save_steps=False,savefile="temp.pickle"):
		"""
		Runs Quantum Espresso calculation for ith slab in self.slabs

		(str) cmd          : terminal command for running Quantum Espresso
		(str) input_string : generated with the QE's calc object
		(bool) save_steps  : if true then will save the simulation after each step in the file
							 specified in "savefile" string  
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
		
		if save_steps:
			self.save(savefile)
		
		try:
			out  = calc.run(cmd, input_string, saveout, outfile, savecoords, coordsfile)
			if out.jobdone:
				self.dat[i]["output"] = out
				self.dat[i]["status"] = "calculated"
			else:
				print("!!! Job not done")
		except:
			print("!!! Error while calculating slab %i energy"%(i+1))
			
		if save_steps:
			self.save(savefile)
				
	def save(self,fname="sim.dat"):
		"""
		Saves the simulation state to a binary file (pickle)
		"""
		with open(fname,"wb") as f:
			pickle.dump(self.dat, f)
			
	def load(self,fname="sim.dat"):
		"""
		Loads simulation from a binary file (pickle)
		"""
		with open(fname,"rb") as f:
			self.dat = pickle.load(f)
			
	def _set_qe_params(self,i):
		params = ["saveinp","inpfile","saveout",
				  "outfile","savecoords","coordsfile"]
		for p in params:
			if p not in self.dat[i].keys():
				self.dat[i][p] = None
