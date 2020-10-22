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
import numpy as np
import os
import csv
import json

class Sim(object):
	"""
    Classe usada para rodar os cálculos com as superfícies
	"""
	def __init__(self):
		self.slabs      = list()
		self.energies   = list()
        
	def create_slabs(self,*params):
		"""
		(list) params : lista de dicionarios de parametros
		"""
		if (type(params[0])==list or type(params[0])==tuple):
			params = params[0]
			
		for i in range(len(params)):
			iparams = params[i].copy()
			surface = iparams.pop("surface")
			molecule = iparams.pop("molecule")
			cell_vector = None
			extents = None
			if "cell_vector" in iparams.keys():
				cell_vector = iparams.pop("cell_vector")
			elif "extents" in iparams.keys():
				extents = iparams.pop("extents")
					
			molecule.set_molecule_parameters(iparams)
			
			islab = Slab(surface,molecule,cell_vector=cell_vector,extents=extents)
		
			self.slabs.append(islab.copy())
			
	def build_slabs(self,silent=False):
		for slab in self.slabs:
			slab.build(silent=silent)
		if not silent:
			print("All slabs built!")
