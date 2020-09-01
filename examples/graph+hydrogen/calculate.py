#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 30 12:26:31 2020

    GraphHydrogen
    -----------
    This module does the calculations for the GraphHydrogen structure using some
    simulation json file.

@author: gabriel
"""
from slabCalc import *

simfile = "simulation.json"
cmd = "mpirun -np 8 pw.x" # Quantum Espresso run command

tsim = Sim()
tsim.load(simfile)

tsim.calculate(cmd,savein=True,saveout=True,savecoords=True)

for slab in tsim.simdata:
	ssite = slab["ssite"]
	if ssite==[1.23, 1.4202816622064793, 0.0]:
		slab["ssite label"] = "hexagon center"
	elif ssite==[2.46, 0.7101408311032397, 0.0]:
		slab["ssite label"] = "above atom"
	elif ssite==[2.46, 1.4202816622064793, 0.0]:
		slab["ssite label"] = "above bond"

keys = tsim.simdata[0].keys()
keys.remove("input_string")
keys.remove("ssite")
keys.remove("axis")
keys.remove("angle")
tsim.write_results(keys,fname="sim_results.csv")

