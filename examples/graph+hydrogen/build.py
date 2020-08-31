#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 14:55:31 2020

    GraphHydrogen
    -----------
    This module builds the GraphHydrogen structure

@author: gabriel
"""
from slabCalc import *
from slabCalc.simulation import *
import numpy as np
import os
import json

tsim = Sim()
# Structures files
sfile = "graphene.xyz"
mfile = "hydrogen.xyz"

# Building parameters
cell_vector = [[2.46,0, 0],[0,4.2608449866194382,0],[0,0,0]]

slab = Slab(sfile,mfile,cell_vector)

ssite1 = center_of_sym([slab.getAtom(i).coords() for i in [2,5]]) # above hexagon center
ssite2 = slab.getAtom(5).coords() # above atom
ssite3 = center_of_sym([slab.getAtom(i).coords() for i in [5,6]]) # above bond
ssites =  [ ssite1, ssite2, ssite3] 
print(ssites)

msite = [slab.molecule.centersym()]

adists = list(np.arange(2,5,.5))

# Params dict
params = {"ssites":ssites,"msites":msite,"adists":adists}

# Building Routine
tsim.build(slab,inpmodel="./inputs/test_espresso.in",
           params=params,savein=True,save_slabs=True)

if "parameters.json" not in os.listdir("."):
    with open("parameters.json","w") as f:
    	json.dump(params,f)
    
tsim.save(fname="simulation.json")