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
import os
import csv
import json

class Sim(object):
    """
    Classe usada para rodar os cálculos com as superfícies
    """
    simdata = list()
    ssites = list()
    msites = list()
    adists = list()
    angles = list()
    axis = list()

    def build(self,slab,inpmodel=None,params=None,savein=False,save_slabs=False):
        """
        Constrói a superfície e armazena suas informações num dicionário
		"""

        self.set_params(params)

        if save_slabs and "structures" not in os.listdir("."):
            os.mkdir("./structures")
            
        if savein and "inputs" not in os.listdir("."):
            os.mkdir("inputs")

        ind = 0 # for identificating slabs

        print("\n*** Building Slabs ***")
        for s in self.ssites:
            for m in self.msites:
                for d in self.adists:
                    for a in self.angles:
                        for ax in self.axis:
                            print("\nssite = "),
                            print(s)
                            print("msite = "),
                            print(m)
                            print("adist = %f"%d)
                            print("angle = %f"%a)
                            print("axis = "),
                            print(ax)

                            slab_data = dict()
                            slab_data["ssite"]=s
                            slab_data["msite"]=m
                            slab_data["adist"]=d
                            slab_data["angle"]=a
                            slab_data["axis"]=ax
                            slab_data["index"]=ind
                            slab_data["label"]="slab%i"%ind
                            slab_data["energy"]=None
                            slab_data["status"]="not calculated"
                            slab.label=slab_data["label"]
							
                            print("Working on structure %s ..."%slab.label)
                            slab.build(s,m,d,rot=True,axis=ax,angle=a)
                            if save_slabs:
                                slab.writexyz("./structures/%s.xyz"%slab.label)
                            istring = slab.write_qe_input(inpmodel=inpmodel,
                                                          saveinp=savein,
                                                          inpfile="./inputs/%s.in"%slab.label) 
                                                          #,saveinp=True,inpfile="%s.in"%slab.label)
                            slab_data["input_string"] = istring
                            self.simdata.append(slab_data)
                            slab.clear()
                            print("%s built with success\n"%slab.label)
                            ind+=1
        print("*** Finish building routine ***\n")

    def calculate(self,cmd,savein=False,saveout=False,savecoords=False):
        """
        Calcula a energia das superfícies no dicionário simdata
        """
        print("\n*** Starting Energy Calculation Routine ***")

        for i in range(len(self.simdata)):
            slab = self.simdata[i]

            infile    = "slab%i.in"%(slab["index"]) if savein else None
            outfile   = "slab%i.out"%(slab["index"]) if saveout else None
            coordfile = "slab%i_coords.xyz"%(slab["index"]) if savecoords else None

            print("\nProcessing slab%i"%slab["index"])

            try:
                if slab["status"]=="calculated":
                    print("Energy already calculated")
                    continue
                else:
                    if slab["status"]=="not calculated":
                        print("Calculating from scratch")
                        option="fs"
                    elif slab["status"]=="paused":
                        print("Continuing from previous calculation")
                        option="r"
                    else:
                        print("Invalid status for slab%i."%slab["index"])
                        continue
                    run(cmd,self,i,option,savein,infile,saveout,outfile,savecoords,coordfile)
                    print("\n")
                    
            except KeyError:
                print("Please set the slab%i status."%slab["index"])
            self.save()
            
        remove_temp()
        print("*** Finished Energy Calculation Routine ***\n")

    def write_results(self,keys=None,fname=None):
        """
        Escreve o estado da simulação em um arquivo csv (fname)
        
        É possível fazer o controle de quais propriedades da superfície você
        deseja que sejam escritas através do parâmetro keys.
        """
        if fname==None:
            fname="results.csv"
        if keys==None:
            keys = self.simdata[0].keys()
            try:
                keys.remove("input_string")
            except:
                pass
            
        res = list()
        for i in range(len(self.simdata)):
            idict = dict()
            for key in keys:
                idict[key]=self.simdata[i][key]
            res.append(idict)
            
        with open(fname,"w") as f:
            writer = csv.DictWriter(f,fieldnames=keys)
            writer.writeheader()
            for slab in res:
                writer.writerow(slab)

    def save(self,fname=None):
        """
		Salva o dicionário simdata em um arquivo json
        """
        if fname==None:
            fname = "simdata.json"
        with open(fname,"w") as f:
            json.dump(self.simdata,f)

    def load(self,fname=None):
        """
        Carrega um dicionário simdata de um arquivo json
        """
        if fname==None:
            fname="simdata.json"
        with open(fname,'rb') as f:
            self.simdata = json.load(f)
                
    def set_params(self,params):
        """
        Configura os parâmetros da simulação através de um dicionário contendo
        os parâmetros.
        
        (dict) params : dicionário de parâmetros
        """
        if not params==None:
            self.ssites = params["ssites"]
            self.msites = params["msites"]
            self.adists = params["adists"]
            
            if "angles" in params.keys():
                self.angles = params["angles"]
                self.axis   = params["axis"]
            else:
                self.angles = [0]
                self.axis   = [[0,0,1]]

