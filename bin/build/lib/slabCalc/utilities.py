#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

    SlabCalc Utilities
    ------------------
    Funções usadas em outros módulos do pacote SlabCalc
    
    Este módulo também pode ser usado para checar os cálculos em andamento
    usando o script no terminal:
        python slabCalc.utilities -check <nome do arquivo json contendo o simdata>
"""
import qe.pw as pw
import os
import json

def run(cmd,simulation,slab_index,option="fs",savein=False,infile=None,saveout=False,outfile=None,savecoords=False,coordfile=None):
	"""
	Roda os cálculos do qe.pw para a superfície
	"""
	slab = simulation.simdata[slab_index]
	calc = read_from_string(slab["input_string"])
	set_restart_mode(calc,option)
	calc.control["prefix"] = "\"%s\""%slab["label"]
	istring = calc.build_input(savein,infile)
	try:
		simulation.save("temp.json")
		out = calc.run(cmd,istring,saveout,outfile,savecoords,coordfile)
		print(out.jobdone)
		if out.jobdone:
			slab["energy"] = out.energy[-1]
			slab["status"] = "calculated"
			print("Energy of slab%i calculated"%slab["index"])
		else:
			slab["status"] = "paused"
	except:
		print("Error calculating energy of slab %i"%slab["index"])

def set_restart_mode(calc,option):
	"""
	Muda o parâmetro restart_mode do input do Quantum Espresso conforme
	o parâmetro option
	
	fs -> from_scratch
	r  -> restart

	"""
	if option=="fs":
		calc.control["restart_mode"]="\"from_scratch\""
	elif option=="r":
		calc.control["restart_mode"]="\"restart\""

def read_from_string(istring):
	"""
	Configura um objeto pw calc a partir de uma string
	"""
	calc = pw.calc()

	with open("ifile.in","w") as f:
		f.write(istring)

	calc.read_input("ifile.in")
	os.remove("ifile.in")

	return calc

def remove_temp():
	"""
	Remove arquivos temporários escritos durante os cálculos
	"""
	if "temp.json" in os.listdir("."):
		os.remove("temp.json")

def check_status(fname="temp.json"):
	"""
	Função para checar o estado da simulação
	"""
	with open(fname,"rb") as f:
		data = json.load(f)
	print("  SLAB    STATUS\n")
	for slab in data:
		print("%5d     %s"%(slab["index"],slab["status"]))

if __name__=="__main__":
	import sys
	n = len(sys.argv)
	if n>1:
		if sys.argv[1]=="-check":
			if n>2:
				check_status(sys.argv[2])
			else:
				check_status()

