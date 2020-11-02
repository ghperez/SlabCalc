from slabCalc import *
import qe.pw as pw
import os

PSEUDO_DIR  = os.path.join("..","..","pseudopotentials","")
CALC_FROM_INPUT = False
INPUT_MODEL = "input_model"
CMD = "mpirun -np 8 pw.x"
SAVEOUT = True
OUTFILE = "calc.out"
SAVECOORDS = True
COORDSFILE = "final_coords.xyz"

PREFIX = "graph+hydrogen"

a = 2.46 # surface cell parameter in angstroms
n,m = 2,1 # surface repetition numbers

def set_calc(slab):
	calc = pw.calc()
	
	if CALC_FROM_INPUT:
		calc.read_input(INPUT_MODEL)
	else:
		""" SET THE CALCULATION PARAMETERS HERE """
		# &CONTROL
		calc.control["calculation"]   = "\"relax\""
		calc.control["restart_mode"]  = "\"from_scratch\""
		calc.control["pseudo_dir"]    = "\""+PSEUDO_DIR+"\""
		calc.control["prefix"]        = "\"simple_example\""
		calc.control["outdir"]        = "\"outputs\""
		calc.control["etot_conv_thr"] = "1.0E-5"
		calc.control["forc_conv_thr"] = "1.0D-4"
		calc.control["nstep"]         = 200
		
		# &SYSTEM
		calc.system["ibrav"]       = 0
		calc.system["celldm(1)"]   = a*1.889725989 # cell parameter in bohr
		#calc.system["nat"]         =
		#calc.system["ntyp"]        =
		calc.system["nspin"]       = 1
		calc.system["occupations"] = "\"smearing\""
		calc.system["degauss"]     = 0.02
		calc.system["ecutwfc"]     = 32
		calc.system["ecutrho"]     = 320.0
		calc.system["smearing"]    = "\"gaussian\""
		calc.system["input_dft"]   = "\"vdW-DF\""
		
		# &ELECTRONS
		calc.electrons["mixing_beta"]      = 0.5
		calc.electrons["electron_maxstep"] = 100
		calc.electrons["conv_thr"]         = "1.0D-6"
		calc.electrons["diagonalization"]  = "\"david\""
		calc.electrons["startingpot"]      = "\"atomic\""
		calc.electrons["startingwfc"]      = "\"atomic+random\""
		
		# &IONS
		calc.ions["ion_dynamics"] = "\"bfgs\""
		
		# CELL_PARAMETERS
		calc.cell_parameters_units = "alat"
		calc.v1 = [n,0,0]
		calc.v2 = [0,m*sqrt(3),0]
		calc.v3 = [0,0,10]
		
		# ATOMIC_SPECIES
		calc.pseudopotential = ["C.pbe-rrkjus.UPF","H.pbe-rrkjus_psl.1.0.0.UPF"]
		
		# ATOMIC_POSITIONS
		calc.atomic_positions_units = "angstrom"
		
		# K_POINTS
		calc.k_points_type = "automatic"
		calc.nk1 = 4
		calc.nk2 = 4
		calc.nk3 = 1
		calc.sk1 = 0
		calc.sk2 = 0
		calc.sk3 = 0
							
	slab.set_qe_calculation(calc)
	
	return calc
	
if __name__=="__main__":
	
	slab = Slab()
	slab.readfromxyz("%s.xyz"%PREFIX)
	calc = set_calc(slab)
	istring = calc.build_input(True,"%s.in"%PREFIX)
	out = calc.run(CMD, istring, SAVEOUT, OUTFILE, SAVECOORDS, COORDSFILE)
	
	if out.jobdone:
		print(45*"="+"\n")
		print("Final Energy = %f eV"%out.energy[-1])
		print("\n"+45*"=")
	else:
		print("Something went wrong")
