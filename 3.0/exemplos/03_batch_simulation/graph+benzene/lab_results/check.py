from slabCalc.simulation import *

sim = Sim()
sim.load("sim.dat")

for i in sim.dat:
	print(i["status"])
