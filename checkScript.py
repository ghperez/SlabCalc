from slabCalc.simulation import *

sim = Sim()
sim.load(input(">>> Enter loadfile: "))

print(2*" "+"%25s | %10s | %15s | %15s"%("Label", "Angle", "Status", "Final Energy"))
print(12*" "+65*"-")
for slab in sim.dat:
	label = slab["label"]
	angle = slab["angle"]
	sts = slab["status"]
	str_out = 2*" "+"%25s | %10f | %15s"%(label, angle, sts)
	if "output" in slab.keys():
		fenergy = slab["output"].energy[-1]
		str_out+=" | %15f"%(fenergy)
	print(str_out)
