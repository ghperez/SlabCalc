# SlabCalc3.0 Log

## TO-DOS
- fix rotational parameters configuration in *mol3D()*
- change *copymol3D()* so it can copy the new properties introduced 10/11/2020

## 11/10/2020
- **slabs are created using molecule and surface passed as mol3D object**
- removed labels
- set extents to [0,0,0] if extents nor cell_vector were passed
- changed method *write_qe_input()*'s name to *set_qe_calculation()*
- now *set_qe_calculation()* receives only some qe package calculation object and configure it
- copyslab() and clear() adapted to changes
- new mol3D object methods and properties (see documentation)
	- **molecule parameters** need to be set before building routine
	- implemented *set_molecule_parameters()*
	- implemented *get_coord_between()*
	
## 12/10/2020
- made minimal example building routine
- fixed rotational parameters configuration
- implemented *copy_molecule_parameters()* method and overrid method *copymol3D()* 

