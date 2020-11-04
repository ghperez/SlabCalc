# SlabCalc3.0 Log

## TO-DOS
- not creating temp files (?)
- implement script for checking the calculations
- fix the documentation
- follow the state of Titanio's and Lab's calculations

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
- made minimal example building routine (simple example)
- fixed rotational parameters configuration
- implemented *copy_molecule_parameters()* method and overrid method *copymol3D()*
- updated *Slab()* object so one can create an instance without passing any parameter (new way to load slabs using xyz files)
- written the minimal example build and calculating routine

## 13/10/2020
- reviewed *simulation* module and raised ideas for improving it

## 15/10/2020
- implemented *create_slabs()* method that replaced *build()* of *Sim()* object
- implemented *build_slabs()* to build the slabs created using *create_slabs()* method 
- started writing the build script for *batch_simulation* example

## 19/10/2020
- updated Slab()'s method *copy*, now it returns another Slab() object
- fixed problem during copying slab's molecules
- finished build routine of *batch_simulation* example

## 22/10/2020
- Fixed and tested slab's *clear* method
- tryed to run the simple simple but it raised an error

## 24/10/2020
- fixed simple example input error
- started running simple example in the lab's computer
- changed the return of Slab's *set_qe_calculation* method

## 25/10/2020
- implemented *Sim's* methods *run_qe_calculations* and *calculate_slabs* for running QE calculations and getting the energies
- implemented *Slab's* method set_qe_parameters
- changed the way *Sim's* *create_slabs* method read dictionaries (remove .pop structure)

## 01/11/2020
- implemented checkpoints funcionality using pickle module (methods *save()* and *load()*)
- checkpoints with pickle doesn't work :( (TypeError: can't pickle SwigPyObject objects)
- reported error while creating slabs using Sim
- fixed error when copying slab's molecules

## 02/11/2020
- used *deepcopy()* function from library **copy** to clone *pw.calc()* objects
- used serialization of dictionaries with pickle to save the simulations
- implemented *save()* and *load()* methods for *Sim*
- split Sim's *run_qe_calculations* into *set_qe* and *run_qe* methods

- fixed bug when reading molecules in *Slab*'s \_\_init\_\_ method
- prepared&run simple and batch examples for graphene+hydrogen in lab's computer (results at examples folder)
- prepared&run batch example for graphene+benzene in Titanio cluster (waiting for results) 