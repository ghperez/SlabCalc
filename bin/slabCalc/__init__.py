#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    SlabCalc
    --------
    Conjunto de objetos para o estudo de superfícies 2D utilizando métodos
    computacionais.

autor:       Gabriel
instituição: Universidade Federal do ABC (UFABC)
"""

from molSimplify.Classes.mol3D import*
from molSimplify.Scripts.cellbuilder_tools import*
from molSimplify.Scripts.cellbuilder import*
import qe.pw as pw
import os

class Slab(mol3D,object):
    """
    Classe da superfície

    subclasse da classe mol3D, para mais informações leia a documentação do molSimplify
    """
    def __init__(self,sfile,mfile,cell_vector=None,extents=None,label="slab"):
        """
        Método Construtor

        Propriedades
              (str)  sfile         :  arquivo xyz da superfície
              
              (str)  mfile         :  arquivo xyz da molécula
              
             (list)  scell_vector  :  vetores de célula, usados para determinar
                                      a extensão da superfície. Você pode tanto
                                      passar este parâmetro ou o parâmetro extents
                                      
             (list)  extents       :  um vetor de coordenadas dos limites da superfície
             
           (string)  label         :  nome para rotular os arquivos
        """
        super(Slab, self).__init__()
        self.label = label
        # READING STRUCTURES FROM XYZ FILE
        self.readfromxyz(sfile)
        self.sfile = sfile
        self.molecule = mol3D()
        self.molecule.readfromxyz(mfile)
        self.mfile = mfile

        if cell_vector!=None:
            self.extents = find_extents_cv(cell_vector)
        elif extents!=None:
            self.extents = extents

    def build(self,ssite,msite,adist,rot=False,axis=None,angle=None,rpoint=None):
        """
        Insere a molécula sob a superfície

            (list) ssite : ponto de alocação na superfície
            (list) msite : ponto que vai ser alinhado na molécula
           (float) adist : distância de alinhamento em angstrons

           Rotational Parameters
           ---------------------
            (bool) rot    : if set to true then the molecule will be rotate
           (float) angle  : ângulo
            (list) axis   : eixo
            (list) rpoint : ponto por onde passa o eixo de rotação
        """
        # INITIALISING THE COMBINED SLAB
        cslab = mol3D()
        cslab.copymol3D(self)
        # INITIALISING THE PAYLOAD MOLECULE
        payload = mol3D()
        payload.copymol3D(self.molecule)
        # ROTATION
        if rot==True:
            if rpoint==None:
                rpoint = payload.centersym()
            try:
                rotate_around_axis(payload,rpoint,axis,angle)
                print("Molecule Rotated")
            except:
                print("Error during molecule rotation")
        # TRANSLATION
        trans_vec = vecdiff(ssite,msite)
        payload.translate(trans_vec)
        try:
            payload.translate([0, 0, self.extents[2]+adist])
        except:
            print("Missing extents information\nUsing z=0 as reference to align distance")
            payload.translate([0, 0, adist])
        # PLACEMENT METHOD
        self.combine(payload)

        self.sanity = self.sanitycheck(silence=True)
        if self.sanity[0]:
            print("\n"+25*"=")
            print("Molecule might be overlapping!\nThe minimum distance is %f.\nStructure: %s"%(self.sanity[1],self.label))
            print(25*"="+"\n")

        print("Slab built")

    def write_qe_input(self,inpmodel=None,saveinp=None,inpfile=None):
        """
        Retorna um input do Quantum Espresso em formato de string

        calc : objeto da biblioteca QE
        """
        calc = pw.calc()
        if inpmodel!=None:
            calc.read_input(inpmodel)

        if saveinp and inpfile==None:
            if "inputs_qe" not in os.listdir("."):
                os.mkdir("inputs_qe")
            inpfile = "./inputs_qe/%s.in"%self.label

        #calc.x,calc.y,calc.z,calc.atomic_species,calc.atomic_mass,calc.atom_type= ([],[],[],[],[],[])
        calc.system["nat"] = len(self.atoms)

        for atom in self.atoms:
            calc.x.append(atom.coords()[0])
            calc.y.append(atom.coords()[1])
            calc.z.append(atom.coords()[2])
            calc.atom_type.append(atom.symbol())
            if atom.symbol() not in calc.atomic_species:
                calc.atomic_mass.append(atom.mass)
                calc.atomic_species.append(atom.symbol())

        calc.system["ntyp"] = len(calc.atomic_species)

        return calc.build_input(saveinp,inpfile)

    def copyslab(self,other):
        """
        Método para fazer cópias das superfícies
        """
        self.sfile = other.sfile
        self.mfile = other.mfile
        self.extents = other.extents
        self.label = other.label
        self.molecule = other.molecule
        self.copymol3D(other)

    def clear(self,new_label=None):
        """
        Reinicia a superfície
        """
        if new_label==None:
            new_label=self.label
        self.__init__(self.sfile,self.mfile,extents=self.extents,label=new_label)

        print("Slab reset")

"""
    Funções úteis para obter os pontos de alocação:

        - center_of_sym([molecule.getAtom(i).coords() for i in inds])
            # retorna as coordenadas do ponto entre os átomos de índicies
            especificados em "inds"
        - mol3D_object.centersym()
            # retorna o centro de simetria de um objeto mol3D
"""