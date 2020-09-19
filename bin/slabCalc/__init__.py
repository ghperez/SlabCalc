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
    def __init__(self,surface,molecule,cell_vector=None,extents=None,label="slab"):
        """
        Método Construtor

        Propriedades
           (object)  surface       :  objeto mol3D da superfície
              
           (object)  molecule      :  objeto mol3D da molécula
              
             (list)  scell_vector  :  vetores de célula, usados para determinar
                                      a extensão da superfície. Você pode tanto
                                      passar este parâmetro ou o parâmetro extents
                                      
             (list)  extents       :  um vetor de coordenadas dos limites da superfície
             
           (string)  label         :  nome para rotular os arquivos
        """
        super(Slab, self).__init__()
        self.label = label
        # READING STRUCTURES FROM XYZ FILE
        self.molecule = molecule
        self.surface = surface

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
            (bool) rot    : se Verdadeiro então a molécula vai ser rotacionada
           (float) angle  : ângulo
            (list) axis   : eixo
            (list) rpoint : ponto por onde passa o eixo de rotação
        """
        # INITIALISING THE COMBINED SLAB
        cslab = mol3D()
        cslab.copymol3D(self.surface)
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
        # COMBINING MOLECULE AND SURFACE INTO SLAB
        self.copymol3D(cslab)
        self.combine(payload)

        self.sanity = self.sanitycheck(silence=True)
        if self.sanity[0]:
            print("\n"+25*"=")
            print("Molecule might be overlapping!\nThe minimum distance is %f.\nStructure: %s"%(self.sanity[1],self.label))
            print(25*"="+"\n")

        print("Slab built")

    def write_qe_input(self,calc,saveinp=None,inpfile=None):
        """
        Retorna um input do Quantum Espresso em formato de string confome um
        um objeto calc pré-configurado
        
        (object) calc      : objeto da biblioteca QE
        """
        if saveinp and inpfile==None:
            inpfile = "%s.in"%self.label

        for atom in self.atoms:
            calc.x.append(atom.coords()[0])
            calc.y.append(atom.coords()[1])
            calc.z.append(atom.coords()[2])
            calc.atom_type.append(atom.symbol())
            if atom.symbol() not in calc.atomic_species:
                calc.atomic_mass.append(atom.mass)
                calc.atomic_species.append(atom.symbol())

        calc.system["nat"] = len(self.atoms)
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
            
        Leitura e escrita de arquivos xyz:
            mol3D.writexyz()
            
"""
