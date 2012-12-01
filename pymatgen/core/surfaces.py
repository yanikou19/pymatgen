'''
Created on Nov 2, 2012

@author: WenhaoMIT
'''

from pymatgen.core.structure import Structure


class SlabStructure(Structure):
    
    def __init__(self, lattice, species, coords, millerindex, termination, numlayers, relax=1, vacuum=10, offstoich=0,adsorbates=0, surfunitcell=[1,1]):
        Structure.__init__(self,lattice, species, coords)
        self.millerindex = millerindex
        self.term = termination
        self.numlayers = numlayers
        self.relax = relax
        self.vacuum = vacuum
        self.offstoich = offstoich
        self.adsorbates = adsorbates
        self.surfunitcell= surfunitcell
        
    def __str__(self):
        output = Structure.__str__(self)
        output += "\nMiller Index : {}".format(self.millerindex)
        return output
        
'''
slab = SlabStructure([[10,0,0],[0,10,0],[0,0,10]], ["Fe","Fe"], [[0,0,0], [0.5,0.5,0.5]], [1,1,1],1,1)
print slab


from pymatgen.io.vaspio.vasp_input import Poscar
p = Poscar(slab)
print p
'''