'''
Created on Dec 1, 2012

@author: WenhaoMIT
'''
import unittest
import numpy as np
from pymatgen.surface.surfacebasistransform import SurfaceOrientedBasisTransform
import unittest
import os
import numpy as np
from pymatgen.surface.makesurface import VaspSurface
from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
import numpy as np
from pymatgen.core.sites import PeriodicSite
from pymatgen.io.vaspio.vasp_input import Poscar
from pymatgen.io.cifio import CifParser
from pymatgen.core.structure_modifier import StructureEditor
from pymatgen.io.vaspio.vasp_output import Vasprun
from pymatgen.surface.surftool import Surftool

test_dir = os.path.join(os.path.dirname(__file__))

class Surftooltest(unittest.TestCase):
        
    def setUp(self):
        self.s=SurfaceOrientedBasisTransform()
        
    def testbuildbasis(self):
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR'))
        bulkstruct = p.struct
        maxindex=np.array([1,1,3])
        sconv=Surftool().getconventional(bulkstruct)
        pg=Surftool().getpg(sconv)
        listofsobats=self.s.transformbasis(sconv,maxindex, pg)
        for sobat in listofsobats:
            print Poscar(sobat.structure())
            
        
        
            
               
if __name__ == '__main__':
    unittest.main()