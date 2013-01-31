'''
Created on Dec 20, 2012

@author: WenhaoMIT
'''

from collections import defaultdict
import unittest
import os
import numpy as np
import numpy.linalg as npl
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
from pymatgen.surface import surfacebasistransform as bb
from pymatgen.surface.surfacebasistransform import Sobat

test_dir = os.path.join(os.path.dirname(__file__))

class makesurfacetest(unittest.TestCase):
        
    def setUp(self):
        self.vs=VaspSurface()
        
    def testmakesurface(self):
        '''Initialize from a POSCAR'''
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR'))
        bulkstruct = p.struct
        maxindex=np.array([2])
        surfdict=self.vs.StandardSurfDict(bulkstruct,maxindex)
        
        
        ''' Calculate surface energy from GULP for all these surfaces '''
        for key in surfdict.keys():
            for surf in surfdict[key]:
                E=self.vs.binoxi_gulp_surface_energy(surf)
                print key+"  "+str(E)

               
if __name__ == '__main__':
    unittest.main()