'''
Created on Dec 20, 2012

@author: WenhaoMIT
'''


import json
import random
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.serializers.json_coders import PMGJSONDecoder
from pymatgen.core.operations import SymmOp
from pymatgen.core.structure_modifier import StructureEditor
from pymatgen.core.structure_modifier import SupercellMaker
from pymatgen.io.smartio import read_structure
from pymatgen.analysis.structure_fitter import StructureFitter, \
    shear_invariant, sqrt_matrix
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
from pymatgen.surface import surfacebasistransform as bb
from pymatgen.surface.surfacebasistransform import Sobat

test_dir = os.path.join(os.path.dirname(__file__))

class makesurfacetest(unittest.TestCase):
        
    def setUp(self):
        self.vs=VaspSurface()
        
    def testplayground(self):
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR'))
        bulkstruct = p.struct
        news = BVAnalyzer().get_valences(bulkstruct)
        print news
        
    def testsame(self):
        p = Poscar.from_file(os.path.join(test_dir, 'P1'))
        P1 = p.struct
        p = Poscar.from_file(os.path.join(test_dir, 'P2'))
        P2 = p.struct
        print StructureMatcher().fit(P1,P2)        
               
if __name__ == '__main__':
    unittest.main()