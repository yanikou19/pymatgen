'''
Created on Nov 20, 2012

@author: weichen
'''

import unittest

from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.structure import Structure, Molecule
from pymatgen.core.lattice import Lattice

from pymatgen.surface.surface import Surface
from pymatgen.surface.adsorption import Adsorption,Adsorbate
import numpy as np
from pymatgen.transformations.standard_transformations import RotationTransformation


class SurfaceTests(unittest.TestCase):
    def setUp(self):
        self.pt = Element("Pt")
        self.au = Element("Au")
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.5, 0.5, 0])
        self.lattice = Lattice([[3.92, 0.00, 0.00],
                                [0.00, 3.92, 0.00],
                                [0.00, 0.00, 3.92]])
        self.struct = Structure(self.lattice, [self.au, self.pt], coords)
        
        self.h = Element("H")
        self.o = Element("O")
        self.acoords = list()
        self.acoords.append([0, 0, 0])
        self.acoords.append([0.76, 0.59, 0])
        self.acoords.append([-0.76,0.59, 0])
        
        #print self.struct

    def testSurfceGenerator(self):
        slab = Surface.surfaceGenerator(self.struct, (1,1,1), 5, 10)
        ads=Adsorbate([self.o, self.h, self.h], self.acoords)
        print slab
        print ads
        
        (lengths, angels) = slab.lattice.lengths_and_angles
        self.assertAlmostEqual(lengths[0], 5.5437, places=4)
        self.assertAlmostEqual(lengths[2], 29.0529, places=4)
        self.assertAlmostEqual(angels[0], 90)
        self.assertAlmostEqual(angels[2], 120)
        
        s=Adsorption.add_adsorbate_on_type(slab, ads, ['Au'], 2, (0,0), 0)
        print s       
if __name__ == '__main__':
    unittest.main()