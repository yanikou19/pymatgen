'''
Created on Nov 20, 2012

@author: weichen
'''

import unittest

from pymatgen.core.periodic_table import Element, Specie
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice

from pymatgen.surface.surface import Surface

class SurfaceTests(unittest.TestCase):
    def setUp(self):
        self.pt = Element("Pt")
        coords = list()
        coords.append([0, 0, 0])
        coords.append([0.5, 0.5, 0])
        self.lattice = Lattice([[3.92, 0.00, 0.00],
                                [0.00, 3.92, 0.00],
                                [0.00, 0.00, 3.92]])
        self.struct = Structure(self.lattice, [self.pt, self.pt], coords)

    def testSurfceGenerator(self):
        slab = Surface.surfaceGenerator(self.struct, (1,1,1), 5, 10)
        (lengths, angels) = slab.lattice.lengths_and_angles
        self.assertAlmostEqual(lengths[0], 5.5437, places=4)
        self.assertAlmostEqual(lengths[2], 29.0529, places=4)
        self.assertAlmostEqual(angels[0], 90)
        self.assertAlmostEqual(angels[2], 120)   

if __name__ == '__main__':
    unittest.main()