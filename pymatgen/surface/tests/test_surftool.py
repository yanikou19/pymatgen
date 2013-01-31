'''
Created on Nov 30, 2012

@author: WenhaoMIT
'''

import unittest
import numpy as np
import numpy.linalg as npl
from numpy import pi
import math
from pymatgen.surface.surftool import Surftool
from pymatgen.io.vaspio.vasp_input import Poscar
from pymatgen.io.cifio import CifParser
import os


test_dir = os.path.join(os.path.dirname(__file__))


class Surftooltest(unittest.TestCase):
    
    def setUp(self):
        self.s=Surftool()
        
    def testangle(self):
        v1=np.array([1, 0, 0])
        v2=np.array([0, 1, 0])
        angle=self.s.ang(v1, v2)
        self.assertAlmostEqual(angle, pi/2, places=4)
    
    def testget2vectsinplane(self):
        basis=[[1, 0.00, 0.00],
               [0.00, 1, 0.00],
               [0.00, 0.00, 1]]
        maxindex=np.array([1, 2, 3])
        twovects,P=self.s.get2vectsinplane(basis, maxindex)
    
    def testtreflection(self):
        v=np.array([1,1,1])
        self.s.treflection(v)
        
    def testtrotation(self):
        v=np.array([1,1,1])
        self.s.trotation(120,v)
        
    def testgetmillerfrom2v(self):
        v2=np.array([ 0.,0.,13.11054233])
        v1=np.array([ 0.,4.80504221,  0.])
        basis=np.array([[0, 4.8050422055, 0],
                        [4.1612866760, -2.4025211027, 0],
                        [0, 0, 13.1105423350]])

    def testgetconventional(self):
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR'))
        bulkstruct = p.struct
        sconv=self.s.getconventional(bulkstruct)
        
    def testgetpg(self):
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR'))
        bulkstruct = p.struct
        sconv=self.s.getconventional(bulkstruct)
        pg=self.s.getpg(sconv)
        
    def testgetsymfam(self):
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR'))
        bulkstruct = p.struct
        sconv=self.s.getconventional(bulkstruct)
        pg=self.s.getpg(sconv)
        basis=sconv.lattice.matrix
        C=self.s.getsymfam(basis, np.array([1, 1, 3]), pg)
        
        
        

            
if __name__ == '__main__':
    unittest.main()