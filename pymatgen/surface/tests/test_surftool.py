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
        v1=np.array([2,1,0])
        v2=np.array([1,0,1])
        basis=[[1, 0.00, 0.00],
               [0.00, 1, 0.00],
               [0.00, 0.00, 1]]
        self.s.getmillerfrom2v(basis, v1, v2)
        
    def testgetsymfam(self):
        pg=np.array([[1.0000,0,0,0,0],
                    [0,180.0000,0,-1.0000,0],
                    [0,180.0000,-0.8660,0.5000,0],
                    [0,180.0000,0.8660,0.5000,0],
                    [2.0000,180.0000,-0.8660,-0.5000,0],
                    [2.0000,180.0000,0.0866,0.0500,0],
                    [2.0000,180.0000,0,1.0000,0],
                    [3.0000,120.0000,0,0,1.0000],
                    [3.0000,240.0000,0,0,1.0000],
                    [-3.0000,120.0000,0,0,1.0000],
                    [-3.0000,240.0000,0,0,1.0000]])
        basis=np.array([[0, 4.8050422055, 0],
                        [4.1612866760, -2.4025211027, 0],
                        [0, 0, 13.1105423350]])
        C=self.s.getsymfam(basis, np.array([1, 1, 3]), pg)
        print C
        
if __name__ == '__main__':
    unittest.main()