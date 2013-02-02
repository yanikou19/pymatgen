'''
Created on Jan 31, 2013

@author: wenhao
'''
import numpy as np
import unittest
from pymatgen.surface.wulff import wulff
from collections import defaultdict
from pymatgen.io.vaspio.vasp_input import Poscar
import os

test_dir = os.path.join(os.path.dirname(__file__))

class WulffTest(unittest.TestCase):
    def setUp(self):
        self.indexSEdict=dict()
        self.SE="""0 0 1 1
0 1 -2 1.05
1 1 -3 1.06
1 1 0 1.12
1 0 0 1.124
3 1 0 1.8
1 0 -1  1.07"""
        surfdict=defaultdict()
        for line in self.SE.split("\n"):
            I=line.split()
            strI="["+str(int(I[0]))+","+str(int(I[1]))+","+str(int(I[2]))+"]"
            E=float(I[3])
            surfdict[strI]=E
        self.indexSEdict=surfdict
        p = Poscar.from_file(os.path.join(test_dir, 'Wulfftestposcar'))
        self.convstruct = p.struct
        
    def test_wulff(self):
        print wulff(self.convstruct, self.indexSEdict).getEAbvWulff()
        
if __name__ == '__main__':
        unittest.main()
        
        
        
