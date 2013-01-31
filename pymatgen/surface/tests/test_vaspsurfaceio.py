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
from pymatgen.surface.vaspsurfaceio import SurfaceIO, CreateSurfaceInputs

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

class vaspsurfaceiotest(unittest.TestCase):
        
    def setUp(self):
        
        self.vs=VaspSurface()
        self.sio=SurfaceIO()
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR'))
        bulkstruct = p.struct
        maxindex=np.array([2])
        surfdict=self.vs.StandardSurfDict(bulkstruct,maxindex)        
        self.surfdict=surfdict   
            
    def testCreateSurfaceInputs(self):
        name="Al2O3"
        homedirectory=test_dir
        slabdict=self.surfdict
        CreateSurfaceInputs(name, homedirectory, slabdict)
        
        
    '''
    
    def testsurfincar(self):
        slabs=self.vs.make_primitive_slab(self.listofsurfs)
        slab=slabs[0]
        surfincar=self.sio.surfincar(slab,True)
        print surfincar
        
    def testrelaxlayers(self):
        slabs=self.vs.make_primitive_slab(self.listofsurfs)
        slab=slabs[0]
        seldyn=self.sio.relaxslab(slab,0,False)
        print Poscar(slab,comment=None,selective_dynamics=seldyn)
    ''' 
        

if __name__ == '__main__':
    unittest.main()