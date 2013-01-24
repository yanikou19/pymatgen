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
        slablength=15
        vacuumlength=10
        utrfix=0
        shouldcheckpolar=1
        '''Creates a list of surfaces with above parameters'''
        listofsurfs=self.vs.makesurface(bulkstruct,maxindex,slablength,vacuumlength,utrfix,shouldcheckpolar)
        
        
        '''Create a dictionary of surfaces - key is index, value is a list of surfaces with that index'''
        print len(listofsurfs)
        indices=np.zeros((len(listofsurfs),3))
        for ii in range(0,len(listofsurfs)):
            indices[ii]=listofsurfs[ii].indices
        indices=np.array([np.array(x) for x in set(tuple(x) for x in indices)])
        I = np.argsort(indices[:,0])
        Millers=[]
        for row in range(0,len(indices)):
            Millers.append(indices[row])
        surfdict=defaultdict(list)
        for I in Millers:
            indexlistofsurfs=[]
            for surf in listofsurfs:
                if str(surf.indices)==str(I):
                    indexlistofsurfs.append(surf)
            surfdict[str(I)]=indexlistofsurfs
            
        '''Eliminate polar and duplicate surfaces in those lists'''
            
        for key in surfdict.keys():
            print key
            listofsurfs=surfdict[key]
            listnopolarsurfs=self.vs.checkpolar(listofsurfs)
            print len(listnopolarsurfs)
            listnoduplicates=self.vs.checksame(listnopolarsurfs)
            print len(listnoduplicates)
            surfdict[key]=listnoduplicates
        
        ''' Calculate surface energy from GULP for all these surfaces '''
        for key in surfdict.keys():
            for surf in surfdict[key]:
                E=self.vs.binoxi_gulp_surface_energy(surf)
                print key+"  "+str(E)

               
if __name__ == '__main__':
    unittest.main()