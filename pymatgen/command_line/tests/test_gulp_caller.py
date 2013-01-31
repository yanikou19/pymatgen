'''
Created on Jan 22, 2013

@author: wenhao
'''
import unittest

from nose.exc import SkipTest

from pymatgen.command_line.gulp_caller import get_binoxi_gulp_energy,Tersoff_pot
from pymatgen.core.structure import Lattice, Structure
from pymatgen.core.periodic_table import Element
from pymatgen.util.io_utils import which
from pymatgen.io.vaspio.vasp_input import Poscar
import os

test_dir = os.path.join(os.path.dirname(__file__))
gulp_present = which('gulp')

class GulpCallerTest(unittest.TestCase):

    def setUp(self):
        if not gulp_present:
            print "No gulp"
            raise SkipTest("gulp not present. Skipping...")
        p = Poscar.from_file(os.path.join(test_dir, 'POSCAR'))
        self.structure = p.struct
            
    def test_make_tersoff_gulpinput(self):
        print get_binoxi_gulp_energy(self.structure)


if __name__ == '__main__':
    unittest.main()
