from pymatgen.core.structure import Structure
from pymatgen.io.aseio import AseAtomsAdaptor
import ase.lattice.general_surface as gensurface

class Surface(Structure):
    '''
    pymatgen Surface object
    '''

    def __init__(self, lattice, species, coords, indices, vac):
        '''
        Create a surface slab.
        Args:
            lattice, species, coords:
                as defined in pymatgen Structure class
            indices:
                surface indices
            vac:
                vacuum spacing of the slab
        '''
        Structure.__init__(self, lattice, species, coords)
        
        #TODO: add hcp support
        if len(indices) == 3:
            self._indices = indices
        elif len(indices) == 4:
            self._indices = indices
            self._ishcp == True
        else:
            raise ValueError("Invalid surface indices!") 
        
        self._vacuum = vac
                       
        
    @property
    def indices(self):
        return self._indices
        
    @property
    def vacuum(self):
        return self._vacuum
        
    def __str__(self):
        outs = "%s surface slab with %s Ang vacuum\n" % (self.indices, self.vacuum)
        return outs + Structure.__str__(self)
       
    @staticmethod
    def fromstruct(struct, indices, vac):     
        return Surface(struct.lattice, struct.species, struct.frac_coords, indices, vac) 
        
    @property
    def to_dict(self):
        '''
        dict representation of Surface
        '''
        d = super(Surface,self).to_dict
        d["indices"] = list(self._indices)
        d["vacuum"] = self._vacuum
        return d
     
    @staticmethod               
    def surfaceGenerator(bulkstruct, indices, layers, vac, tole=1e-10):
        asebulk = AseAtomsAdaptor.get_atoms(bulkstruct)
        asesurface = gensurface.surface(asebulk, indices, layers, vacuum=vac, tol=tole)
        struct = AseAtomsAdaptor.get_structure(asesurface)
        slab = Surface.fromstruct(struct, indices, vac)
        return slab 