from pymatgen.core.structure import SiteCollection, Structure
from pymatgen.core.lattice import Lattice
from pymatgen.io.aseio import AseAtomsAdaptor
import ase.lattice.general_surface as gensurface
from pymatgen.core.sites import PeriodicSite
import numpy as np
from pymatgen.serializers.json_coders import MSONable
from pymatgen.core.structure_modifier import StructureEditor

class Surface(Structure):
    '''
    pymatgen Surface object
    '''

    def __init__(self, lattice, species, coords, indices, vac, onesided=False, site_properties=None):
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
        self._onesided = onesided
                              
    @property
    def indices(self):
        return self._indices
        
    @property
    def vacuum(self):
        return self._vacuum
                   
    def __str__(self):
        outs = "{0} surface slab with {1} Ang vacuum\n".format(self.indices, self.vacuum)
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
        '''
        Generate a surface object with positive z lattice
        '''
        asebulk = AseAtomsAdaptor.get_atoms(bulkstruct)
        asesurface = gensurface.surface(asebulk, indices, layers, vacuum=vac, tol=tole)
        struct = AseAtomsAdaptor.get_structure(asesurface)
        zcard=struct.lattice.matrix[2][2]
        if zcard < 0:
            newmatrix=struct.lattice.matrix
            newmatrix[2][2]*=-1
            slab_normal=StructureEditor(struct)        
            slab_normal.modify_lattice(Lattice(newmatrix))
            return Surface.fromstruct(slab_normal.modified_structure, indices,vac) 
        else:
            return Surface.fromstruct(struct,indices,vac)
    
    def getSurfaceSites(self, onesided=False, subsurface=False):
        '''
        Get surface and subsurface sites
        '''
        if self._onesided != onesided:
            print "Slab has different number of surface layer!"      
        zcoords=self.zCoords
        sites = lambda layer: [s for s in self.sites if s.coords[2] == zcoords[layer]]
        if not subsurface:
            if onesided:
                return SurfaceSiteCollection([SurfaceSite(i) for i in sites(-1)])
            else:
                return SurfaceSiteCollection([SurfaceSite(i) for i in sites(-1)+sites(0)])
        if subsurface:
            if onesided:
                return SurfaceSiteCollection([SurfaceSite(i) for i in sites(-1)]+\
                       [SurfaceSite(i,atsurface=False) for i in sites(-2)])
            else:
                return SurfaceSiteCollection([SurfaceSite(i) for i in sites(-1)+sites(0)]+\
                       [SurfaceSite(i,atsurface=False) for i in sites(-2)+sites(1)])
      
    @property                 
    def zCoords(self):
        return np.unique(np.array([s._coords[2] for s in self.sites]), return_index=False)
                        
    @property
    def surfacecomposition(self):
        surfacesites=self.getSurfaceSites(onesided=self._onesided)
        return surfacesites.composition
                   
class SurfaceSite(PeriodicSite):
    '''
    Site class for atoms at surfaces and subsurfaces
    '''
    def __init__(self, site, atsurface=True):
        for i in site.__dict__:
            setattr(self,i,site.__dict__[i])
        self._atsurface = atsurface
    
    def __repr__(self):
        if self._atsurface:
            outs = ["Surface Site"]
        else:
            outs = ["Subsurface Site"]
        outs.append("abc : (%0.4f, %0.4f, %0.4f)" % tuple(self._fcoords))
        for k, v in self._species.items():
            outs.append("element    : %s" % k.symbol)
            outs.append("occupation : %0.2f" % v)
        return "\n".join(outs)

class SurfaceSiteCollection(SiteCollection):
    def __init__(self, site):
        self._sites=[]
        for i in site:
            if isinstance(i, SurfaceSite):
                self._sites.append(i)
            else:
                raise ValueError("Trying to add a site that is not a SurfaceSite!")  
    @property
    def sites(self):
        return self._sites
    
    def get_distance(self, i, j):pass
    