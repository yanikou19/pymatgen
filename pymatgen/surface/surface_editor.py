'''
Created on Dec 5, 2012

@author: weichen
'''
from pymatgen.surface.surface import Surface
from pymatgen.core.structure_modifier import StructureEditor

class SurfaceEditor(StructureEditor):
    def __init__(self, surface):
        StructureEditor.__init__(self, surface)
        self._original_surface = Surface.fromstruct(self._original_structure, surface.indices, surface.vacuum)
        self._indices = surface.indices
        self._vac = surface.vacuum
    
    @property
    def modified_surface(self):
        coords = [site.frac_coords for site in self._sites]
        species = [site.species_and_occu for site in self._sites]
        indices = self._indices
        vac = self._vac
        props = {}
        if self._sites[0].properties:
            for k in self._sites[0].properties.keys():
                props[k] = [site.properties[k] for site in self._sites]
        return Surface(self._lattice, species, coords, indices, vac,
                         site_properties=props)
    
    @property
    def original_surface(self):
        return self._original_surface