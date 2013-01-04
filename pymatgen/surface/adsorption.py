'''
Created on Dec 3, 2012

@author: weichen
'''
from pymatgen.surface.surface import Surface
from pymatgen.core.sites import Site
from pymatgen.surface.surface_editor import SurfaceEditor
from pymatgen.core.structure import Molecule
from pymatgen.core.periodic_table import Element, Specie
import numpy as np
from pymatgen.core.operations import SymmOp
import numpy.linalg as linalg
from scipy.optimize import fmin_slsqp


class Adsorption(object):
    '''
    classdocs
    '''
    def __init__(self, slab, adsorbates, position):
        '''
        Constructor
        '''
        self._slab = slab
        self._adsorbates = adsorbates
        self._position = position
        
    @staticmethod
    def add_adsorbate_on_type(slab, adsorbate, surfaceatom=[], height=2.0, coords=(0,0), siteindex=0):
        '''
        add an adsorbate at surface atoms
        
        Args:
            surfaceatom: List of surface atoms to attach adsorbates
            coords: displacement for adsorbates in the x-y plane
            siteindex: index for the site in a molecular adsorbate that is closest to the surface
        '''
        surfacesites=slab.getSurfaceSites(subsurface=False)
        specielist=[site.species_string for site in surfacesites]
        print specielist
        for i in surfaceatom:
            if i not in specielist:
                raise ValueError("No {0} atom found in the surface".format(surfaceatom))
        zcoords=slab.zCoords
        s=SurfaceEditor(slab)
        for site in surfacesites:
            if site.species_string in surfaceatom:
                if site.coords[2] == zcoords[0]:
                    ad_coords=[site.coords[0]+coords[0],site.coords[1]+coords[1],site.coords[2]-height]
                    up=True 
                elif site.coords[2] == zcoords[-1]:
                    ad_coords=[site.coords[0]+coords[0],site.coords[1]+coords[1],site.coords[2]+height]
                    up=False
                if isinstance(adsorbate, str):
                    s.append_site(adsorbate, ad_coords, True)
                if isinstance(adsorbate, Adsorbate):
                    ads = Adsorption.move_molecular_adsorbates(adsorbate, siteindex, ad_coords, up)
                    Adsorption.insert_molecular_adsorbates(ads, s)       
        return s.modified_surface
    
    @staticmethod
    def move_molecular_adsorbates(adsorbates,siteindex,position,up=True):
        if not isinstance(adsorbates, Adsorbate):
            raise ValueError("Adsorabates input is not a valid type")
        vec = adsorbates.molecule_adosrbate_optimizer(siteindex)[0]
        rotation = Adsorbate.rotation_matrix_to_z(vec)
        print rotation.rotation_matrix
        if up==True:
            pos = np.array(position)
        else:
            pos = np.array(position) * (-1)
        indexcoords = adsorbates[siteindex].coords   
        for i in adsorbates.sites:
            print i._coords-adsorbates[siteindex].coords
            i._coords=np.dot(rotation.rotation_matrix, (i.coords-indexcoords))+pos
            print i._coords
        return adsorbates
    
    @staticmethod
    def insert_molecular_adsorbates(sites, editor):
        for i in sites:
            editor.append_site(i.specie, i.coords, coords_are_cartesian=True)
            
    #TODO
    @staticmethod
    def hcpsite(): pass
    @staticmethod
    def fccsite(): pass
    @staticmethod
    def bridgesite(): pass

class Adsorbate(Molecule):
    
    def __init__(self,species,coords):
        Molecule.__init__(self,species,coords)

    @staticmethod
    def rotation_angle(x,y=0,len=0):
        '''
        Get the counter-clockwise rotation angle from (1,0) to (x,y) in a 2D plane
        Optional len if only positive rotation is possible
        ''' 
        if y==0 and len==0:
            raise ValueError("Either give a value for y or the length") 
        deg=np.arccos(x/((x**2+y**2)**0.5))
        if len != 0:
            return np.arccos(x/len)
        if y>=0:
            return deg
        else:
            return deg*(-1)
        
    @staticmethod
    def rotation_matrix_to_z(rotation_vec):
        '''(
        Generate the rotation matrix to rotate a vector to the z direction
        
        Args:
            rotation_vec:
                the vector to be rotated
            translation_vec:
                Optional translation vector to move the vector before the rotation
        '''
        deg1=Adsorbate.rotation_angle(rotation_vec[0], rotation_vec[1])*(-1)
        to_xz = SymmOp.from_axis_angle_and_translation((0,0,1), deg1, angle_in_radians=True)              
        length = (rotation_vec[0]**2+rotation_vec[1]**2+rotation_vec[2]**2)**0.5
        deg2=Adsorbate.rotation_angle(rotation_vec[2],len=length)*(-1)
        to_z = SymmOp.from_axis_angle_and_translation((0,1,0), deg2, angle_in_radians=True)
        return np.dot(to_z,to_xz)
    
   
    @staticmethod
    def length(vec):
        '''
        Get the length of vector
        '''
        return linalg.norm(vec)
    
    @staticmethod
    def unitvec(vec):
        '''
        Get the unit vector of a vector
        '''
        return vec/Adsorbate.length(vec) 
    
    @staticmethod
    def vector_angle(vec1,vec2):
        '''
        Get the angle between two vectors
        '''
        return np.arccos(np.dot(vec1,vec2)/(Adsorbate.length(vec1)*Adsorbate.length(vec2)))
    
    @staticmethod
    def gen_vec(p1,p2):
        '''
        Get the vector connecting two points
        '''
        return np.array([p1[i]-p2[i] for i in range(3)])
    
    def initialguess(self, points):
        return np.average(np.array([Adsorbate.unitvec(points[i]) for i in points]), axis=0)
                                   
    def angelsum(self, var, siteindex, vectopoint):
        '''
        Target function to get the sum of angles between lattice vectors to the vector
        '''
        anglelist=[]
        vectar=Adsorbate.gen_vec(var,self.sites[siteindex].coords)
        for i in range(len(self.sites)):
            if i != siteindex:
                anglelist.append(Adsorbate.vector_angle(vectopoint[i],vectar))
        return np.sum(np.array(anglelist)**2)
        
    def constr(self, var, *args):
        return Adsorbate.length(var)-1  
                 
    def molecule_adosrbate_optimizer(self, siteindex, accu=1e-5, eps=1e-4):
        '''
        Find the best direction to rotate the molecule for a surface
        Args:
            siteindex: the index of the site to be closest to the surface
        '''
        #generate vector with the reference point in molecule
        vectopoint=dict()
        
        for i in range(len(self.sites)):
            if i !=siteindex:
                vectopoint[i]= Adsorbate.gen_vec(self.sites[i].coords, self.sites[siteindex].coords)

        init = self.initialguess(vectopoint)
        print init                   
        try:
            res=fmin_slsqp(self.angelsum, init, args=[siteindex, vectopoint], iter=1e5, 
                           f_eqcons=self.constr, iprint=2, acc=accu, epsilon=eps, full_output=True)
            print res
            if res[3] != 0:
                raise Exception(res[-1])
        except Exception as e:
            print "Error optimizing the structure"+str(e)       
        else:
            return res