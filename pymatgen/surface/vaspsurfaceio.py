'''
Created on Jan 27, 2013

@author: wenhao
'''
import os
import numpy as np

from pymatgen.io.vaspio.vasp_input import Incar, Poscar, Kpoints, Potcar,PotcarSingle
from pymatgen import Composition, Structure, zopen
from pymatgen.io.io_abc import VaspInput
from pymatgen.io.vaspio_set import MITVaspInputSet, MITHSEVaspInputSet, \
    MaterialsProjectVaspInputSet, MITGGAVaspInputSet
    

class CreateSurfaceInputs():
    def __init__(self, name, homedirectory, slabdict,relaxlayers=3,zonly=True):
        '''
        Create a directory for a surface dictionary
        Folder name: Structure name
        Folder for each index
        Each index:
        - 1 Bulk calculation
        - All slabs: 
          -With and without relaxation
        Info files for each calculation
        
        Implements the "k-points algorithm" for rapid convergence. 
        '''
        
        surfacedirectory=homedirectory+"/"+name
        

        if not os.path.exists(surfacedirectory):
            os.makedirs(surfacedirectory)
            
        for index in slabdict.keys():
            listofsurfs=slabdict[index]
            if len(listofsurfs) > 0:
                indexdirectory=surfacedirectory+"/"+index
                if os.path.exists(indexdirectory):
                    raise SystemError('Index directory exists')
                os.makedirs(indexdirectory)
                
                listofsurfs=slabdict[index]
                bulkdirectory=indexdirectory+"/BULK/"
                os.makedirs(bulkdirectory)
                
                '''Create Bulk input files'''
                
                bulk=listofsurfs[0].sobat.structure()
                poscar=Poscar(bulk).write_file(bulkdirectory+"POSCAR")
                Kpts=Kpoints().gamma_automatic((7,7,7), (0,0,0))
                Kpts.write_file(bulkdirectory+"KPOINTS")
                INCAR=SurfaceIO().surfincar(bulk,"static")
                INCAR.write_file(bulkdirectory+"INCAR")
                surfparamset = MaterialsProjectVaspInputSet()
                #POTCAR = surfparamset.get_potcar(bulk)
                #POTCAR.write_file(bulkdirectory+"POTCAR")
                
                for slabnum in range(0,len(listofsurfs)):
                    slabdirectory=indexdirectory+"/slab"+str(slabnum)
                    os.makedirs(slabdirectory)
                    slab=listofsurfs[slabnum]
                    for relaxyn in ["relax", "static"]:
                        slabrelaxdirectory=slabdirectory+"/"+relaxyn
                        os.makedirs(slabrelaxdirectory)
                        if relaxyn=="relax":
                            seldyn=SurfaceIO().relaxslab(slab, relaxlayers, zonly)
                            poscar=Poscar(slab,comment=None,selective_dynamics=seldyn).write_file(slabrelaxdirectory+"/POSCAR")
                        else:
                            poscar=Poscar(slab).write_file(slabrelaxdirectory+"/POSCAR")
                        Kpts=Kpoints().gamma_automatic((7,7,1), (0,0,0))
                        Kpts.write_file(slabrelaxdirectory+"/KPOINTS")
                        INCAR=SurfaceIO().surfincar(slab,relaxyn)
                        INCAR.write_file(slabrelaxdirectory+"/INCAR")
                        #POTCAR = surfparamset.get_potcar(bulk)
                        #POTCAR.write_file(slabrelaxdirectory+"POTCAR")
                    

class SurfaceIO():
    def relaxslab(self,slab,layers,zonly):
        '''
        Create a selective dynamics POSCAR to relax the outermost layers of a surface:
        
        For N<0, that is the number of layers to relax
        For N>0, that is the 'distance' into the slab to relax
        
        zonly = set whether or not relaxation is vertical only or in all directions
        '''
        layers=float(layers)
        atoms=np.array(slab.frac_coords)
        uniqterm=np.unique(np.around(atoms[:,2],decimals=5))
        seldyns=[]
        if layers>0:
            relaxlen=layers/slab.lattice.c
            maxz=np.max(uniqterm)-relaxlen
            minz=np.min(uniqterm)+relaxlen
        elif layers<0:
            maxz=uniqterm[layers]
            minz=uniqterm[-layers-1]
        else: 
            return 
        
        for site in slab.sites:
            if site.frac_coords[2] <= maxz and site.frac_coords[2] >= minz:
                seldyn=[False, False, False]
            else:
                if zonly==True:
                    seldyn=[False, False, True]
                else:
                    seldyn=[True, True, True]
            seldyns.append(seldyn)
        return seldyns    
        
    def surfincar(self,slab,relax):
        '''Ideally, this should just take a slab and return a reasonable surface INCAR object for the slab.
        
        Major things:
        ISIF=2
        ISMEAR=1
        SIGMA=0.05
        ^^ The SIGMA=0.05 eV value is suggested from Da Silva et al, Surface Science (2006)
        Lowest SIGMA is better, but too low challenges self-consistent calculations. If cannot converge, then increase SIGMA
        Varying SIGMA does not change the surface energy or work function much, but does affect the interlayer spacings a lot. 
        
        If relaxation is on
        IBRION=2
        NSW=51
        SYMPREC=1E-8
        
        For static, 
        IBRION=-1
        SYMPREC=1E-6
        NSW=0
    
        For Oxides:
        LDA+U as necessary
    
        For Magnetic:
        Return an INCAR with reasonable Magnetic states
        '''
        surfparamset = MaterialsProjectVaspInputSet()
        incar = surfparamset.get_incar(slab)
        incar["ISIF"] = "2"
        incar["ISMEAR"] = "1"
        incar["SIGMA"] = "0.05"
        
        if relax=="relax":
            incar["IBRION"] = "2"
            incar["SYMPREC"] = "1E-8"
            incar["NSW"] = "51"
        elif relax=="static": 
            incar["IBRION"] = "-1"
            incar["SYMPREC"] = "1E-6"
            incar["NSW"] = "0"
            
        return incar
    
    def slabkpoints(self):
        '''
        Since there is 1 kpoint in the z direction, always use a Gamma-point grid
        
        The size in the other two direction is based on size
        
        Distribution ... doesn't matter too much, as long as slab and bulk are the same.
        Let's just go 7x7x1 on all slabs and 7x7x7 on all bulks for now ...
        No huge difference in cost and will be safest. 

        '''
        kpts=Kpoints().gamma_automatic((7,7,1), (0,0,0))
        return kpts
    
    def bulksurfkpts(self):
        '''
        Since there is 1 kpoint in the z direction, always use a Gamma-point grid
        
        The size in the other two direction is based on size
        
        Distribution ... doesn't matter too much, as long as slab and bulk are the same.
        Let's just go 7x7x1 on all slabs and 7x7x7 on all bulks for now ...
        No huge difference in cost and will be safest. 

        '''
        kpts=Kpoints().gamma_automatic((7,7,7), (0,0,0))
        return kpts
        