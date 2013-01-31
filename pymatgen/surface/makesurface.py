from collections import defaultdict
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.command_line.aconvasp_caller import run_aconvasp_command
from pymatgen.io.vaspio.vasp_input import Poscar
from pymatgen.command_line.gulp_caller import get_binoxi_gulp_energy
from pymatgen.core.structure import Structure
from pymatgen.surface import surftool as st
from pymatgen.core.structure import SiteCollection, Structure
from pymatgen.core.lattice import Lattice
from pymatgen.surface import surface
from pymatgen.surface import surfacebasistransform as bb
from pymatgen.surface.surfacebasistransform import Sobat
from pymatgen.core.sites import PeriodicSite

import numpy as np
import copy
import numpy.linalg as npl
from numpy import pi
import fractions
import math 
import time

'''
Created on Dec 20, 2012

@author: WenhaoMIT
'''


class VaspSurface():
    def StandardSurfDict(self,bulkstruct,maxindex,slablength=15,vacuumlength=10,utrfix=0):
        '''Create a dictionary of surfaces - key is index, value is a list of surfaces with that index'''
        '''Eliminate polar, duplicate surfaces, return primitive slabs'''
        
        listofsurfs=self.makesurface(bulkstruct,maxindex,slablength,vacuumlength,utrfix)
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
                if all(surf.indices[a]==I[a] for a in range(0,3)):
                    indexlistofsurfs.append(surf)
            strI="["+str(int(I[0]))+","+str(int(I[1]))+","+str(int(I[2]))+"]"
            surfdict[strI]=indexlistofsurfs
            
        '''Eliminate polar and duplicate surfaces in those lists'''
            
        for key in surfdict.keys():
            print key
            listofsurfs=surfdict[key]
            listnopolarsurfs=self.checkpolar(listofsurfs)
            print len(listnopolarsurfs)
            listnoduplicates=self.checksame(listnopolarsurfs)
            print len(listnoduplicates)
            '''Make primitive slabs'''
            primitives=self.make_primitive_slab(listnoduplicates)
            surfdict[key]=primitives
        return surfdict

    
    def makesurface(self, bulkstruct,maxindex,slablength,vacuumlength,utrfix):
        '''Makes a list of surfaces from a bulk structure, given the initial parameters
        
        bulkstruct=Bulk structure. Should be conventional unit cell!! Can use surftools.SurfTools().getconventional() for this
        maxindex=np.array of index, [1, 1, 1] for example
        slablength: If slablength is positive, that's the minimum length of slab in angstroms. If negative, that's number of layers desired
        vacuumlength: Vacuumlength
        utrfix: Bad variable name, will change later. Basically, fix termination or not? Fix = 1, enumerate all = 0
        
        '''
        
        ListofSurfs=[]
        if vacuumlength==0:
            seldyn=0
            utrfix=1
            shouldcheckpolar=0
        
        atoms=bulkstruct.sites
        basis=bulkstruct.lattice.matrix
        
        atomcoords=np.array(bulkstruct.frac_coords)

        
        sconv=st.Surftool().getconventional(bulkstruct)
        pg=st.Surftool().getpg(sconv)
        
        origB=basis
        listofSobats=bb.SurfaceOrientedBasisTransform().transformbasis(sconv, maxindex, pg)
                
        for sobat in listofSobats:
            oldsobat=copy.deepcopy(sobat)
            atoms=sobat.atomcoords
            uniqterm=np.unique(np.around(atoms[:,2],decimals=5))
            utr=len(uniqterm)
            
            if utrfix!=0:
                utr=1
                
            basis=sobat.standardbasis
            
            if slablength > 0:
                numcells=np.ceil(slablength/npl.norm((basis[2])))
            elif slablength < 0:
                numcells=-slablength
            elif slablength==0:
                numcells=1
            basis[2]=numcells*basis[2]
                
            primvol=np.abs(np.dot(origB[2],np.cross(origB[0],origB[1])))
            newvol=sobat.newbasisvolume()
            volmulti=round(newvol/primvol)
            
            offang=st.Surftool().ang(basis[2], np.cross(basis[0],basis[1]))
            vacuumzproject=vacuumlength/math.cos(offang)
            vacbasisvect=basis[2]+vacuumzproject*basis[2] / npl.norm(basis[2])
            
            dilationfactor=npl.norm(basis[2])/npl.norm(vacbasisvect)
            newbasis=basis
            newbasis[2]=vacbasisvect
            
            for utii in range(0,utr):
                
                species=[]
                surfatoms=[]
                
                termatoms=atoms
                for elem in range(0,len(termatoms)):
                    termatoms[elem][2]=(termatoms[elem][2]-uniqterm[utii])%1
                ext=2
                numcells=int(numcells)
                
                for NN in range(0,numcells):
                    extatoms=termatoms
                    newextatoms=np.zeros((len(extatoms),4))
                    fltnumcells=float(numcells)
                    fltNN=float(NN)
                    for elem in range(0,len(extatoms)):
                        newextatoms[elem][0]=(extatoms[elem][0])
                        newextatoms[elem][1]=(extatoms[elem][1])
                        newextatoms[elem][ext]=(extatoms[elem][ext]+fltNN)/(fltnumcells);
                        newextatoms[elem][3]=(extatoms[elem][3])
                    surfatoms.append(newextatoms)
                    
                onespecies=bulkstruct.species
                for MM in range(0,len(onespecies)):
                    for NN in range(0,int(volmulti)):                    
                        species.append(onespecies[MM])
                
                cc=0
                tmpsurfatoms=np.zeros((len(surfatoms)*len(surfatoms[0]),4))
                for PAx in range(0,len(surfatoms)):
                    for PAy in range(0,len(surfatoms[PAx])):
                        tmpsurfatoms[cc]=surfatoms[PAx][PAy]
                        cc=cc+1
                
                surfatoms=tmpsurfatoms
                surfatoms=surfatoms[surfatoms[:,3].argsort(),]
                
                for elem in range(0,len(surfatoms)):
                    surfatoms[elem][2]=surfatoms[elem][2]*dilationfactor
                
                surfatoms=np.round(surfatoms,6)
                surfatoms=surfatoms[:,0:3]
                Surf=surface.Surface(newbasis, species, surfatoms,sobat.millerindex,vacuumlength,numcells,oldsobat)
                
                ListofSurfs.append(Surf)
        
        return ListofSurfs
    
    
    def make_primitive_slab(self,listofslabs):
        ''' For a list of slabs, makes a list of primitive slabs '''
        
        listofprimsurfs=[]
        for slab in listofslabs:
            primsurf=run_aconvasp_command(["aconvasp", "--prim"],slab)
            primslab_string = ""
            for line in primsurf[0].split("\n"):
                primslab_string = primslab_string + line + "\n"
            primslabstruct=Poscar.from_string(primslab_string).struct
            
            primslabstruct_coords=np.array(primslabstruct.frac_coords)
            
            
            #for atom in primslabstruct.sites:
            #    primslabstruct_coords.append(atom.frac_coords)
            
            '''bulk
            %   You don't want to do bulk
            %   Because of the terminations, this doesn't work well at all
            %   The safest way is probably to just use the bulk unit cell generated from surfacebasistransform
            %   Even if the surface unit cell is not primitive, it won't affect the KPOINTS calculation
            %   Because the symmetry should be the same, and the k-points sampling will be fine.
            %   The following was a waste of an afternoon
            
                
                
            novacv3=(primslabstruct.lattice.c-slab.vacuum)
            extensionfactor=primslabstruct.lattice.c/novacv3
            
            surfatoms=[]
            for atom in primslabstruct.sites: 
                surfatoms.append(atom.frac_coords)
            
            for elem in range(0,len(surfatoms)):
                    surfatoms[elem][2]=surfatoms[elem][2]*extensionfactor
            newv3=novacv3/float(slab.numlayers)
            maxatomz=newv3/novacv3
            bulkatoms=[]
            for atom in surfatoms:
                if atom[2]<maxatomz-0.003:
                    atom[2]=atom[2]/maxatomz
                    bulkatoms.append(atom)
            
            npbulkatoms=np.zeros((len(bulkatoms),3))
            for ii in range(0,len(bulkatoms)):
                npbulkatoms[ii]=bulkatoms[ii]
                
            if slab.num_sites % len(bulkatoms) != 0:
                raise SystemError("Tolerance should be smaller")
            
            multi=slab.num_sites/len(bulkatoms)
            
            elamt=slab.composition.get_el_amt_dict()
            
            for el in elamt:
                elamt[el]=elamt[el]/multi
                
            bulkspecies=[]
            for el in elamt:
                for ii in range(0,int(elamt[el])):
                    bulkspecies.append(el)
            
            Bulklattice=np.array(primslabstruct.lattice.matrix)
            Bulklattice[2]=Bulklattice[2]*newv3/npl.norm(Bulklattice[2])
            primbulk=Structure(Bulklattice,bulkspecies,bulkatoms)
            newsobat=bb.SurfaceOrientedBasisTransform().sobatfromprimbulk(primbulk, slab.indices)
            
            primslab=surface.Surface(primslabstruct.lattice,primslabstruct.species,primslabstruct_coords,slab.indices,slab.vacuum,slab.numlayers,newsobat)
            '''
            primslab=surface.Surface(primslabstruct.lattice,primslabstruct.species,primslabstruct_coords,slab.indices,slab.vacuum,slab.numlayers,slab.sobat)
            listofprimsurfs.append(primslab)
            
        return listofprimsurfs
                
    
    
    def checkpolar(self,listofsurfs):
        ''' For a list of slabs, returns non-polar slabs. 
        Checks for oxidations states using bond valence sum analyzer'''
        
        listnopolarsurfs=[]
        oxidations=BVAnalyzer().get_valences(listofsurfs[0])
        for surf in listofsurfs:
            B=surf.lattice.matrix
            species=surf.species
            A=np.array(surf.frac_coords)
            
            zcenter=B[2]*(np.max(A[:,2])+np.min(A[:,2]))/2
            normal=np.cross(B[0],B[1])/npl.norm(np.cross(B[0], B[1]))
            center=st.Surftool().proj(normal,zcenter)
            
            DP=0
            for ii in range(0,len(A)):
                r=A[ii]
                q=oxidations[ii]
                DP=DP+q*(st.Surftool().proj(normal,r)-center)
            if npl.norm(DP)<1E-3:
                listnopolarsurfs.append(surf)
        return listnopolarsurfs
            
            
    def binoxi_gulp_surface_energy(self,slab):
        '''Given a slab, returns the surface energy calculated from a Tersoff Potential via GULP'''
        
        bulk=slab.sobat.structure()
        Ebulk=get_binoxi_gulp_energy(bulk)
        N=float(slab.num_sites)/float(bulk.num_sites)
        Eslab=get_binoxi_gulp_energy(slab)
        lattice=slab.lattice.matrix
        A=npl.norm(np.cross(lattice[0],lattice[1]))
        SE=16.02177*(Eslab-N*Ebulk)/(2*A)
        return SE
            
    def checksame(self,listofslabs):
        '''Eliminates duplicate slabs from a list of slabs'''
        listduplicates=[]
        Eslabs=[]
        for slab in listofslabs:
            Eslabs.append(get_binoxi_gulp_energy(slab))
        
        for bbx in range(0,len(listofslabs)-1):
            for bby in range(bbx+1,len(listofslabs)):
                sameyn=False
                
                Eone=Eslabs[bbx]
                Etwo=Eslabs[bby]
                
                if np.abs(Eone-Etwo)<0.00001:
                    slabone=listofslabs[bbx]
                    slabtwo=listofslabs[bby]
                    sameyn=StructureMatcher().fit(slabone,slabtwo)
                    
                if sameyn==True:
                    listduplicates.append(bby)
        listuniquesurfs=[]
        for a in range(0,len(listofslabs)):
            if a not in listduplicates:
                listuniquesurfs.append(listofslabs[a])
        return listuniquesurfs