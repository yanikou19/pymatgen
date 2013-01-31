'''
Created on Dec 1, 2012

@author: WenhaoMIT
'''
from pymatgen.io.vaspio.vasp_input import Poscar
from pymatgen.command_line.aconvasp_caller import run_aconvasp_command
from pymatgen.surface import surftool as st
from pymatgen.core.structure import SiteCollection, Structure
from pymatgen.core.lattice import Lattice
import numpy as np
import numpy.linalg as npl
from numpy import pi
import fractions
import math 
import time

class Sobat(object):
    def __init__(self, threepts,twovects,v3,standardbasis, millerindex,atomcoords,species,bulkstructure):
        self.threepts=threepts
        self.twovects=twovects
        self.v3=v3
        self.standardbasis=standardbasis
        self.millerindex=millerindex
        self.atomcoords=atomcoords
        '''bulkstructure should be conventional cell'''
        self.species=species
        self.bulkstructure=bulkstructure
        
    def display(self):
        print "*************************************************************************"
        print "Miller Index: "
        print self.millerindex
        print "\n Three Points (w.r.t Original Cell)\n"
        print self.threepts
        print "\n Basal Vectors v1 and v2 (w.r.t Original Cell) \n "
        print self.twovects
        print "\n Maximally Orthogonal v3 (w.r.t Original Cell) \n "
        print self.v3
        print "\n POSCAR BASIS (w.r.t Standard Basis Gram Schmidt around v1\n"
        print self.standardbasis
        print "\n New Basis Volume \n"
        print self.newbasisvolume()
        print "\n Covariantly Transformed Atomic Coordinates \n"
        print self.atomcoords
        print
        print
    
    def newbasisvolume(self):
        POSBAS=self.standardbasis
        POSBASVOL=np.abs(np.dot(POSBAS[0],np.cross(POSBAS[1],POSBAS[2])))
        return POSBASVOL
    
    def v1xv2(self):
        return np.cross(self.twovects[0],self.twovects[1])
    
    def structure(self):
        lattice = self.standardbasis
        ii=0
        atomcoordsnoindex=np.zeros((len(self.atomcoords),3))
        for coord in self.atomcoords:
            atomcoordsnoindex[ii]=coord[0:3]
            ii=ii+1
        return Structure(lattice, self.species,atomcoordsnoindex)
    
    '''
    Another failed attempt to get the primitive bulk in 2D
        
    def getprimbulk(self):
        lattice = self.standardbasis
        lattice[2]=lattice[2]*2
        ii=0
        atomcoordsnoindex=np.zeros((len(self.atomcoords),3))
        for coord in self.atomcoords:
            coord[2]=coord[2]*0.5
            atomcoordsnoindex[ii]=coord[0:3]
            ii=ii+1
        bulkstruct=Structure(lattice, self.species,atomcoordsnoindex)
        
        primsurf=run_aconvasp_command(["aconvasp", "--prim"],bulkstruct)
        primslab_string = ""
        for line in primsurf[0].split("\n"):
            primslab_string = primslab_string + line + "\n"
        primslabstruct=Poscar.from_string(primslab_string).struct
        
        primbulkatoms=[]
        for atom in primslabstruct.sites:
            primbulkatoms.append(atom.frac_coords)
        
        lattice = primslabstruct.lattice.matrix
        lattice[2]=lattice[2]*0.5
        ii=0
        atomcoordsnoindex=np.zeros((len(primbulkatoms),3))
        for coord in primbulkatoms:
            coord[2]=coord[2]*2
            atomcoordsnoindex[ii]=coord[0:3]
            ii=ii+1
        return Structure(lattice, primslabstruct.species,atomcoordsnoindex)
        '''
        
        
class SurfaceOrientedBasisTransform():
    
    def sobatfromprimbulk(self,structure,maxindex):
        '''From a given structure, return a sobat'''
        
        L=np.array(structure.lattice.matrix)
        P=np.zeros((3,3))
        P[1]=L[0]
        P[2]=L[1]
        twovects=np.zeros((2,3))
        twovects[0]=L[0]
        twovects[1]=L[1]
        atomcoords=[]
        for site in structure.sites:
            atomcoords.append(site.frac_coords)
        newsobat=Sobat(P,twovects,L[2],L,maxindex,atomcoords,structure.species,structure)
        return newsobat
        
    
    def transformbasis(self, convstruct, maxindex, pg):
        '''Reorients a CONVENTIONAL unit cell such that its basal (001) plane lies in the miller index of desire
        Requires a symmetry file called pg. Can be obtained by the function getpg in surftools
        Maxindex is a np.array. Can be a Miller Index (np.array([1,1,1]), or the dimensions of a supercell to find surfaces in. 
            For example, maxindex=np.array([2]) creates a 2x2x2 supercell of bravais points, and finds all surfaces possible in there.
        '''
        
        
        
        atomsites=convstruct.sites
        basis=convstruct.lattice.matrix
        
        atoms=[]
        for atom in atomsites: 
            atoms.append(atom.frac_coords)
        
        ListofSobats=[]
        if len(maxindex)==1:
            millera=maxindex
            millerb=maxindex
            millerc=maxindex
            vector=0
        elif len(maxindex)==3:
            vector=1
            twovects,P=st.Surftool().get2vectsinplane(basis, maxindex)
            v3=np.array([0,0,0])
            standardbasis=np.eye(3,3)
            atomcoords=np.array([])
            species=[]
            S3A2V=Sobat(P,twovects,v3,standardbasis,maxindex,atomcoords,species,convstruct)
            ListofSobats.append(S3A2V)
        else: 
            raise StandardError('Warning, maxindex must be a number or a vector, exiting')
            return
        
        rr=1
        primvol=np.abs(np.dot(basis[0],(np.cross(basis[1],basis[2]))))
        '''tic'''
        
        if len(maxindex)==1:
            
            bravaissup=[]
            for ii in range(0,millera):
                for jj in range(0,millerb):
                    for kk in range(0,millerc):
                        x=basis[0]*ii
                        y=basis[1]*jj
                        z=basis[2]*kk
                        v=x+y+z
                        bravaissup.append(v)
            
            M=len(bravaissup)
            POS=[]
            v1xv2check=[]
            millercheck=[]
            counter=0;
            for ii in range(1,M+1):
                for jj in range (ii+1,M+1):
                    for kk in range(jj+1,M+1):
                        ncr=[ii ,jj ,kk]
                        POS.append(np.array(ncr))
        
        
            for ii in range(0,len(POS)):
                
                Setof3Atoms=np.array([bravaissup[POS[ii][0]-1], bravaissup[POS[ii][1]-1], bravaissup[POS[ii][2]-1]])
                
                for aa in range(0,3):
                    
                    v1=Setof3Atoms[np.mod(aa+1,3)]-Setof3Atoms[np.mod(aa,3)]
                    v2=Setof3Atoms[np.mod(aa+2,3)]-Setof3Atoms[np.mod(aa,3)]
                    T=st.Surftool().ang(v1, v2)
                    if T <= pi/2:
                        bb=aa
                        break
                    
                threepts=np.array([Setof3Atoms[np.mod(bb,3)], Setof3Atoms[np.mod(bb+1,3)], Setof3Atoms[np.mod(bb+2,3)]])    
                v1xv2=np.cross(v1,v2)
                
                if npl.norm(v1xv2)==0:
                    continue
                
                if v1xv2[0] < 0:
                    vtemp=v1
                    v1=v2
                    v2=vtemp
                    v1xv2=v1xv2*-1
                    toswitch=1
                elif v1xv2[0]==0:
                    if v1xv2[1]<0:
                        vtemp=v1
                        v1=v2
                        v2=vtemp
                        v1xv2=v1xv2*-1
                        toswitch=1
                    elif v1xv2[1]==0:
                        if v1xv2[2]<0:
                            vtemp=v1
                            v1=v2
                            v2=vtemp
                            v1xv2=v1xv2*-1
                            toswitch=1
                else:
                    continue
                
                
                millerindex=st.Surftool().getmillerfrom2v(basis, v1, v2)
                C=st.Surftool().getsymfam(basis,millerindex,pg)
                millerindex=C[-1]
                
                duplicate=0
                replace=0
                
                for mm in range(0,len(millercheck)):
                    if (millercheck[mm]==millerindex).all():
                        duplicate=1
                        if npl.norm(v1xv2)<npl.norm(v1xv2check[mm]):
                            replace=mm
                            v1xv2check[mm]=v1xv2
                            millercheck[mm]=millerindex
                
                if vector==0:
                    if duplicate==0:
                        v1xv2check.append(v1xv2)
                        millercheck.append(millerindex)
                        counter=counter+1
                        twovects=np.array([v1,v2])
                        v3=[0,0,0]
                        standardbasis=np.eye(3,3)
                        atomiccoords=np.array([[0,0,0,0]])
                        species=[]
                        S3A2V=Sobat(threepts,twovects,v3,standardbasis,millerindex,atomiccoords,species,convstruct)
                        ListofSobats.append(S3A2V)
                        
                    elif duplicate==1 and replace !=0:
                        atomiccoords=np.array([[0,0,0,0]])
                        S3A2V=Sobat(threepts,twovects,v3,standardbasis,millerindex,atomiccoords,species,convstruct)
                        ListofSobats[replace]=S3A2V
                        
        
        '''
        Check if this part is correct:            
        for jj in range(0,len(ListofSobats)):
            X=ListofSobats[jj]
            X.display();
        '''
                    
        '''
        Find third vector
        '''
        
        counter=len(ListofSobats)
        
        x=basis[0]
        y=basis[1]
        z=basis[2]
        
        for ss in range(0,counter):
            S=ListofSobats[ss]
            lastT=100
            S7=np.cross(S.twovects[0],S.twovects[1])
            v1=S.twovects[0]
            v2=S.twovects[1]
            for ii in range(-1,2):
                for jj in range(-1,2):
                    for kk in range(-1,2):
                        if ii==0 and jj==0 and kk==0:
                            continue
                        
                        v3temp=ii*x+jj*y+kk*z
#                       print v3temp
                        T=st.Surftool().ang(S7, v3temp)
                        Tv1=st.Surftool().ang(v1,v3temp)
                        Tv2=st.Surftool().ang(v2,v3temp)
                        if T<= pi/2 and Tv1 <= pi/2 and Tv2 <= pi/2 and Tv1 > 0 and Tv2 >0:
                            if T<lastT:
                                lastT=T
                                v3=v3temp
            
            ListofSobats[ss].v3=v3
#            ListofSobats[ss].display()

        
        '''
        Graham-Schmidt Transform Standard Basis
        '''
        
        for tt in range(0,counter):
            S=ListofSobats[tt]
            v1=S.twovects[0]
            v2=S.twovects[1]
            v3=S.v3

            u1=v1
            u2=v2-st.Surftool().proj(u1,v2)
            u3=v3-st.Surftool().proj(u1,v3)-st.Surftool().proj(u2,v3)
            
            b1=v1
            b2=v2
            b3=v3
            f1hat=u1/npl.norm(u1)
            f2hat=u2/npl.norm(u2)
            f3hat=u3/npl.norm(u3)
            
            POSBAS=np.array([[st.Surftool().sproj(f1hat,b1),0,0],
                 [st.Surftool().sproj(f1hat,b2),st.Surftool().sproj(f2hat,b2),0],
                 [st.Surftool().sproj(f1hat,b3),st.Surftool().sproj(f2hat,b3),st.Surftool().sproj(f3hat,b3)]])
            
            
            if np.abs(st.Surftool().ang(f1hat,f2hat)-pi/2)>=0.0001 or np.abs(st.Surftool().ang(f3hat,f2hat)-pi/2)>=0.0001 or np.abs(st.Surftool().ang(f1hat,f3hat)-pi/2)>=0.0001:
                raise SystemError('ERROR: Gram-schmidt Failure') 
        
            if st.Surftool().sproj(f2hat,b1) >= 0.001 or st.Surftool().sproj(f3hat,b1)>0.001 or st.Surftool().sproj(f3hat,b2)>0.001:
                raise SystemError('ERROR: Non-bottom diagonalized POSCAR matrix')
            
            POSBAS=np.array([[st.Surftool().sproj(f1hat,b1),0,0],
                             [st.Surftool().sproj(f1hat,b2),st.Surftool().sproj(f2hat,b2),0],
                             [st.Surftool().sproj(f1hat,b3),st.Surftool().sproj(f2hat,b3),st.Surftool().sproj(f3hat,b3)]])   
            
            
            S.standardbasis=POSBAS
            
            Ta=st.Surftool().ang(b1,b2)-st.Surftool().ang(POSBAS[0],POSBAS[1])
            Tb=st.Surftool().ang(b2,b3)-st.Surftool().ang(POSBAS[1],POSBAS[2])
            Tc=st.Surftool().ang(b1,b3)-st.Surftool().ang(POSBAS[0],POSBAS[2])
            
            if Ta>=0.001 or Tb>=0.001 or Tc>=0.001:
                raise SystemError('ERROR: Angles in POSCAR inconsistent with Basis!')
            
            La=npl.norm(b1)-npl.norm(POSBAS[0])
            Lb=npl.norm(b2)-npl.norm(POSBAS[1])
            Lc=npl.norm(b3)-npl.norm(POSBAS[2])
            
            if La>0.0001 or Lb>=0.0001 or Lc>=0.0001:
                raise SystemError('ERROR: Lengths in POSCAR inconsistent with Basis!')
            
            
        ListofAtoms=[]
        
        for mm in range(0,counter):
            S=ListofSobats[mm]
            Sxyz=np.array([[0,0,0],[0,0,0]])
            Ra=len(atoms)
            O=S.threepts[0]
            
            NewB=np.zeros((3,3))
            NewB[0]=S.twovects[0]
            NewB[1]=S.twovects[1]
            NewB[2]=S.v3
            
            C=np.dot(NewB,npl.inv(basis))
            
            v1=NewB[0]
            v2=NewB[1]
            v3=NewB[2]
            
            millerv1=np.dot(npl.inv((np.transpose(basis))),np.transpose(v1))
            millerv2=np.dot(npl.inv((np.transpose(basis))),np.transpose(v2))
            millerv3=np.dot(npl.inv((np.transpose(basis))),np.transpose(v3))
            
            for xx in range(0,2):
                for yy in range(0,2):
                    for zz in range(0,2):
                        M=millerv1*xx+millerv2*yy+millerv3*zz
                        for aa in range(0,3):
                            if M[aa]>Sxyz[0][aa]:
                                Sxyz[0][aa]=round(M[aa])
                            if M[aa]<Sxyz[1][aa]:
                                Sxyz[1][aa]=round(M[aa])
            
            R=len(atoms)
            labelatoms=np.zeros((R,4))
            for ii in range(0,R):
                label=np.array([atoms[ii][0],atoms[ii][1],atoms[ii][2], ii])              
                labelatoms[ii]=label
            
            NewAtomsList=[]
            for ii in range(0,R):
                for x in range(Sxyz[1][0],Sxyz[0][0]+1):
                    for y in range(Sxyz[1][1],Sxyz[0][1]+1):
                        for z in range(Sxyz[1][2],Sxyz[0][2]+1):
                            zz=labelatoms[ii]+[x,0,0,0]+[0,y,0,0]+[0,0,z,0]
                            NewAtomsList.append(zz)
            
            NewAtoms=np.zeros((len(NewAtomsList),4))
            for NAx in range(0,len(NewAtomsList)):
                NewAtoms[NAx]=NewAtomsList[NAx]
                
            
            C=np.dot(NewB,npl.inv(basis))
            C=np.append(C,[[0],[0],[0]],1)
            C=np.append(C,[[0,0,0,1]],0)
            
            POSATOMS=[]
            Rr=len(NewAtoms)
            NewAtoms=np.transpose(np.dot(npl.inv(np.transpose(C)),np.transpose(NewAtoms)))
            
#            print NewAtoms
            
            '''
            Remove Excess Atoms
            '''
            for Nz in range(0,Rr):
                XA=NewAtoms[Nz][0]
                XB=NewAtoms[Nz][1]
                XC=NewAtoms[Nz][2]
                XL=NewAtoms[Nz][3]
                if XA>=-0.001 and XA<0.999 and XB>=-0.001 and XB<0.999 and XC>=-0.001 and XC<0.999:
                    newatom=np.array([[XA,XB,XC,XL]])
                    POSATOMS.append(newatom)
            
            AtomicCoords=np.zeros((len(POSATOMS),4))
            for PAx in range(0,len(POSATOMS)):
                AtomicCoords[PAx]=POSATOMS[PAx]
            
            POSBASVOL=S.newbasisvolume()/primvol
            RPA=len(AtomicCoords)

            if np.abs(RPA/Ra-POSBASVOL)>=0.01:
                raise SystemError("ERROR: Number of atoms inconsistent with the size of basis!")
            
            S.atomcoords=AtomicCoords
            species=[]
            volmulti=round(S.newbasisvolume()/primvol)
            onespecies=convstruct.species
            for MM in range(0,len(onespecies)):
                for NN in range(0,int(volmulti)):                    
                    species.append(onespecies[MM])
            S.species=species
        
        return ListofSobats
        
        
        
        