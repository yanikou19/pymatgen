'''
Created on Nov 30, 2012

@author: WenhaoMIT
'''

import numpy as np
import numpy.linalg as npl
from numpy import pi
import fractions
import math 
from pymatgen.command_line.aconvasp_caller import run_aconvasp_command
from pymatgen.io.vaspio.vasp_input import Poscar
from pymatgen.io.cifio import CifParser
from pymatgen.symmetry.finder import SymmetryFinder
from pymatgen.core.structure import Lattice, Structure
from pymatgen.core.periodic_table import Element

class Surftool():
    
    def __init__(self):
        '''A host of tools important for creating surfaces'''
        
    def ang(self, v1, v2):
        if npl.norm(v1)==0 or npl.norm(v2)==0:
            raise StandardError('One of your vectors has length zero')
        x=np.dot(v1,v2)/(npl.norm(v1)*npl.norm(v2))
        a= np.float('%.12f' %x)
        angle=np.arccos(a)
        return angle
    
    def get2vectsinplane(self,basis, maxindex):
        #% The two returned vectors are in cartesian coordinates
        if len(maxindex) != 3:
            raise StandardError('Error need 3 elements in your maxindex')
    
        zeroind=np.where(maxindex==0)[0]
        numzero = len(zeroind)
        
        h = float(maxindex[0])
        k = float(maxindex[1])
        l = float(maxindex[2])
       
        if numzero == 0:
            hd=fractions.Fraction(1/h).limit_denominator(100).denominator
            kd=fractions.Fraction(1/k).limit_denominator(100).denominator
            ld=fractions.Fraction(1/l).limit_denominator(100).denominator
            lst=[hd,kd,ld]
            multfact=1
            for ii in lst:
                multfact=multfact*ii / fractions.gcd(multfact,ii)

            p1=np.array([multfact/h,0,0])
            p2=np.array([0,multfact/k,0])
            p3=np.array([0,0,multfact/l])
            
            P=np.array([p1,p2,p3])
            for aa in range(0,3):
                v1=P[np.mod(aa+1,3),:]-P[np.mod(aa,3),:]
                v2=P[np.mod(aa+2,3),:]-P[np.mod(aa,3),:]
                T=self.ang(v1,v2);
                if T<=pi/2:
                    twovects=np.array([v1,v2]);
                    break;
            
        elif numzero==1:
            ind=list()
            P=list()
            for jj in range(0,3):
                if jj==zeroind:
                    p1=np.array([0, 0, 0])
                    p1[jj]=1
                    v1=p1
                else:
                    ind.append(jj)
                    P.append(maxindex[jj])
            
            ad=fractions.Fraction(1/P[0]).limit_denominator(100).denominator
            bd=fractions.Fraction(1/P[1]).limit_denominator(100).denominator
            lst=[ad,bd]
            multfact=1
            for ii in lst:
                multfact=multfact*ii / fractions.gcd(multfact,ii)
                
            points=np.zeros((2,3));
            for mm in range(0,2):
                pointtemp=np.array([0, 0, 0])
                pointtemp[ind[mm]]=multfact/P[mm]
                points[mm,:]=pointtemp
            
            v2=points[1,:]-points[0,:]
            
            twovects=np.array([v1,v2])
            P=np.array([p1, points[0],points[1]])
            
        elif numzero==2:
            maxindex=maxindex/npl.norm(maxindex);
            b = []
            for i in maxindex:
                b.append(abs(float(i)))
            maxindex=b        
            if maxindex == [1, 0, 0]:
                twovects=np.array([[0, 0, 1],[0, 1, 0]])
            elif maxindex==[0, 1, 0]:
                twovects=np.array([[1, 0, 0],[0, 0, 1]])
            elif maxindex==[0, 0, 1]:
                twovects=np.array([[1, 0, 0],[0, 1, 0]])

            P=np.array([twovects[0], twovects[1], [0, 0,0]]);
            
        twovects[0]=twovects[0]/fractions.gcd(fractions.gcd(twovects[0,0], twovects[0,1]),twovects[0,2])
        twovects[1]=twovects[1]/fractions.gcd(fractions.gcd(twovects[1,0], twovects[1,1]),twovects[1,2])
        twovects=np.dot(twovects,basis);
        return twovects,P 
        
    def treflection(self,v):
        n=v/npl.norm(v);
        a=n[0]; b=n[1]; c=n[2];
        R=np.eye(3,3)+np.array([[-2*a**2, -2*a*b, -2*a*c],[-2*a*b, -2*b**2, -2*b*c],[-2*a*c, -2*b*c, -2*c**2]])
        return R

    def trotation(self,ang, axis):
        u=axis/npl.norm(axis)
        ang=np.radians(ang);
        ux=u[0]; uy=u[1]; uz=u[2];
        usubx=np.array([[0, -uz, uy],[uz, 0, -ux],[ -uy, ux, 0]])
        uxu=np.array([[ux**2, ux*uy, ux*uz],[ux*uy, uy**2, uy*uz],[ux*uz, uy*uz, uz**2]])
        Rot=np.eye(3,3)*math.cos(ang)+math.sin(ang)*usubx+(1 - math.cos(ang))*uxu
        return Rot
    
    def getmillerfrom2v(self,basis,v1,v2):
        TM=np.eye(3,3)
        millerv1=np.dot(TM,np.transpose(v1))
        millerv2=np.dot(TM,np.transpose(v2))
        
        millerv1=np.dot(npl.inv(np.transpose(basis)),millerv1)
        millerv2=np.dot(npl.inv(np.transpose(basis)),millerv2)
        
        millerindex=np.transpose(np.cross(millerv1,millerv2))
        
        md=np.zeros((3,1))
        
        for ni in range(0,3):
            md[ni]=fractions.Fraction(millerindex[ni]).limit_denominator(10).denominator
        
        MLCM=1
        for ii in md:
            MLCM=MLCM*ii / fractions.gcd(MLCM,ii)
        
        millerindex=millerindex*MLCM;
        roundedmillerindex=np.array([ round(elem,1) for elem in millerindex])
        diffmill=roundedmillerindex-millerindex
        if npl.norm(diffmill)>0.001:
            raise SystemError('EGREGIOUS ERROR, MILLER INDEX IS FRACTIONAL')
        
        for elem in range(0,3):
            if roundedmillerindex[elem]==-0:
                roundedmillerindex[elem]=0
        
        millerindex=roundedmillerindex
        MGCD=np.abs(fractions.gcd(fractions.gcd(millerindex[0],millerindex[1]),millerindex[2]))
        millerindex=millerindex/MGCD
        return millerindex    
        
    def getsymfam(self,basis, mindex, pg):
        A,useless=self.get2vectsinplane(basis,mindex)
        v1=A[0]
        v2=A[1]
        millerindex=self.getmillerfrom2v(basis, v1, v2)
        
        if millerindex.all()==mindex.all()*-1:
            v2=A[0]
            v1=A[1]
        C=[]
        pgR=len(pg)
        for pgx in range(0,pgR):
            TM=np.eye(3,3)
            T=pg[pgx,0]
            ang=pg[pgx,1]
            axis=np.array(pg[pgx,2:5])
            
            if T==0:
                TM=np.dot(TM,self.treflection(axis))
            elif np.abs(T)==1:
                if T==-1:
                    TM=TM*-1
                else:
                    TM=np.eye(3,3)
            else:
                TM=np.dot(TM,self.trotation(ang,axis))
                if T<0:
                    TM=TM*-1
                
            millerv1=np.dot(TM, np.transpose(v1))   
            millerv2=np.dot(TM, np.transpose(v2))
            millerindex=np.array(self.getmillerfrom2v(basis, millerv1, millerv2))
            exists=False
            for Ccheck in C:
                c=0
                for ind in range(0,3):
                    if Ccheck[ind]==millerindex[ind]:
                        c+=1
                if c==3:
                    exists=True                    
            if exists==False:
                C.append(millerindex)    
        
        for pgx in range(0,pgR):
            T=pg[pgx,0]
            if T==0:
                ang=pg[pgx,1]
                axis=np.array(pg[pgx,2:5])
                for Cx in range(0,len(C)):
                    A,useless=self.get2vectsinplane(basis,C[Cx])
                    v1=A[0]
                    v2=A[1]
                    TM=np.eye(3,3)
                    TM=np.dot(TM,self.treflection(axis))
                    millerv1=np.dot(TM, np.transpose(v1))
                    millerv2=np.dot(TM, np.transpose(v2))
                    millerindex=self.getmillerfrom2v(basis, millerv1, millerv2)
                    exists=False
                    for Ccheck in C:
                        if (Ccheck==millerindex).all():
                            exists=True                    
                    if exists==False:
                        C.append(millerindex)
                        
        retC=np.zeros((len(C),3))
        for Cx in range(0,len(C)):
            retC[Cx]=C[Cx]
        a=retC
        a=a[a[:,2].argsort(),]
        a=a[a[:,1].argsort(),]
        a=a[a[:,0].argsort(),]
        return a    
    
    def proj(self,u,v):
        p=np.dot(v,u)/np.dot(u,u)*u    
        return p
    
    def sproj(self,u,v):
        p=np.dot(v,u)/np.dot(u,u)*npl.norm(u)    
        return p
    
    def gramschmidt(self,B):
        b1=B[0]
        b2=B[1]
        b3=B[2]
        f1=b1
        f2=b2-self.proj(f1,b2)
        f3=b3-self.proj(f1,b3)-self.proj(f2,b3)
        
        f1hat=f1/npl.norm(f1)
        f2hat=f2/npl.norm(f2)
        f3hat=f3/npl.norm(f3)
        
        F=np.array([f1hat,f2hat,f3hat])
        return F

    def getconventional(self,bulkstruct):
        output=run_aconvasp_command(["aconvasp", "--sconv"],bulkstruct)
        if "failed" in output[1]:
            print output
            raise SystemError("Aconvasp failed")
        poscar_string=" "
        for line in output[0].split("\n"):
            poscar_string = poscar_string + line + "\n"
        
        conv=Poscar.from_string(poscar_string).struct
        
        if "HEX" in poscar_string:
            output1=run_aconvasp_command(["aconvasp", "--abccar"], conv)
            hex_string=[]
            basis=[]
            for line in output1[0].split("\n"):
                hex_string.append(line)
            for line in hex_string[2].split():
                basis.append(line)
            a=float(basis[0]);b=float(basis[1]);c=float(basis[2]);aa=float(basis[3]);ab=float(basis[4]);ac=float(basis[5])
            
            if np.abs(a-b)<0.001 and aa==ab and aa-90==0 and ac-120==0:
                tbx=0.866025*b
                tby=-0.5*b
                hexlat = Lattice([[0, a, 0], [tbx, tby, 0], [0, 0, c]])
                sconv=Structure(hexlat,conv.species,conv.frac_coords)
                
            else:
                print "Check Hex Error"
        else:
            sconv=conv
            
        return sconv
                
    def getpg(self,convstruct):
        output=run_aconvasp_command(["aconvasp", "--pgroup"],convstruct)
        pg_string=[]
        for line in output[0].split("\n"):
            if "theta" in line:
                fixline=line.replace(")"," )").replace("(","( ").replace(" -I","I").replace("=","= ").replace("m","0").replace("unity I","unityI")
                splitfixline=fixline.split()
                pgline=np.array([float(splitfixline[11]),float(splitfixline[3]),float(splitfixline[6]),float(splitfixline[7]),float(splitfixline[8])])
                pg_string.append(pgline)
        pg=np.zeros((len(pg_string),5))
        for PAx in range(0,len(pg_string)):
            for PAy in range(0,len(pg_string[PAx])):
                pg[PAx][PAy]=pg_string[PAx][PAy]
        return pg
            
        