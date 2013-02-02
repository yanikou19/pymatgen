'''
Created on Jan 31, 2013

@author: wenhao
'''

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
import numpy.linalg as npl
from pyhull.convex_hull import ConvexHull, qconvex

class wulff(object):
    def __init__(self, convstruct, indexSEdict):
        self.convstructure=convstruct
        self.indexSEdict=indexSEdict
        basis=np.array(convstruct.lattice.matrix)
        pg=st.Surftool().getpg(convstruct)
        
        '''Generate all equivalent surfaces from the symmetry of the lattice'''
        famsurfdict=dict()
        allsurfaces=[]
        for key in indexSEdict.keys():
            index=st.Surftool().indexfromstr(key)
            indexfam=st.Surftool().getsymfam(basis, index, pg)
            indexdict=dict()
            for ind in range(0,len(indexfam)):
                I=indexfam[ind]
                strI="["+str(int(I[0]))+","+str(int(I[1]))+","+str(int(I[2]))+"]"
                E=float(indexSEdict[key])
                indexdict[strI]=E
                allsurfappend=[np.array([int(I[0]),int(I[1]),int(I[2])]),E,index]
                allsurfaces.append(allsurfappend)
            famsurfdict[key]=indexdict
            
        '''Make Wulff Shape'''
            
        '''First, Make Gamma(n)'''
        gn=[]
        dualgn=[]
        
        for surf in allsurfaces:
            millerindex=surf[0]
            gamma=surf[1]
            twovects,P=st.Surftool().get2vectsinplane(basis, millerindex)
            v1=twovects[0]
            v2=twovects[1]
            mindex=st.Surftool().getmillerfrom2v(basis, v1, v2)
            thedot=np.dot(mindex,millerindex)/npl.norm(mindex)/npl.norm(millerindex)
            if 0.99<thedot and thedot < 1.01:
                normal=np.cross(v1,v2)
            elif -1.01<thedot and thedot < -0.99:
                normal=np.cross(v2,v1)
            else:
                raise SystemError("Gamma(n) error")
            
            normal=gamma*normal/(npl.norm(normal))
            gn.append(normal)
            
            dualnormal=normal/(npl.norm(normal))/gamma
            dualgn.append(dualnormal)
            
        K=qconvex("i",dualgn)
        K.pop(0)
        
        offhull=[]
        onhull=[]
        uniquek=str(K).replace("[","").replace("'","").replace(",","").replace("]","").split()
        uniquek=np.unique(np.array(uniquek))
        
        for sdgn in range(0,len(dualgn)):
            contains=1
            for suk in range(0,len(uniquek)):
                uk=int(uniquek[suk])
                if uk==sdgn:
                    contains=0
            if contains==1:
                offhullapp=[allsurfaces[sdgn][2],sdgn,dualgn[sdgn],[]]
                offhull.append(offhullapp)
            else:
                onhullapp=[allsurfaces[sdgn],sdgn,dualgn[sdgn]]
                onhull.append(onhullapp)
        
        surfaceabovehull=[]
        oldinds=[]
        if len(offhull)>0:
            finduniques=[]
            for off in offhull:
                finduniques.append(off[0])
            finduniques=np.array(finduniques)
            uniquesurfaces=np.array([np.array(x) for x in set(tuple(x) for x in finduniques)]) 
            for asius in range(0,len(uniquesurfaces)):
                for off in offhull:
                    if not any((off[0] == x).all() for x in oldinds):
                        if uniquesurfaces[asius].all()==off[0].all():
                            I=uniquesurfaces[asius]
                            strI="["+str(int(I[0]))+","+str(int(I[1]))+","+str(int(I[2]))+"]"
                            surfaceabovehull.append(off)
                            oldinds.append(off[0])
                      
        '''Determine dual space convex hull points, also calculate energy above the hull'''
        
        abcd=[]
        wulffverts=[]
        dabovehull=np.ones((len(surfaceabovehull),1))*1000
        ptonhull=np.zeros((len(surfaceabovehull),3))
        for ind in range(0,len(K)):
            Kind=K[ind].split()
            p1=dualgn[int(Kind[0])]
            p2=dualgn[int(Kind[1])]
            p3=dualgn[int(Kind[2])]
            for jj in range(0,len(surfaceabovehull)):
                linevector=np.array(surfaceabovehull[jj][2])
                PT=self.intersectLinePlane(p1,p2,p3,linevector)
                D=self.distancePoints3d(PT,linevector)
                if D<=dabovehull[jj]:
                    dabovehull[jj]=D
                    ptonhull[jj]=PT
            chv1=p3-p1
            chv2=p2-p1
            N=np.cross(chv1,chv2)
            if np.dot(N,p1)-np.dot(N,p2)<1E-3 and np.dot(N,p1)-np.dot(N,p3)<1E-3:
                d=np.dot(N,p1)
                if np.sign(d)==-1:
                    N=np.cross(chv2,chv1)
                    d=np.dot(N,p1)
            else:
                raise SystemError("plane error")
        
            dist=d/npl.norm(N)
            abcdappend=[N,d,dist]
            abcd.append(abcdappend)
            RN=N/npl.norm(N)/dist
            wulffverts.append(RN)
        
        
        for ii in range(0,len(surfaceabovehull)):
            surfaceabovehull[ii][3]=dabovehull[ii]
        
        eabvwulff=[]
        for eabvsurf in range(0,len(surfaceabovehull)):
            N=surfaceabovehull[eabvsurf][2]
            PT=ptonhull[eabvsurf]
            normal=N/npl.norm(N)
            if np.abs((np.dot(normal,N)/npl.norm(normal)/npl.norm(N))-1)<1E-6:
                gamma=npl.norm(normal)/npl.norm(N)
                hullgamma=npl.norm(normal)/npl.norm(PT)
                eabvhull=gamma-hullgamma
                eabvwulffappend=[surfaceabovehull[eabvsurf][0],eabvhull]
                
                eabvwulff.append(eabvwulffappend)
            else:
                raise SystemError("Eabvwulff error")
        self.eabvwulff=eabvwulff
        self.wulff=qconvex("FA",wulffverts)
        
    def getEAbvWulff(self):
        Eabvwulffdict=dict()
        for e in self.eabvwulff:
            I=e[0]
            strI="["+str(int(I[0]))+","+str(int(I[1]))+","+str(int(I[2]))+"]"
            Eabvwulffdict[strI]=e[1]
        return Eabvwulffdict
    
    def intersectLinePlane(self,p1,p2,p3,linevector):
        #lb=[0,0,0]
        #p1, p2 and p3 are three points on the plane
        xa=linevector[0];ya=linevector[1];za=linevector[2]
        x1=p1[0];y1=p1[1];z1=p1[2] 
        x2=p2[0];y2=p2[1];z2=p2[2]
        x3=p3[0];y3=p3[1];z3=p3[2]
        
        mata=np.array([[xa,x2-x1,x3-x1],
                       [ya,y2-y1,y3-y1],
                       [za,z2-z1,z3-z1]])
        
        tuv=np.dot(npl.inv(mata),np.array([[xa-x1],[ya-y1],[za-z1]]))
        t=tuv[0]
        return linevector*(1-t)
        
        
    def distancePoints3d(self,p1,p2):
        return np.sqrt((p2[0]-p1[0])**2+(p2[1]-p1[1])**2+(p2[2]-p1[2])**2)