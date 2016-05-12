import time
import scipy
import sys
import numpy as np
class SurfaceModule():


        handles = []
        ### pos of center, dir of normal, dir of tangential, dir
        ### of binormal, extension in tangential dir, extension in
        ### binormal dir
        surfaces = []
        def __init__(self, surfaces_in):
                self.surfaces = surfaces_in

        def GetRelevantSurfaces(self, p):

                dbest = float('Inf')
                ibest = -1
                for i in range(0,self.surfaces.shape[0]):
                        di = float('Inf')
                        srfc = self.surfaces[i,:,:]
                        center = srfc[0,:]
                        dnormal = srfc[1,:]
                        dtangential = srfc[2,:]
                        dbinormal = srfc[3,:]
                        ext_t = srfc[4,0]
                        ext_o = srfc[5,0]

                        dp = p - center

                        if np.dot(dp,dnormal)<0:
                                continue

                        di = np.linalg.norm(dp)
                        if di < dbest:
                                dbest = di
                                ibest = i
                print "closest surface:",ibest,dbest
                return ibest

        def SampleSurface(self, M, k, env):
                srfc = self.surfaces[k,:,:]
                center = srfc[0,:]
                dnormal = srfc[1,:]
                dtangential = srfc[2,:]
                dbinormal = srfc[3,:]
                ext_t = srfc[4,0]
                ext_o = srfc[5,0]

                k = 0
                while k < M:
                        rt = np.random.uniform(-ext_t,ext_t)
                        ro = np.random.uniform(-ext_o,ext_o)

                        p = rt*dtangential + ro*dbinormal + center

                        self.handles.append(env.env.plot3(points=p,
                                           pointsize=0.05,
                                           colors=np.array(((1.0,0.0,0.0,0.8))),
                                           drawstyle=1))
                        k = k+1
