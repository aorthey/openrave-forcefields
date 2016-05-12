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

        def GetNearestPointOnSurface(self, p, k):
                srfc = self.surfaces[k,:,:]
                center = srfc[0,:]
                dnormal = srfc[1,:]
                dtangential = srfc[2,:]
                dbinormal = srfc[3,:]
                ext_t = srfc[4,0]
                ext_o = srfc[5,0]

                ### project p onto plane
                dp = p - center
                dn = np.dot(dp, dnormal)
                dplane = dp - dn*dnormal
                print dn,dplane

                dt = np.dot(dplane,dtangential)
                do = np.dot(dplane,dbinormal)

                if abs(dt) > ext_t:
                        dplane = dplane - (dt-ext_t)*dtangential
                if abs(do) > ext_o:
                        dplane = dplane - (do-ext_o)*dbinormal

                return center + dplane

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
                        pn = p + 0.2*dnormal

                        A = env.env.drawarrow(p1=p,p2=pn,linewidth=0.005,color=np.array((1,0.5,0)))
                        self.handles.append(A)
                        #self.handles.append(env.env.plot3(points=p,
                        #                   pointsize=0.05,
                        #                   colors=np.array(((1.0,0.0,0.0,0.8))),
                        #                   drawstyle=1))
                        k = k+1

        def GetNearestContactTransform(self, env, T, k):
                srfc = self.surfaces[k,:,:]
                center = srfc[0,:]
                dnormal = srfc[1,:]
                dtangential = srfc[2,:]
                dbinormal = srfc[3,:]

                ex = np.array((1,0,0))
                ey = np.array((0,1,0))
                ez = np.array((0,0,1))

                c = self.GetNearestPointOnSurface(T[0:3,3], k)
                R = self.GetRotationFromTo( ex, ey, ez, dtangential, dbinormal, dnormal)

                T[0:3,0:3] = R
                T[0:3,3] = c

                DEBUG = 1
                if DEBUG:
                        V0 = np.dot(T,np.array((0,0,0,1)))[0:3]
                        V1 = np.dot(T,np.array((1,0,0,1)))[0:3]
                        V2 = np.dot(T,np.array((0,1,0,1)))[0:3]
                        V3 = np.dot(T,np.array((0,0,1,1)))[0:3]
                        
                        A = env.env.drawarrow(p1=V0,p2=V1,linewidth=0.03,color=np.array((1,0,1)))
                        self.handles.append(A)
                        A = env.env.drawarrow(p1=V0,p2=V2,linewidth=0.03,color=np.array((1,0,1)))
                        self.handles.append(A)
                        A = env.env.drawarrow(p1=V0,p2=V3,linewidth=0.03,color=np.array((1,0,1)))
                        self.handles.append(A)


                return T


        def GetRotationFromTo( self, v1,v2,v3, w1,w2,w3):
                B = np.outer(w1,v1.T)+np.outer(w2,v2.T)+np.outer(w3,v3.T)
                #B = np.outer(v3,w3.T)
                U,S,V = np.linalg.svd(B)
                du = np.linalg.det(U)
                dv = np.linalg.det(V)
                M = np.diag((1,1,du*dv))
                R = np.dot(np.dot(U,M),V.T)
                return R



