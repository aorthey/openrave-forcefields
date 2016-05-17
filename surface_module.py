import time
import scipy
import sys
import numpy as np
from util import Rax
from math import pi
class SurfaceModule():

        OFFSET_CONTACT_FROM_SURFACE_BOUNDARY = 0.05
        #OFFSET_CONTACT_TO_PLANE = 0.008
        OFFSET_HAND_CONTACT_TO_PLANE = 0.02
        OFFSET_FOOT_CONTACT_TO_PLANE = 0.001
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

        def checkNormalization(self, v):
                if not(np.linalg.norm(v)-1 < 1e-5):
                        print v
                        print np.around(np.linalg.norm(v),decimals=10)
                        print "ERROR in normalization"
                        sys.exit(0)

        def GetRotationFromTo( self, v1,v2,v3, w1,w2,w3):
                self.checkNormalization(v1)
                self.checkNormalization(v2)
                self.checkNormalization(v3)
                self.checkNormalization(w1)
                self.checkNormalization(w2)
                self.checkNormalization(w3)

                V1 = np.array((v1,v2,v3))
                W1 = np.array((w1,w2,w3))

                ### solve V1 = R*W1 for R
                R, res, rank, s = np.linalg.lstsq(W1, V1)

                if abs(np.linalg.det(R) - 1 ) > 1e-5:
                        print np.linalg.det(R)
                        print "determinant not 1"
                        sys.exit(1)

                return R

        def GetNearestPointOnSurface(self, dist_to_surface, p, k, env=None):
                ### do not make a contact at the boundary

                srfc = self.surfaces[k,:,:]
                center = srfc[0,:]
                dnormal = srfc[1,:]
                dtangential = srfc[2,:]
                dbinormal = srfc[3,:]
                ext_t = srfc[4,0]-self.OFFSET_CONTACT_FROM_SURFACE_BOUNDARY
                ext_o = srfc[5,0]-self.OFFSET_CONTACT_FROM_SURFACE_BOUNDARY

                ### project p onto plane
                dp = p - center
                dn = np.dot(dp, dnormal)
                dplane = dp - dn*dnormal

                ### add small offset to not be exactly in the plane (otherwise
                ### there will be a collision between robot and env)
                dt = np.dot(dplane,dtangential)
                do = np.dot(dplane,dbinormal)

                if dt > ext_t:
                        dplane = dplane - (dt-ext_t)*dtangential
                if dt < -ext_t:
                        dplane = dplane - (dt+ext_t)*dtangential

                if do > ext_o:
                        dplane = dplane - (do-ext_o)*dbinormal
                if do < -ext_o:
                        dplane = dplane - (do+ext_o)*dbinormal

                #A = env.env.drawarrow(p1=center+dp,p2=center+dplane,linewidth=0.01,color=np.array((1,1,1)))
                #self.handles.append(A)

                ### TODO: make it a cone
                #A = env.env.drawarrow(p1=center+dplane,p2=center+dplane+0.3*dnormal,linewidth=0.01,color=np.array((1,0,1)))
                #self.handles.append(A)

                center = center + dist_to_surface*dnormal
                dplane = dplane + dist_to_surface*dnormal

                return center + dplane

        def GetNearestContactTransformLeftHand(self, env, T, k):
                return self.GetNearestContactTransformHand(env, T, k, -pi/2)

        def GetNearestContactTransformRightHand(self, env, T, k):
                return self.GetNearestContactTransformHand(env, T, k, pi/2)

        def GetNearestContactTransformHand(self, env, T, k, thetaNormal):
                dist_to_surface = self.OFFSET_HAND_CONTACT_TO_PLANE
                T = self.GetNearestContactTransform(dist_to_surface, env, T, k)

                srfc = self.surfaces[k,:,:]
                center = srfc[0,:]
                dnormal = srfc[1,:]
                dtangential = srfc[2,:]
                dbinormal = srfc[3,:]

                R = T[0:3,0:3]
                R = np.dot(Rax(pi/2, dbinormal),R)
                R = np.dot(Rax(thetaNormal, dnormal),R)
                T[0:3,0:3]=R

                return T

        def GetNearestContactTransformFoot(self, env, T, k):
                dist_to_surface = self.OFFSET_FOOT_CONTACT_TO_PLANE
                T = self.GetNearestContactTransform(dist_to_surface, env, T, k)
                return T

        def GetNearestContactTransform(self, dist_to_surface, env, T, k):
                srfc = self.surfaces[k,:,:]
                center = srfc[0,:]
                dnormal = srfc[1,:]
                dtangential = srfc[2,:]
                dbinormal = srfc[3,:]

                ex = np.array((1,0,0))
                ey = np.array((0,1,0))
                ez = np.array((0,0,1))

                c = self.GetNearestPointOnSurface(dist_to_surface, T[0:3,3], k, env)
                R = self.GetRotationFromTo( ex, ey, ez, dtangential, dbinormal, dnormal)

                T[0:3,3] = c
                T[0:3,0:3] = R

                DEBUG = False
                if DEBUG:
                        V0 = np.dot(T,np.array((0,0,0,1)))[0:3]
                        size = 0.5
                        V1 = np.dot(T,np.array((size,0,0,1)))[0:3]
                        V2 = np.dot(T,np.array((0,size,0,1)))[0:3]
                        V3 = np.dot(T,np.array((0,0,size,1)))[0:3]
                        
                        lw = 0.03
                        A = env.env.drawarrow(p1=V0,p2=V1,linewidth=lw,color=np.array((1,0,0)))
                        self.handles.append(A)
                        A = env.env.drawarrow(p1=V0,p2=V2,linewidth=lw,color=np.array((0,1,0)))
                        self.handles.append(A)
                        A = env.env.drawarrow(p1=V0,p2=V3,linewidth=lw,color=np.array((0,0,1)))
                        self.handles.append(A)
                        lw = 0.02
                        A = env.env.drawarrow(p1=V0,p2=V0+dtangential,linewidth=lw,color=np.array((1,0.3,0.3)))
                        self.handles.append(A)
                        A = env.env.drawarrow(p1=V0,p2=V0+dbinormal,linewidth=lw,color=np.array((0.3,1,0.3)))
                        self.handles.append(A)
                        A = env.env.drawarrow(p1=V0,p2=V0+dnormal,linewidth=lw,color=np.array((0.3,0.3,1)))
                        self.handles.append(A)


                return T

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


