import abc
import sys
from openravepy import *
from numpy import array,pi
from math import cos,sin
from util import Rz, Rax
import math
import numpy as np

class ForceEnvironment():
        __metaclass__ = abc.ABCMeta
        #######################################################################
        ZPOS_ARROW = 0.1
        ViewerName = 'qtcoin'
        ##'bullet', 'ode'
        PhysicsEngineName = 'bullet'
        ##'bullet', 'ode', 'pqp', 'fcl_'
        CollisionCheckerName = 'ode'
        FORCE_FIELD_MIN_SPACING = 0.3
        FORCE_FIELD_MAX_SPACING = 0.8
        FORCE_FIELD_PT_SIZE = 4
        FORCE_FIELD_ARROW_SIZE = 0.02
        #FORCE_FIELD_COLOR = np.array((0.0,0.0,0.0))
        FORCE_FIELD_COLOR = np.array((1.0,1.0,1.0,0.4))
        COLOR_CONTACT_PATCH = np.array((0,1,0,1))
        #######################################################################
        env_xml=''
        robot_xml=''
        env=''
        robot=''
        handles = []
        force_handles = []
        extra_handles = []
        static_handles = []
        cells = None
        forces = None

        def ResetForceHandles(self):
                force_handles = []

        def AddForceHandles(self,H):
                if H is None:
                        return
                else:
                        for h in H:
                                self.force_handles.append(h)

        def __init__(self):
                self.env=Environment()

                self.env.SetViewer(self.ViewerName)
                self.env.Reset()
                with self.env:
                        self.env.StopSimulation() 
                        self.physics = RaveCreatePhysicsEngine(self.env, self.PhysicsEngineName)
                        self.physics.SetGravity(array((0,0,-9.80)))
                        self.env.SetPhysicsEngine(self.physics)
                        self.cc = RaveCreateCollisionChecker(self.env, self.CollisionCheckerName)
                        self.env.SetCollisionChecker(self.cc)
                        #self.env.SetForces( self.GetForces() )
                        #self.recorder = RaveCreateModule(self.env,'viewerrecorder')
                        #self.env.AddModule(self.recorder,'')

        def setrobotenv(self,robot_xml,env_xml):
                with self.env:
                        self.env_xml = env_xml
                        self.robot_xml = robot_xml
                        self.env.Add(self.env.ReadRobotXMLFile(self.robot_xml))
                        self.env.Load(self.env_xml)

        def GetCellsAll(self, verbose=True):
                with self.env:
                        B = self.env.GetBodies()[1]
                        if verbose:
                                print B.GetLinks()
                        return B.GetLinks()

        @abc.abstractmethod
        def GetCells(self):
                pass

        @abc.abstractmethod
        def GetForces(self):
                pass

        @abc.abstractmethod
        def RobotGetInitialPosition(self):
                pass

        @abc.abstractmethod
        def RobotGetGoalPosition(self):
                pass
        def MakeRobotTransparent(self, alpha):
                assert( 0 <= alpha <= 1)
                with self.env:
                        self.robot = self.env.GetRobots()[0]
                        for link in self.robot.GetLinks():
                                for geom in link.GetGeometries():
                                        geom.SetTransparency(1.0-alpha)

        def MakeRobotInvisible(self):
                self.MakeRobotTransparent(0.0)
        def MakeRobotVisible(self):
                self.MakeRobotTransparent(1.0)

        def AddRandomDisturbance(self):
                maxforce = 10
                r = (numpy.random.rand(2)-0.5)
                F = np.array((r[0]*maxforce,r[1]*maxforce,0))
                self.AddForceToWorld(F)

        def AddForceToWorld(self, F):
                zero = np.array((0,0,0))
                W = self.env.GetKinBody('world')
                for link in W.GetLinks():
                        P = link.GetLocalCOM()
                        link.SetForce(F,P,False)
                        link.SetTorque(zero, False)

        def GetRobot(self):
                with self.env:
                        self.robot = self.env.GetRobots()[0]

                        [xi,yi,zi,ti,dxi,dyi,dzi,dti]=self.RobotGetInitialPosition()
                        [xg,yg,zg,tg,dxg,dyg,dzg,dtg]=self.RobotGetGoalPosition()
                        #self.robot.SetDOFLimits((-10,10),(-5,5),(-1,1))
                        #self.robot.SetDOFLimits((-10,-8,-0.15,-4*pi),(10,8,0.15,4*pi))
                        self.robot.SetDOFLimits((-10,-10,-0.2,-2*pi),(10,10,0.2,2*pi))
                        self.robot.SetDOFValues((xi,yi,zi,ti))
                        self.robot.SetDOFVelocities((dxi,dyi,dzi,dti))
                        self.robot.SetDOFVelocityLimits([10.0,10.0,0.0,5.0])
                        self.robot.SetDOFAccelerationLimits([5.0,5.0,0.0,5.0])


                        self.handles.append(self.env.plot3(points=array(((xg,yg,zg))),
                                           pointsize=0.08,
                                           colors=array(((1.0,0.0,0.0,0.8))),
                                           drawstyle=1))

                        self.handles.append(self.env.plot3(points=array(((xi,yi,zg))),
                                           pointsize=0.08,
                                           colors=array(((0.0,1.0,0.0,0.8))),
                                           drawstyle=1))

                        d = 0.3
                        self.handles.append(self.env.drawlinestrip(points=np.array(((xi,yi,zi),(xi+d*cos(ti),yi+d*sin(ti),zi))),
                                           linewidth=5,
                                           colors=np.array(((0.0,0.0,0.0,0.8)))))
                        self.handles.append(self.env.drawlinestrip(points=np.array(((xg,yg,zg),(xg+d*cos(tg),yg+d*sin(tg),zg))),
                                           linewidth=5,
                                           colors=np.array(((0.0,0.0,0.0,0.8)))))
                        return self.robot

        def DrawBorderAroundCell(self, cell):
                handle = []
                with self.env:
                        G1 = cell.GetGeometries()[0]
                        B = G1.GetBoxExtents()
                        T = G1.GetTransform()
                        h = self.DrawBoxMesh(T,B)
                        handle.append(h)
                        return handle

        ### transform matrix T \in [3x3]
        ### box extend vector B \in [3x1]
        def DrawBoxMesh(self, T, B):
                        ########################################################
                        ## visualize extend of force constraint box
                        ########################################################
                        P = array(((T[0,3]-B[0],T[1,3]-B[1],B[2]+T[2,3]+0.01), \
                                (T[0,3]+B[0],T[1,3]-B[1],B[2]+T[2,3]+0.01), \
                                (T[0,3]+B[0],T[1,3]+B[1],B[2]+T[2,3]+0.01), \
                                (T[0,3]-B[0],T[1,3]+B[1],B[2]+T[2,3]+0.01),\
                                (T[0,3]-B[0],T[1,3]-B[1],B[2]+T[2,3]+0.01),\
                                ## go up
                                (T[0,3]-B[0],T[1,3]-B[1],T[2,3]-B[2]+0.01), \
                                ## next point
                                (T[0,3]+B[0],T[1,3]-B[1],T[2,3]-B[2]+0.01), \
                                (T[0,3]+B[0],T[1,3]-B[1],T[2,3]+B[2]+0.01), \
                                (T[0,3]+B[0],T[1,3]-B[1],T[2,3]-B[2]+0.01), \
                                ## next point
                                (T[0,3]+B[0],T[1,3]+B[1],T[2,3]-B[2]+0.01), \
                                (T[0,3]+B[0],T[1,3]+B[1],T[2,3]+B[2]+0.01), \
                                (T[0,3]+B[0],T[1,3]+B[1],T[2,3]-B[2]+0.01), \
                                ## next point
                                (T[0,3]-B[0],T[1,3]+B[1],T[2,3]-B[2]+0.01),\
                                (T[0,3]-B[0],T[1,3]+B[1],T[2,3]+B[2]+0.01),\
                                (T[0,3]-B[0],T[1,3]+B[1],T[2,3]-B[2]+0.01),\
                                ## back start
                                (T[0,3]-B[0],T[1,3]-B[1],T[2,3]-B[2]+0.01)))

                        h=self.env.drawlinestrip(points=P,linewidth=5.0,colors=array(((0.5,0.5,0.5,0.5))))
                        return h

        def GetZeroForce(self):
                return np.zeros(self.forces[0].shape)

        collision_handler = []

        def CheckCollisionAtX(self, X):
                self.robot.SetDOFValues(X)

                floorDim = len(self.cells)
                outercells = self.GetCellsAll(verbose=False)[floorDim:]

                #self.collision_handler = []
                for i in range(0,len(outercells)):
                        robotIsInCollision = self.env.CheckCollision(self.robot, outercells[i])
                        if robotIsInCollision:
                                print "COLLISION:",X[0:3],outercells[i]
                                A = self.DrawCollision(X[0:3])
                                self.collision_handler.append(A)
                                return True
                return False

        def DrawCollision(self, X):
                return self.env.plot3(points=X,
                                   pointsize=0.05,
                                   colors=array(((1.0,0.0,0.0,1.0))),
                                   drawstyle=1)


        def DrawLeftHandContactPatchFromTransform(self, T, cpatch=None):
                self.DrawHandContactPatchFromTransform(T, -pi/2, cpatch)
        def DrawRightHandContactPatchFromTransform(self, T, cpatch=None):
                self.DrawHandContactPatchFromTransform(T, pi/2, cpatch)

        def DrawHandContactPatchFromTransform(self, T, thetaNormal, cpatch=None):

                Thand = np.eye(4)
                Thand[0:3,3] = T[0:3,3]

                e0 = np.array((0,0,0,1))
                ex = np.array((1,0,0,1))
                ey = np.array((0,1,0,1))
                ez = np.array((0,0,1,1))

                fpos = np.dot(T,e0)[0:3]
                fdir1 = np.dot(T,ex)[0:3] - fpos
                fdir2 = np.dot(T,ey)[0:3] - fpos
                fdir3 = np.dot(T,ez)[0:3] - fpos

                R = T[0:3,0:3]
                R = np.dot(Rax(pi/2, fdir3),R)
                R = np.dot(Rax(thetaNormal, fdir2),R)
                Thand[0:3,0:3]=R

                self.DrawFootContactPatchFromTransform(Thand, cpatch)


        def DrawFootContactPatchFromTransform(self, T, cpatch=None):
                ### geometry of robot feet (escher)
                FOOT_WIDTH = 0.07
                FOOT_LENGTH = 0.12
                lw = 8

                e0 = np.array((0,0,0,1))
                ex = np.array((1,0,0,1))
                ey = np.array((0,1,0,1))
                ez = np.array((0,0,1,1))

                fpos = np.dot(T,e0)[0:3]
                fdir1 = np.dot(T,ex)[0:3] - fpos
                fdir2 = np.dot(T,ey)[0:3] - fpos
                fdir3 = np.dot(T,ez)[0:3] - fpos

                #offset = 0.001*lw
                #foffset = fpos + offset*fdir3
                foffset = fpos

                X=np.zeros((3,7))
                X[:,0] = foffset + FOOT_LENGTH*fdir1 - FOOT_WIDTH*fdir2
                X[:,1] = foffset + FOOT_LENGTH*fdir1 + FOOT_WIDTH*fdir2
                X[:,2] = foffset - FOOT_LENGTH*fdir1 + FOOT_WIDTH*fdir2
                X[:,3] = foffset - FOOT_LENGTH*fdir1 - FOOT_WIDTH*fdir2
                X[:,4] = foffset + FOOT_LENGTH*fdir1 - FOOT_WIDTH*fdir2
                X[:,5] = foffset + FOOT_LENGTH*fdir1
                X[:,6] = foffset

                if cpatch is None:
                        cpatch = self.COLOR_CONTACT_PATCH

                h = self.env.drawlinestrip(points=X.T,
                                linewidth = lw,
                                colors=cpatch)
                self.static_handles.append(h)

        def DrawNormalVectors(self, posN, dirN):
                N = posN.shape[0]
                naxis = np.array((0,0,1))
                lw = 0.01
                for i in range(0,N):
                        A = self.env.drawarrow(p1=posN[i],p2=posN[i]+0.2*dirN[i],linewidth=lw,color=naxis)
                        self.static_handles.append(A)

        def DrawTangentialVectors(self,posN,tN,oN):
                N = posN.shape[0]
                ctaxis = np.array((0,1,0))
                coaxis = np.array((1,0,0))
                lw = 0.01
                for i in range(0,N):
                        A = self.env.drawarrow(p1=posN[i],p2=posN[i]+0.2*tN[i],linewidth=lw,color=ctaxis)
                        self.static_handles.append(A)
                        A = self.env.drawarrow(p1=posN[i],p2=posN[i]+0.2*oN[i],linewidth=lw,color=coaxis)
                        self.static_handles.append(A)

        def DrawSingleCoordinateFrame(self, Pf, Nf, Tf, Bf):
                lw = 0.01
                size = 0.2
                alpha = 0.5
                offset = 2*lw
                nframes = Pf.shape[0]
                cnaxis = np.array((0,0,1,alpha))
                coaxis = np.array((0,1,0,alpha))
                ctaxis = np.array((1,0,0,alpha))
                center = Pf+offset*Nf
                A = self.env.drawarrow(p1=center,p2=center+size*Nf,linewidth=lw,color=cnaxis)
                self.static_handles.append(A)
                A = self.env.drawarrow(p1=center,p2=center+size*Tf,linewidth=lw,color=ctaxis)
                self.static_handles.append(A)
                A = self.env.drawarrow(p1=center,p2=center+size*Bf,linewidth=lw,color=coaxis)
                self.static_handles.append(A)

        def DrawCoordinateFrame(self, Pf, Nf, Tf, Bf):
                nframes = Pf.shape[0]
                for i in range(0,nframes):
                        self.DrawSingleCoordinateFrame(Pf[i], Nf[i], Tf[i], Bf[i])

        def DrawFrameFromTransform(self, T, scale=0.1):
                e0 = np.dot(T,np.array((0,0,0,1)))[0:3]
                e1 = np.dot(T,np.array((scale,0,0,1)))[0:3]
                e2 = np.dot(T,np.array((0,scale,0,1)))[0:3]
                e3 = np.dot(T,np.array((0,0,scale,1)))[0:3]
                caxis1 = np.array((1,0,0))
                caxis2 = np.array((0,1,0))
                caxis3 = np.array((0,0,1))
                lw = 0.01
                A = self.env.drawarrow(p1=e0,p2=e1,linewidth=lw,color=caxis1)
                self.static_handles.append(A)
                A = self.env.drawarrow(p1=e0,p2=e2,linewidth=lw,color=caxis2)
                self.static_handles.append(A)
                A = self.env.drawarrow(p1=e0,p2=e3,linewidth=lw,color=caxis3)
                self.static_handles.append(A)

        def DrawAxes(self):
                e0 = np.array((0,0,0))
                e1 = np.array((1,0,0))
                e2 = np.array((0,1,0))
                e3 = np.array((0,0,1))
                caxis1 = np.array((1,0,0))
                caxis2 = np.array((0,1,0))
                caxis3 = np.array((0,0,1))
                lw = 0.02
                A = self.env.drawarrow(p1=e0,p2=e1,linewidth=lw,color=caxis1)
                self.static_handles.append(A)
                A = self.env.drawarrow(p1=e0,p2=e2,linewidth=lw,color=caxis2)
                self.static_handles.append(A)
                A = self.env.drawarrow(p1=e0,p2=e3,linewidth=lw,color=caxis3)
                self.static_handles.append(A)


        def DrawArrow(self, pos, direction, deleteOld=False, lw = None, ca = None):
                if deleteOld:
                        self.extra_handles = []
                if np.linalg.norm(direction)<0.05:
                        return
                if lw is None:
                        lw = self.FORCE_FIELD_ARROW_SIZE
                if ca is None:
                        ca = self.FORCE_FIELD_COLOR
                A = self.env.drawarrow(p1=pos,p2=pos+direction,linewidth=lw,color=ca)
                self.extra_handles.append(A)

        def GetForceAtX(self, X):
                self.robot.SetDOFValues(X)

                F = self.GetZeroForce()
                for i in range(0,len(self.cells)):
                        robotIsInsideCell = self.env.CheckCollision(self.robot, self.cells[i])
                        if robotIsInsideCell:
                                F = F + self.forces[i]

                return F

        def DrawForceArrowsInCell(self, cell, force):
                G1 = cell.GetGeometries()[0]
                B = G1.GetBoxExtents()
                T = G1.GetTransform()
                return self.DrawForceArrowsInBox(T,B,force)

        def DrawForceArrowsInBox(self, T, B, F):
                ## force in box
                Fx = F[0]
                Fy = F[1]
                Fz = F[2]

                ## middle point of box
                mean_x = T[0,3]
                mean_y = T[1,3]
                mean_z = T[2,3]
                ## extend of box
                ext_x = B[0]
                ext_y = B[1]
                ext_z = B[2]

                if np.sqrt(Fx*Fx+Fy*Fy+Fz*Fz) < 0.001:
                        ## no force to display
                        return None

                lx = np.sign(Fx)*min((abs(Fx)/10.0),10.0)*self.FORCE_FIELD_MAX_SPACING
                ly = np.sign(Fy)*min((abs(Fy)/10.0),10.0)*self.FORCE_FIELD_MAX_SPACING
                lz = np.sign(Fz)*min((abs(Fz)/10.0),10.0)*self.FORCE_FIELD_MAX_SPACING

                ########################################################
                ## compute arrows inside of box (assume rectangular cover)
                ########################################################
                dxspacing = min(max(self.FORCE_FIELD_MIN_SPACING,abs(lx)),self.FORCE_FIELD_MAX_SPACING)
                dyspacing = min(max(self.FORCE_FIELD_MIN_SPACING,abs(ly)),self.FORCE_FIELD_MAX_SPACING)
                dzspacing = min(max(self.FORCE_FIELD_MIN_SPACING,abs(lz)),self.FORCE_FIELD_MAX_SPACING)

                Nx = int(np.floor(2*ext_x/dxspacing))
                Ny = int(np.floor(2*ext_y/dyspacing))
                Nz = int(np.floor(2*ext_z/dzspacing))

                Nx = max(Nx,1)
                Ny = max(Ny,1)
                Nz = max(Nz,1)

                xstart = (mean_x - ext_x) + 0.5*dxspacing
                ystart = (mean_y - ext_y) + 0.5*dyspacing
                zstart = (mean_z - ext_z) + 0.5*dzspacing + self.ZPOS_ARROW

                l = np.sqrt(lx*lx+ly*ly+lz*lz)

                handles=[]
                lx = 0.8*lx
                ly = 0.8*ly
                lz = 0.8*lz
                for i in range(0,Nx):
                        for j in range(0,Ny):
                                for k in range(0,Nz):
                                        x = xstart+i*dxspacing
                                        y = ystart+j*dyspacing
                                        z = zstart+k*dzspacing

                                        scale = 0.9
                                        dx = scale*i*dxspacing
                                        dy = scale*j*dyspacing
                                        dz = scale*k*dzspacing

                                        A = self.env.drawarrow(p1=array((x,y,z)),p2=array((x+lx,y+ly,z+lz)),linewidth=self.FORCE_FIELD_ARROW_SIZE,color=self.FORCE_FIELD_COLOR)
                                        handles.append(A)
                return handles

        recorder = None

        def VideoRecordStart(self, fname):
                with self.env:
                        codecs = self.recorder.SendCommand('GetCodecs') # linux only
                        codec = 13 # mpeg4
                        recorder.SendCommand('Start 640 480 30 codec %d timing realtime filename %s\nviewer %s'%(codec,fname,self.env.GetViewer().GetName()))

        def VideoRecordStop(self):
                        self.recorder.SendCommand('Stop') # stop the video
                        self.env.Remove(self.recorder) # remove the recorder


        def DisplayForces(self):
                if self.forces is None:
                        self.forces = self.GetForces()
                if self.cells is None:
                        self.cells = self.GetCells()

                #with self.env:
                assert(len(self.forces)==len(self.cells))

                self.ResetForceHandles()
                for i in range(0,len(self.cells)):
                        h = self.DrawBorderAroundCell(self.cells[i])
                        self.AddForceHandles(h)
                        h = self.DrawForceArrowsInCell(self.cells[i], self.forces[i])
                        self.AddForceHandles(h)

if __name__ == "__main__":
        env = ForceEnvironment()
        xmlenv='environments/the_stream.env.xml'
        xmlrobot='robots/pointrobot.robot.xml'
        env.setrobotenv(xmlrobot,xmlenv)
