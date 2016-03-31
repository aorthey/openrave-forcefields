import abc
from openravepy import *
from numpy import array,pi
import numpy as np

class ForceEnvironment():
        __metaclass__ = abc.ABCMeta
        #######################################################################
        ZPOS_ARROW = 0.05
        ViewerName = 'qtcoin'
        PhysicsEngineName = 'ode'
        CollisionCheckerName = 'ode'
        #######################################################################
        env_xml=''
        robot_xml=''
        env=''
        robot=''
        handles = []
        force_handles = []
        cells = None
        forces = None

        def ResetForceHandles(self):
                force_handles = []

        def PrintForceHandleInfo(self):
                print len(self.force_handles)
                print self.force_handles

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
                        self.physics = RaveCreatePhysicsEngine(self.env, self.PhysicsEngineName)
                        self.physics.SetGravity(array((0,0,-9.81)))
                        self.env.SetPhysicsEngine(self.physics)
                        self.cc = RaveCreateCollisionChecker(self.env, self.CollisionCheckerName)
                        self.env.SetCollisionChecker(self.cc)
                        self.env.SetForces( self.GetForces() )
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

        def GetRobot(self):
                with self.env:
                        self.robot = self.env.GetRobots()[0]

                        [xi,yi,zi,ti,dxi,dyi,dzi,dti]=self.RobotGetInitialPosition()
                        [xg,yg,zg,tg,dxg,dyg,dzg,dtg]=self.RobotGetGoalPosition()
                        #self.robot.SetDOFLimits((-10,10),(-5,5),(-1,1))
                        #self.robot.SetDOFLimits((-10,-8,-0.15,-4*pi),(10,8,0.15,4*pi))
                        self.robot.SetDOFLimits((-10,-10,-0.15,-2*pi),(10,10,0.15,2*pi))
                        self.robot.SetDOFValues((xi,yi,zi,ti))
                        self.robot.SetDOFVelocities((dxi,dyi,dzi,dti))
                        self.robot.SetDOFVelocityLimits([10.0,10.0,0.0,5.0])
                        self.robot.SetDOFAccelerationLimits([5.0,5.0,0.0,5.0])


                        self.handles.append(self.env.plot3(points=array(((xg,yg,zg))),
                                           pointsize=0.15,
                                           colors=array(((1.0,0.0,0.0,0.8))),
                                           drawstyle=1))

                        self.handles.append(self.env.plot3(points=array(((xi,yi,zg))),
                                           pointsize=0.15,
                                           colors=array(((0.0,1.0,0.0,0.8))),
                                           drawstyle=1))
                        return self.robot

        def DrawBorderAroundCell(self, cell):
                handle = []
                with self.env:
                        G1 = cell.GetGeometries()[0]
                        B = G1.GetBoxExtents()
                        T = G1.GetTransform()

                        ########################################################
                        ## visualize extend of force constraint box
                        ########################################################
                        P = array(((T[0,3]-B[0],T[1,3]-B[1],B[2]+T[2,3]), \
                                (T[0,3]+B[0],T[1,3]-B[1],B[2]+T[2,3]), \
                                (T[0,3]+B[0],T[1,3]+B[1],B[2]+T[2,3]), \
                                (T[0,3]-B[0],T[1,3]+B[1],B[2]+T[2,3]),\
                                (T[0,3]-B[0],T[1,3]-B[1],B[2]+T[2,3])))

                        h=self.env.drawlinestrip(points=P,linewidth=5.0,colors=array(((1,1,1,0.5))))
                        handle.append(h)
                        return handle

        def GetZeroForce(self):
                return np.zeros(self.forces[0].shape)

        def CheckCollisionAtX(self, X):
                self.robot.SetDOFValues(X)

                floorDim = len(self.cells)
                outercells = self.GetCellsAll(verbose=False)[floorDim:]

                for i in range(0,len(outercells)):
                        robotIsInCollision = self.env.CheckCollision(self.robot, outercells[i])
                        if robotIsInCollision:
                                print "COLLISION:",X[0:2],outercells[i]
                                return True
                return False

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
                Fx = force[0]
                Fy = force[1]

                ########################################################
                ## compute arrows inside of box
                ########################################################
                bx = T[0,3]
                by = T[1,3]
                N = int(np.ceil(B[0]))+1
                M = int(np.ceil(B[1]))+1
                dxspacing = B[0]/N
                dyspacing = B[1]/M

                if np.sqrt(Fx*Fx+Fy*Fy) < 0.001:
                        return None

                ### lx,ly == direction vector of line
                maxlx = dxspacing-dxspacing/4
                minlx = -maxlx
                maxly = dyspacing-dyspacing/4
                minly = -maxly
                lx = (Fx/10.0)*maxlx
                ly = (Fy/10.0)*maxly

                if lx > maxlx:
                        lx = maxlx
                if lx < minlx:
                        lx = minlx
                if ly > maxly:
                        ly = maxly
                if ly < minly:
                        ly = minly

                if (lx == 0) & (ly==0):
                        lx=lx+0.001

                l = np.sqrt(lx*lx+ly*ly)

                xstart = -B[0]+bx
                ystart = -B[1]+by
                handles=[]
                for i in range(0,2*N-1):
                    for j in range(0,2*M-1):
                        x = xstart+i*dxspacing+dxspacing
                        y = ystart+j*dyspacing+dyspacing
                        A = self.env.drawarrow(array((x,y,self.ZPOS_ARROW)),array((x+lx,y+ly,self.ZPOS_ARROW)),linewidth=0.08*l,color=array((1,0,0)))
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
