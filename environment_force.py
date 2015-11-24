import abc
import openravepy
if not __openravepy_build_doc__:
    from openravepy import *
    from numpy import *

class ForceEnvironment():
        __metaclass__ = abc.ABCMeta
        env_xml=''
        robot_xml=''
        env=''
        robot=''
        handles = []
        cells = None
        forces = None
        def __init__(self):
                self.env=Environment()
                self.env.SetViewer('qtcoin')
                self.env.Reset()
                with self.env:
                        self.physics = RaveCreatePhysicsEngine(self.env,'ode')
                        self.physics.SetGravity(array((0,0,-9.81)))
                        self.env.SetPhysicsEngine(self.physics)
                        #self.recorder = RaveCreateModule(self.env,'viewerrecorder')
                        #self.env.AddModule(self.recorder,'')

        def setrobotenv(self,robot_xml,env_xml):
                with self.env:
                        self.env_xml = env_xml
                        self.robot_xml = robot_xml
                        self.env.Add(self.env.ReadRobotXMLFile(self.robot_xml))
                        self.env.Load(self.env_xml)

        def GetCellsAll(self):
                with self.env:
                        B = self.env.GetBodies()[1]
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
                        [x,y]=self.RobotGetInitialPosition()
                        H = self.robot.GetTransform()
                        H[0,3] = x
                        H[1,3] = y
                        self.robot.SetTransform(H)
                        [xg,yg]=self.RobotGetGoalPosition()
                        self.handles.append(self.env.plot3(points=array(((xg,yg,0.05))),
                                           pointsize=0.15,
                                           colors=array(((0.0,1.0,0.0))),
                                           drawstyle=1))
                        return self.robot

        def DrawBorderAroundCell(self, cell):
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

                        handle=self.env.drawlinestrip(points=P,linewidth=5.0,colors=array(((1,1,1))))
                        return handle

        def DrawForceArrowsInCell(self, cell, force):
                with self.env:
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
                        N = int(math.ceil(B[0]))+1
                        M = int(math.ceil(B[1]))+1
                        dxspacing = B[0]/N
                        dyspacing = B[1]/M

                        if sqrt(Fx*Fx+Fy*Fy) < 0.001:
                                return None

                        ### lx,ly == direction vector of line
                        maxlx = dxspacing-dxspacing/4
                        minlx = -maxlx
                        maxly = dyspacing-dyspacing/4
                        minly = -maxly
                        lx = Fx
                        ly = Fy

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

                        l = sqrt(lx*lx+ly*ly)

                        xstart = -B[0]+bx
                        ystart = -B[1]+by
                        zval=0.1
                        handles=[]
                        for i in range(0,2*N-1):
                            for j in range(0,2*M-1):
                                x = xstart+i*dxspacing+dxspacing
                                y = ystart+j*dyspacing+dyspacing
                                A = self.env.drawarrow(array((x,y,zval)),array((x+lx,y+ly,zval)),linewidth=0.08*l,color=array((1,0,0,0.5)))
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

                assert(len(self.forces)==len(self.GetCells()))

                for i in range(0,len(self.cells)):
                        h = self.DrawBorderAroundCell(self.cells[i])
                        self.handles.append(h)
                        h = self.DrawForceArrowsInCell(self.cells[i], self.forces[i])
                        self.handles.append(h)

if __name__ == "__main__":
        env = ForceEnvironment()
        xmlenv='environments/the_stream.env.xml'
        xmlrobot='robots/pointrobot.robot.xml'
        env.setrobotenv(xmlrobot,xmlenv)