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
        cells = []
        forces = []
        def __init__(self):
                self.env=Environment()
                self.env.SetViewer('qtcoin')
                self.env.Reset()
                self.physics = RaveCreatePhysicsEngine(self.env,'ode')
                self.physics.SetGravity(array((0,0,-9.81)))
                self.env.SetPhysicsEngine(self.physics)

        def setrobotenv(self,robot_xml,env_xml):
                self.env_xml = env_xml
                self.robot_xml = robot_xml
                self.env.Add(self.env.ReadRobotXMLFile(self.robot_xml))
                self.env.Load(self.env_xml)

        def GetCellsAll(self):
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


        def DisplayForces(self):
                if self.forces is None:
                        self.forces = GetForces()
                if self.cells is None:
                        self.cells = GetCells()

                assert(len(self.forces)==len(self.GetCells()))

                for i in range(0,len(self.cells)):

                        G1 = self.cells[i].GetGeometries()[0]
                        B = G1.GetBoxExtents()
                        T = G1.GetTransform()

                        P = array(((T[0,3]-B[0],T[1,3]-B[1],B[2]+T[2,3]), \
                                (T[0,3]+B[0],T[1,3]-B[1],B[2]+T[2,3]), \
                                (T[0,3]+B[0],T[1,3]+B[1],B[2]+T[2,3]), \
                                (T[0,3]-B[0],T[1,3]+B[1],B[2]+T[2,3]),\
                                (T[0,3]-B[0],T[1,3]-B[1],B[2]+T[2,3])))

                        A=self.env.drawlinestrip(points=P,linewidth=5.0,colors=array(((1,1,1))))
                        self.handles.append(A)
                        

                        bx = T[0,3]
                        by = T[1,3]
                        N = int(B[0])
                        M = int(B[1])
                        dxspacing = B[0]/N
                        dyspacing = B[1]/M

                        Fx = self.forces[i][0]
                        Fy = self.forces[i][1]
                        if sqrt(Fx*Fx+Fy*Fy) < 0.001:
                                continue

                        ddx = dxspacing-dxspacing/4
                        dx = ddx*Fx/sqrt(Fx*Fx+0.001)
                        ddy = dyspacing-dyspacing/4
                        dy = ddy*Fy/sqrt(Fy*Fy+0.001)
                        #dx = min(dxspacing-dxspacing/10,sqrt(Fx*Fx))
                        #dy = min(dyspacing-dyspacing/10,sqrt(Fy*Fy))
                        if (dx == 0) & (dy==0):
                                dx=dx+0.001

                        l = sqrt(dx*dx+dy*dy)
                        print dx,dy,N,M,l

                        print dxspacing,dyspacing
                        xstart = -B[0]+bx
                        ystart = -B[1]+by
                        zval=0.1
                        for i in range(0,2*N-1):
                            for j in range(0,2*M-1):
                                x = xstart+i*dxspacing+dxspacing
                                y = ystart+j*dyspacing+dyspacing
                                A = self.env.drawarrow(array((x,y,zval)),array((x+dx,y+dy,zval)),linewidth=0.08*l,color=array((1,0,0,0.5)))
                                self.handles.append(A)





if __name__ == "__main__":
        env = ForceEnvironment()
        xmlenv='environments/the_stream.env.xml'
        xmlrobot='robots/pointrobot.robot.xml'
        env.setrobotenv(xmlrobot,xmlenv)
