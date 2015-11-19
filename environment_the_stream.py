from environment import *
import openravepy

class EnvironmentTheStream(ForceEnvironment):
        def __init__(self):
                ForceEnvironment.__init__(self)
                xmlenv='environments/the_stream.env.xml'
                xmlrobot='robots/pointrobot.robot.xml'
                self.setrobotenv(xmlrobot,xmlenv)

        def GetCells(self):
                C = self.GetCellsAll()
                ## do not use the first link, because it is a background
                self.cells = C[1:]
                return self.cells

        def GetForces(self):
                ##
                F1=numpy.array((1.0,0.0,0.0))
                F2=numpy.array((0.0,-3.0,0.0))
                F3=numpy.array((0.0,0.0,0.0))
                self.forces = [F1,F2,F3]
                return self.forces

        def RobotGetInitialPosition(self):
                return [-4.0,-2.5]

        def RobotGetGoalPosition(self):
                return [4.0,2.0]




if __name__ == "__main__":
        env = EnvironmentTheStream()
