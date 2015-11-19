from environment import *
import openravepy

class EnvironmentTheCounterStream(ForceEnvironment):
        def __init__(self):
                ForceEnvironment.__init__(self)
                xmlenv='environments/the_counterstream.env.xml'
                xmlrobot='robots/pointrobot.robot.xml'
                self.setrobotenv(xmlrobot,xmlenv)

        def GetCells(self):
                C = self.GetCellsAll()
                ## do not use the first link, because it is a background
                self.cells = C[1:]
                return self.cells

        def GetForces(self):
                ##
                self.forces=[]
                self.forces.append(numpy.array((0.0,0.0,0.0)))
                self.forces.append(numpy.array((0.0,0.0,0.0)))
                self.forces.append(numpy.array((3.0,0.0,0.0)))
                self.forces.append(numpy.array((-3.0,0.0,0.0)))
                self.forces.append(numpy.array((0.0,0.0,0.0)))
                return self.forces

        def RobotGetInitialPosition(self):
                return [-2.0,1.0]

        def RobotGetGoalPosition(self):
                return [5.0,-3.0]




if __name__ == "__main__":
        env = EnvironmentTheStream()
