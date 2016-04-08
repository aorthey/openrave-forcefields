from environment_force import *
import numpy as np
import openravepy

class EnvironmentHumanoid(ForceEnvironment):
        def __init__(self):
                ForceEnvironment.__init__(self)
                xmlenv='environments/plane.env.xml'
                xmlrobot='misc/escher_v2.kinbody.urdf'
                self.setrobotenv(xmlrobot,xmlenv)

        def GetCells(self):
                C = self.GetCellsAll()
                self.cells = C[0:1]
                return self.cells

        def GetForces(self):
                ##
                self.forces = np.array((0.0,0.0,0.0))
                self.forces = np.vstack([self.forces,(0.0,0.0,0.0)])
                return self.forces

        def RobotGetInitialPosition(self):
                return [-2.5,0.0,0.1,-pi,0,0,0,0]

        def RobotGetGoalPosition(self):
                return [-4.5,-0.0,0.1,-pi,0,0,0,0]

        def GetRobot(self):
                with self.env:
                        self.robot = self.env.GetRobots()[0]
                        return self.robot

if __name__ == "__main__":
        env = EnvironmentHumanoid()
