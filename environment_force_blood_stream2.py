from environment_force import *
import numpy as np
import openravepy

class EnvironmentBloodStream2(ForceEnvironment):
        def __init__(self):
                ForceEnvironment.__init__(self)
                xmlenv='environments/the_blood_stream2.env.xml'
                xmlrobot='robots/virus.robot.xml'
                self.setrobotenv(xmlrobot,xmlenv)

        def GetCells(self):
                C = self.GetCellsAll()
                self.cells = C[0:3]
                return self.cells

        def GetForces(self):
                self.forces = np.array((0.0,-3.5,0.0))
                self.forces = np.vstack([self.forces,(0.0,0.0,0.0)])
                self.forces = np.vstack([self.forces,(0.0,0.0,0.0)])

                return self.forces

        def RobotGetInitialPosition(self):
                return [0.0,-4.0,0.1,pi/2,0,0,0,0]

        def RobotGetGoalPosition(self):
                return [0.0,3.0,0.1,pi/2,0,0,0,0]

if __name__ == "__main__":
        env = EnvironmentBloodStream2()
        env.DisplayForces()
