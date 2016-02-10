from environment_force import *
import numpy as np
import openravepy

class EnvironmentTheRay(ForceEnvironment):
        def __init__(self):
                ForceEnvironment.__init__(self)
                xmlenv='environments/the_ray.env.xml'
                xmlrobot='robots/pointrobot.robot.xml'
                self.setrobotenv(xmlrobot,xmlenv)

        def GetCells(self):
                C = self.GetCellsAll()
                ## do not use the first link, because it is a background
                self.cells = C[0:3]
                return self.cells

        def GetForces(self):
                ##
                self.forces = np.array((0.0,0.0,0.0))
                self.forces = np.vstack([self.forces,(0.0,-1.5,0.0)])
                self.forces = np.vstack([self.forces,(0.0,0.0,0.0)])
                return self.forces

        def RobotGetInitialPosition(self):
                return [-2.5,0.0,0.15,0,0,0,0,0]
                #return [3.0,3.0]

        def RobotGetGoalPosition(self):
                return [-4.5,-0.0,0.15,0,0,0,0,0]




if __name__ == "__main__":
        env = EnvironmentTheRay()
