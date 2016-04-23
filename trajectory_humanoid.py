import abc
import time
import numpy as np
from trajectory import *

class TrajectoryHumanoid(Trajectory):

        @classmethod
        def from_ravetraj(cls, ravetraj, robot):

                M = ravetraj.GetNumWaypoints()
                N = len(robot.GetActiveDOFValues())
                active_dofs = robot.GetActiveConfigurationSpecification()

                W=np.zeros((N,M))

                for i in range(0,M):
                        W[:,i] = ravetraj.GetWaypoint(i, active_dofs)

                return cls(W)
                ###[qlimL,qlimU]=robot.GetActiveDOFLimits()
