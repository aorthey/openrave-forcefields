import abc
import time
import numpy as np
from trajectory import *

class TrajectoryHumanoid(Trajectory):

        @classmethod
        def from_ravetraj(cls, ravetraj):
                rave_traj = ravetraj
                N = ravetraj.GetNumWaypoints()
                W=[]
                for i in range(0,N):
                        w = np.array((ravetraj.GetWaypoint(i)[0],ravetraj.GetWaypoint(i)[1],ravetraj.GetWaypoint(i)[2],ravetraj.GetWaypoint(i)[3]))
                        W.append((w))
                W = np.array(W).T
                return cls(W)
