import numpy as np
from scipy.interpolate import interp1d,splev,splrep,splprep
from scipy.misc import derivative

class Trajectory():

        traj = []

        def __init__(self, trajectory):
                self.traj = trajectory 

        @classmethod
        def from_ravetraj(cls, ravetraj):
                N = ravetraj.GetNumWaypoints()
                W=[]
                for i in range(0,N):
                        w = np.array((ravetraj.GetWaypoint(i)[0],ravetraj.GetWaypoint(i)[1],ravetraj.GetWaypoint(i)[2]))
                        W.append((w))
                W = np.array(W).T
                return cls.from_waypoints(W)

        @classmethod
        def from_waypoints(cls, W):
                Nwaypoints = W.shape[1]
                if Nwaypoints<=4:
                        degree=1
                else:
                        degree=3
                trajectory,tmp = splprep(W,k=degree)
                return cls(trajectory)

        def evaluateAt(self, t):
                f0 = splev(t,self.traj)
                df0 = splev(t,self.traj,der=1)
                f0 = np.array(f0)
                df0 = np.array(df0)
                return [f0,df0]

