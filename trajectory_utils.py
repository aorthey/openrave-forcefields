import numpy as np
from scipy.interpolate import interp1d,splev,splrep,splprep
from scipy.misc import derivative


def funcEval(tau, t):
        f0 = splev(t,tau)
        df0 = splev(t,tau,der=1)
        f0 = np.array(f0)
        df0 = np.array(df0)
        return [f0,df0]

def GetLengthOfTrajectory(traj):
        dt = 0.05

        dd = 0.0
        for t in np.linspace(0.0, 1.0-dt, num=1000):
                [ft,df0] = funcEval(traj,t)
                [ftdt,df0] = funcEval(traj,t+dt)
                dd = dd + np.linalg.norm(ft-ftdt)
        return dd

def WaypointsToTrajectory(W):
        Nwaypoints = W.shape[1]
        if Nwaypoints<=4:
                degree=1
        else:
                degree=3
        tau,tmp = splprep(W,k=degree)
        return tau

def TrajectoryToWaypoints(tau, N = 100):
        pts = np.zeros((3,N))
        dpts = np.zeros((3,N))
        ctr = 0
        for t in np.linspace(0.0, 1.0, num=N):
                [f0,df0] = funcEval(tau,t)
                if f0.shape[0] == 2:
                        pts[:,ctr] = np.array(((f0[0],f0[1],0.15)))
                        dpts[:,ctr] = np.array(((df0[0],df0[1],0.15)))
                else:
                        pts[:,ctr] = f0
                        dpts[:,ctr] = df0
                ctr = ctr+1
        return [pts,dpts]

