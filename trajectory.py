import numpy as np
from scipy.interpolate import interp1d,splev,splrep,splprep
from scipy.misc import derivative

class Trajectory():

        traj = []
        handle = []
        ptsize = 0.05
        trajectory_color = np.array((0.7,0.2,0.7,0.9))

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

        def evaluate_at(self, t):
                f0 = splev(t,self.traj)
                df0 = splev(t,self.traj,der=1)
                f0 = np.array(f0)
                df0 = np.array(df0)
                return [f0,df0]

        def get_length(self):
                dt = 0.05
                dd = 0.0
                for t in np.linspace(0.0, 1.0-dt, num=1000):
                        [ft,df0] = self.evaluate_at(t)
                        [ftdt,df0] = self.evaluate_at(t+dt)
                        dd = dd + np.linalg.norm(ft-ftdt)
                return dd

        def get_waypoints(self, N = 100):
                pts = np.zeros((3,N))
                dpts = np.zeros((3,N))
                ctr = 0
                for t in np.linspace(0.0, 1.0, num=N):
                        [f0,df0] = self.evaluate_at(t)
                        if f0.shape[0] == 2:
                                pts[:,ctr] = np.array(((f0[0],f0[1],0.15)))
                                dpts[:,ctr] = np.array(((df0[0],df0[1],0.15)))
                        else:
                                pts[:,ctr] = f0
                                dpts[:,ctr] = df0
                        ctr = ctr+1
                return [pts,dpts]

        def draw(self, env, keep_handle=True):
                L = self.get_length()
                dt = 0.05
                Nwaypoints = int(L/dt)

                [W,dW] = self.get_waypoints(N=Nwaypoints)

                tmp_handle = []
                for i in range(0,Nwaypoints):
                        pt = np.array(((W[0,i],W[1,i],0.1)))
                        tmp_handle.append(env.env.plot3(points=pt,
                                           pointsize=self.ptsize,
                                           colors=self.trajectory_color,
                                           drawstyle=1))
                if keep_handle:
                        self.handle = tmp_handle
                else:
                        return tmp_handle

        def draw_delete(self, env):
                self.handle = []
