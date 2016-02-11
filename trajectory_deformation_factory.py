import abc
import time
import numpy as np

from trajectory_utils import *

class TrajectoryDeformation():
        __metaclass__ = abc.ABCMeta

        env = []
        traj_ori = []
        traj_deformed = []
        handle = []

        ptsize = 0.05

        def __init__(self, trajectory, environment):
                self.env = environment
                self.traj_ori = trajectory 
                self.traj_deformed = trajectory 

        @classmethod
        def from_ravetraj(cls, ravetraj, environment):
                N = ravetraj.GetNumWaypoints()
                W=[]
                for i in range(0,N):
                        w = np.array((ravetraj.GetWaypoint(i)[0],ravetraj.GetWaypoint(i)[1],ravetraj.GetWaypoint(i)[2]))
                        W.append((w))
                W = np.array(W).T
                trajectory = WaypointsToTrajectory(W)
                return cls(trajectory, environment)

        @classmethod
        def from_waypoints(cls, W, environment):
                trajectory = WaypointsToTrajectory(W)
                return cls(trajectory, environment)

        @abc.abstractmethod 
        def deform_onestep(self):
                pass

        def deform(self, N_iter = 1):
                for i in range(0,N_iter):
                        self.deform_onestep()


        def GetForcesAtWaypoints(self, W):
                Nwaypoints = W.shape[1]
                F = np.zeros((3,Nwaypoints))
                for i in range(0,Nwaypoints):
                        pt = np.array(((W[0,i],W[1,i],-0.1,0.001)))
                        F[:,i] = self.env.GetForceAtX(pt)
                return F

        def GetFirstInfeasibleWaypoint(self, W, dW, F):
                Nwaypoints = W.shape[1]
                for i in range(0,Nwaypoints):
                        d = np.linalg.norm(F[:,i])
                        pt = np.array(((W[0,i],W[1,i],0.1)))
                        if d > 0.2:
                                return i


        def draw_deformation(self):
                M = 20
                L = GetLengthOfTrajectory(self.traj_ori)
                dt = 0.05
                Nwaypoints = int(L/dt)
                print Nwaypoints

                [W0,dW] = TrajectoryToWaypoints(self.traj_ori,N=Nwaypoints)
                [W1,dW] = TrajectoryToWaypoints(self.traj_deformed,N=Nwaypoints)

                
                for i in range(0,M):
                        k = float(i)/float(M)
                        Wk = (1-k)*W0 + k*W1
                        print k
                        handle_tmp = []
                        for i in range(0,Nwaypoints):
                                pt = np.array(((Wk[0,i],Wk[1,i],0.1)))
                                handle_tmp.append(self.env.env.plot3(points=pt,
                                                   pointsize=self.ptsize,
                                                   colors=np.array(((0.7,0.2,0.7,0.9))),
                                                   drawstyle=1))
                        self.handle = handle_tmp
                        time.sleep(0.05)

        def draw_waypoints(self, Win):
                Nwaypoints = Win.shape[1]

                self.handle = []
                for i in range(0,Nwaypoints):
                        pt = np.array(((Win[0,i],Win[1,i],0.1)))
                        self.handle.append(self.env.env.plot3(points=pt,
                                           pointsize=self.ptsize,
                                           colors=np.array(((0.1,0.7,0.1,0.9))),
                                           drawstyle=1))

        def draw_specific_trajectory(self, tau):
                L = GetLengthOfTrajectory(tau)
                dt = 0.05
                Nwaypoints = int(L/dt)
                print Nwaypoints

                [W,dW] = TrajectoryToWaypoints(tau,N=200)

                Nwaypoints = W.shape[1]

                self.handle = []
                for i in range(0,Nwaypoints):
                        pt = np.array(((W[0,i],W[1,i],0.1)))
                        self.handle.append(self.env.env.plot3(points=pt,
                                           pointsize=self.ptsize,
                                           colors=np.array(((0.7,0.2,0.7,0.9))),
                                           drawstyle=1))

        def draw_trajectory_original(self):
                self.draw_specific_trajectory(self.traj_ori)

        def draw_trajectory_deformed(self):
                self.draw_specific_trajectory(self.traj_deformed)



        #ptsize = 0.05
        #def DrawRedPoint(self,env,X):
        #        self.handles.append(env.env.plot3(points=X,
        #                           pointsize=self.ptsize,
        #                           colors=np.array(((0.8,0.0,0.0,0.9))),
        #                           drawstyle=1))

        #def DrawRedArrow(self,env,X,dX):
        #        A = env.env.drawarrow(X,X+dX,linewidth=0.02,color=np.array((1,0,0)))
        #        self.handles.append(A)

        #def DrawGreenPoint(self,env,X):
        #        self.handles.append(env.env.plot3(points=X,
        #                           pointsize=self.ptsize,
        #                           colors=np.array(((0.0,0.8,0.0,0.9))),
        #                           drawstyle=1))

