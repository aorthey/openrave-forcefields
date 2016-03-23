import abc
import sys
import time
import numpy as np
from trajectory import *
import copy

class Deformation():
        __metaclass__ = abc.ABCMeta

        env = []
        traj_ori = []
        traj_deformed = []
        traj_current = []

        traj_display = []
        handle = []
        forcehandle = []


        def __init__(self, trajectory, environment):
                self.env = environment

                ## traj_ori : original trajectory before any calls
                ## traj_current : the current deformation
                ## traj_deformed : the working copy of current for the next iteration
                ## traj_display : the trajectory which we display at the moment
                self.traj_ori = trajectory 
                self.traj_current = copy.copy(trajectory)
                self.traj_deformed = copy.copy(trajectory)
                self.traj_display = copy.copy(trajectory)

        @abc.abstractmethod 
        def deform_onestep(self):
                pass

        def deform(self, N_iter = 1):
                for i in range(0,N_iter):
                        self.deform_onestep()
                        self.traj_current = copy.copy(self.traj_deformed)

        def GetForcesAtWaypoints(self, W):
                Ndim = W.shape[0]
                Nwaypoints = W.shape[1]
                F = np.zeros((Ndim,Nwaypoints))
                for i in range(0,Nwaypoints):
                        pt = np.array(((W[0,i],W[1,i],-0.1,0.001)))
                        F[0:3,i] = self.env.GetForceAtX(pt)
                return F

        def draw_forces_at_waypoints(self, W, F):
                for i in range(0,W.shape[1]):
                        df = np.linalg.norm(F[:,i])
                        if df > 0.3:
                                nF = 0.3*F[:,i]/df
                        else:
                                nF = F[:,i]
                        W1 = np.array((W[0,i],W[1,i],W[2,i]))
                        W2 = np.array((W[0,i]+nF[0],W[1,i]+nF[1],W[2,i]+nF[2]))
                        P = np.array(((W1),(W2)))
                        ls = self.traj_current.linsize
                        h=self.env.env.drawlinestrip(points=P,linewidth=ls,colors=np.array(((0.8,0.2,0.2,0.9))))
                        self.forcehandle.append(h)


        def GetFirstInfeasibleWaypoint(self, W, dW, F):
                Nwaypoints = W.shape[1]
                for i in range(0,Nwaypoints):
                        d = np.linalg.norm(F[:,i])
                        pt = np.array(((W[0,i],W[1,i],0.1)))
                        if d > 0.2:
                                return i

        def draw_deformation(self):
                M = 20
                L = self.traj_current.get_length()
                dt = 0.05
                Nwaypoints = int(L/dt)

                print self.traj_display.info()
                print self.traj_current.info()

                [W0,dW] = self.traj_display.get_waypoints(N=Nwaypoints) 
                [W1,dW] = self.traj_current.get_waypoints(N=Nwaypoints) 

                self.forcehandle = []
                for i in range(0,M):
                        k = float(i)/float(M)
                        Wk = (1-k)*W0 + k*W1
                        self.traj_current.new_from_waypoints(Wk)
                        self.handle = self.traj_current.draw(self.env, keep_handle=False)
                        if i == 0:
                                raw_input('Press <ENTER> to draw deformation.')
                        time.sleep(0.1)

                self.traj_display = copy.copy(self.traj_current)


        def draw_trajectory_original(self):
                self.traj_ori.draw(self.env)

        def draw_trajectory_deformed(self):
                self.traj_current.draw(self.env)

