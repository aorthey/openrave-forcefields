

import time
import numpy as np
from scipy.interpolate import interp1d,splev,splrep,splprep
from scipy.misc import derivative


class TrajectoryAnalyzer():
        tau = []
        tau_step = []

        handles = []
        Whandle = []
        degree = 3

        def __init__(self, W):
                self.tau = self.WaypointsToTrajectory(W)

        def WaypointsToTrajectory(self,W):
                Nwaypoints = W.shape[1]
                if Nwaypoints<=4:
                        self.degree=1
                else:
                        self.degree=3
                tau,tmp = splprep(W,k=self.degree)
                return tau

        def GetForcesAtWaypoints(self, W, env):
                Nwaypoints = W.shape[1]
                F = np.zeros((3,Nwaypoints))
                for i in range(0,Nwaypoints):
                        pt = np.array(((W[0,i],W[1,i],-0.1,0.001)))
                        F[:,i] = env.GetForceAtX(pt)
                return F

        def TrajectoryToWaypoints(self, tau, N = 100):
                pts = np.zeros((3,N))
                dpts = np.zeros((3,N))
                ctr = 0
                for t in np.linspace(0.0, 1.0, num=N):
                        [f0,df0] = self.funcEval(tau,t)
                        if f0.shape[0] == 2:
                                pts[:,ctr] = np.array(((f0[0],f0[1],0.15)))
                                dpts[:,ctr] = np.array(((df0[0],df0[1],0.15)))
                        else:
                                pts[:,ctr] = f0
                                dpts[:,ctr] = df0
                        ctr = ctr+1
                return [pts,dpts]

        def funcEval(self, tau, t):
                f0 = splev(t,tau)
                df0 = splev(t,tau,der=1)
                f0 = np.array(f0)
                df0 = np.array(df0)
                return [f0,df0]


        ptsize = 0.05
        def DrawRedPoint(self,env,X):
                self.handles.append(env.env.plot3(points=X,
                                   pointsize=self.ptsize,
                                   colors=np.array(((0.8,0.0,0.0,0.9))),
                                   drawstyle=1))
        def DrawGreenPoint(self,env,X):
                self.handles.append(env.env.plot3(points=X,
                                   pointsize=self.ptsize,
                                   colors=np.array(((0.0,0.8,0.0,0.9))),
                                   drawstyle=1))

        
        def GetLengthOfTrajectory(self, tau):
                dt = 0.05

                dd = 0.0
                for t in np.linspace(0.0, 1.0-dt, num=1000):
                        [ft,df0] = self.funcEval(tau,t)
                        [ftdt,df0] = self.funcEval(tau,t+dt)
                        dd = dd + np.linalg.norm(ft-ftdt)

                return dd


        def DrawWaypoints(self, Win, env):
                tau = self.WaypointsToTrajectory(Win)

                L = self.GetLengthOfTrajectory(tau)
                dt = 0.05
                Nwaypoints = int(L/dt)
                print Nwaypoints

                [W,dW] = self.TrajectoryToWaypoints(tau,N=Nwaypoints)

                Nwaypoints = W.shape[1]
                self.Whandle = []
                for i in range(0,Nwaypoints):
                        pt = np.array(((W[0,i],W[1,i],0.1)))
                        self.Whandle.append(env.env.plot3(points=pt,
                                           pointsize=self.ptsize,
                                           colors=np.array(((0.7,0.2,0.7,0.9))),
                                           drawstyle=1))

        def analyze(self,env):

                [W,dW] = self.TrajectoryToWaypoints(self.tau, N = 100)
                Nwaypoints = W.shape[1]
                F = self.GetForcesAtWaypoints(W, env)

                for i in range(0,Nwaypoints):
                        d = np.linalg.norm(F[:,i])
                        pt = np.array(((W[0,i],W[1,i],0.1)))
                        if d > 0.2:
                                self.DrawRedPoint(env, pt)
                        else:
                                self.DrawGreenPoint(env, pt)

        def GetFirstInfeasibleWaypoint(self, W, dW, F):
                Nwaypoints = W.shape[1]
                for i in range(0,Nwaypoints):
                        d = np.linalg.norm(F[:,i])
                        pt = np.array(((W[0,i],W[1,i],0.1)))
                        if d > 0.2:
                                return i

                ##no deform necessary


        def deform(self, env):

                [W,dW] = self.TrajectoryToWaypoints(self.tau, N = 100)

                Nwaypoints = W.shape[1]

                F = self.GetForcesAtWaypoints(W, env)

                i_infeasible = self.GetFirstInfeasibleWaypoint(W, dW, F)

                dW_direction = dW[:,i_infeasible]/(np.linalg.norm(dW[:,i_infeasible]))

                self.DrawRedPoint(env,W[:,i_infeasible])
                self.DrawGreenPoint(env,W[:,i_infeasible]+dW_direction)
                self.DrawWaypoints(W, env)

                qmax = 1.0

                for p in range(0,20):
                        ##velocity applied to each iteration
                        dV = 0.5*p


                        dist_until_captured = dV*dV/(2*qmax)
                        print dist_until_captured,dW_direction
                        #self.DrawGreenPoint(env,W[:,i_infeasible]-dist_until_captured*dW_direction)

                        Wupdate = np.zeros((3,Nwaypoints))
                        dd = 0.0

                        i = i_infeasible
                        while dd < dist_until_captured and i > 1:
                                ##project W[:,i-1] onto W - distW*dWdirection
                                scale = -(i-i_infeasible)
                                dd = scale * dist_until_captured/100
                                #dd = np.linalg.norm(W[:,i_infeasible] - W[:,i-1])
                                W_line = W[:,i_infeasible] - dd*dW_direction
                                Wupdate[:,i-1] = W_line - W[:,i-1]
                                i = i - 1

                        print "projected",-(i-i_infeasible),"points"
                        W = W + 0.1*Wupdate 
                        print W[:,-1]

                        self.DrawWaypoints(W, env)
                        time.sleep(0.1)

