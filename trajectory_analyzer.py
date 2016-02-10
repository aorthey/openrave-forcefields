

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

        def DrawRedArrow(self,env,X,dX):
                A = env.env.drawarrow(X,X+dX,linewidth=0.02,color=np.array((1,0,0)))
                self.handles.append(A)

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
                Nwaypoints = Win.shape[1]

                self.Whandle = []
                for i in range(0,Nwaypoints):
                        pt = np.array(((Win[0,i],Win[1,i],0.1)))
                        self.Whandle.append(env.env.plot3(points=pt,
                                           pointsize=self.ptsize,
                                           colors=np.array(((0.7,0.2,0.7,0.9))),
                                           drawstyle=1))

        def draw(self, env, keep=False):
                L = self.GetLengthOfTrajectory(self.tau)
                dt = 0.05
                Nwaypoints = int(L/dt)
                print Nwaypoints

                [W,dW] = self.TrajectoryToWaypoints(self.tau,N=200)

                Nwaypoints = W.shape[1]

                if not keep:
                        self.Whandle = []
                for i in range(0,Nwaypoints):
                        pt = np.array(((W[0,i],W[1,i],0.1)))
                        self.Whandle.append(env.env.plot3(points=pt,
                                           pointsize=self.ptsize,
                                           colors=np.array(((0.7,0.2,0.7,0.9))),
                                           drawstyle=1))

        def draw_deformation(self, env, tau0, tau1):

                M = 20

                L = self.GetLengthOfTrajectory(self.tau)
                dt = 0.05
                Nwaypoints = int(L/dt)
                print Nwaypoints

                [W0,dW] = self.TrajectoryToWaypoints(tau0,N=Nwaypoints)
                [W1,dW] = self.TrajectoryToWaypoints(tau1,N=Nwaypoints)

                
                for i in range(0,M):
                        k = float(i)/float(M)
                        Wk = (1-k)*W0 + k*W1
                        print k
                        self.Whandle = []
                        for i in range(0,Nwaypoints):
                                pt = np.array(((Wk[0,i],Wk[1,i],0.1)))
                                self.Whandle.append(env.env.plot3(points=pt,
                                                   pointsize=self.ptsize,
                                                   colors=np.array(((0.7,0.2,0.7,0.9))),
                                                   drawstyle=1))
                        time.sleep(0.1)




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

        def gaussian(self, a, b, c, t):
                return a*np.exp(-((t-b)*(t-b))/(2*c*c))

        def lambda_position(self, t0, tidx):
                a = 0.25
                c = 15.0
                return self.gaussian(a,0,c,abs(t0-tidx))

        def lambda_endpoint(self, t0, endT):
                a = 1.0
                c = 0.5
                return self.gaussian(a,0,c,abs(t0-endT))

        def lambda_stretch(self, t0, tidx):
                ## a has to be 1.0 to make sure that we stretch the last
                ##point the correct amount
                a = 1.0 
                c = 10.0
                return self.gaussian(a,0,c,abs(t0-tidx))

        ## idx: index at which waypoint will be shifted
        ## W  : array of waypoints
        ## dF : direction of force applied at W[:,idx]
        ## 
        ## return: a force direction to each waypoint
        def GetForcePosition(self, idx, W, dF):

                N = W.shape[1]
                dW = np.zeros((3,N))
                for i in range(0,N):
                        dW[:,i] = self.lambda_position(i, idx)*dF

                return dW

        def GetForceCellCorrection(self, W, Cells, dF):

                N = W.shape[1]
                dW = np.zeros((3,N))
                for i in range(0,N):
                        dW[:,i] = self.lambda_position(i, idx)*dF
                return dW

        def GetForceStretch(self, idx, W, dF):

                assert(idx>=1)

                N = W.shape[1]
                dW = np.zeros((3,N))

                for i in range(0,idx):
                        dW[:,i] = self.lambda_stretch(i, idx)*dF
                return dW

        ### remove any updates on the endpoint, keep them static
        def GetForceEndpointCorrection(self, Wupdate, W, Winit, Wgoal):

                N = Wupdate.shape[1]
                dW = np.zeros((3,N))

                assert(W.shape[1]==Wupdate.shape[1])

                for i in range(0,N):
                        dW[:,i] = dW[:,i] + self.lambda_endpoint(i, 0)*(-Wupdate[:,0])
                        dW[:,i] = dW[:,i] + self.lambda_endpoint(i, N-1)*(-Wupdate[:,-1])

                if np.linalg.norm(dW[:,-1])>0.001:
                        print "update correction not successful!"
                        sys.exit(0)

                dWinit = Winit[0:3] - W[:,0]
                dWgoal = Wgoal[0:3] - W[:,-1]

                for i in range(0,N):
                        dW[:,i] = dW[:,i] + self.lambda_endpoint(i, 0)*(dWinit)
                        dW[:,i] = dW[:,i] + self.lambda_endpoint(i, N-1)*(dWgoal)

                Wcorrection = Wupdate + dW
                return Wcorrection


        def deform_onestep(self, env):
                Ninsert = 50
                lambda_insert = 0.03

                traj_old = self.tau
                [Wori,dWori] = self.TrajectoryToWaypoints(self.tau, N = 100)
                Nwaypoints = Wori.shape[1]
                F = self.GetForcesAtWaypoints(Wori, env)

                idx_invalid = self.GetFirstInfeasibleWaypoint(Wori, dWori, F)
                dW_invalid = dWori[:,idx_invalid]/(np.linalg.norm(dWori[:,idx_invalid]))
                print "############################################################"
                print "INVALID POINT FOUND"
                print idx_invalid
                print Wori[:,idx_invalid], dW_invalid

                #self.DrawRedPoint(env, Wori[:,idx_invalid])
                self.DrawRedArrow(env, Wori[:,idx_invalid], dW_invalid)

                Wdown = Wori[:,0:idx_invalid]
                Wup = Wori[:,idx_invalid:]

                W = np.zeros((3,Ninsert+Nwaypoints))

                W[:,0:idx_invalid] = Wdown
                W[:,idx_invalid+Ninsert:] = Wup

                Wmiddle = np.zeros((3,Ninsert))

                for i in range(0,Ninsert):
                        Wstretch = (Ninsert-(i))*lambda_insert*dW_invalid
                        Wmiddle[:,i] = Wori[:,idx_invalid] - Wstretch

                W[:,idx_invalid:idx_invalid+Ninsert] = Wmiddle

                Wupdate_stretch = self.GetForceStretch(idx_invalid, W, -(Ninsert+1)*lambda_insert*dW_invalid)
                Wupdate_position = self.GetForcePosition(idx_invalid+Ninsert, W, -F[:,idx_invalid])

                Wupdate = Wupdate_position+Wupdate_stretch

                Winit=env.RobotGetInitialPosition()
                Wgoal=env.RobotGetGoalPosition()
                Wupdate_correction = self.GetForceEndpointCorrection(Wupdate, W, Winit, Wgoal)

                W = W + Wupdate_correction

                print "UPDATE WAYPOINTS:"
                print np.around(Wupdate_correction[:,0],decimals=2)
                print np.around(Wupdate_correction[:,-1],decimals=2)

                #print np.around(W[:,:idx_invalid+Ninsert+2].T,decimals=2)
                #print np.around(W[:,idx_invalid+Ninsert-2:].T,decimals=2)
                self.tau = self.WaypointsToTrajectory(W)
                [f0,df0] = self.funcEval(self.tau,0)
                print "FUNCTION WAYPOINTS VS ORIGINAL WAYPOINTS"
                print np.around(f0,decimals=2), W[:,0]
                [f1,df1] = self.funcEval(self.tau,1)
                print np.around(f1,decimals=2), W[:,-1]

                #self.DrawWaypoints(W, env)
                print "############################################################"
                #for i in range(0,20):
                        #self.DrawWaypoints(W, env)
                        #time.sleep(0.05)
                traj_new = self.tau


                return [traj_old, traj_new]

                #W = self.TrajectoryToWaypoints(self.tau, N=220)
                #print np.around(np.array(W).T,decimals=2)

