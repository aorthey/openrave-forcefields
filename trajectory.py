import cvxpy as cvx
from cvxpy import Variable, Problem, Minimize, norm, SCS, OPTIMAL
from cvxpy import norm
from numpy import sqrt,sin,cos,pi

import abc
import time
import numpy as np
from util import Rz
from scipy.interpolate import interp1d,splev,splrep,splprep
from scipy.misc import derivative
from pylab import plot,title,xlabel,ylabel,figure
import pylab as plt

from util_force import *
from util_mvc import *
from util_mvc_approx import *

class Trajectory():
        __metaclass__ = abc.ABCMeta

        rave_traj = []
        traj = []
        waypoints = []
        handle = []

        ##drawing parameters
        ptsize = 0.03
        linsize = 1.5
        FONT_SIZE = 20
        dVECTOR_LENGTH = 0.5
        trajectory_color = np.array((0.7,0.2,0.7,0.9))

        def __init__(self, trajectory, waypoints_in = []):
                self.traj = trajectory 
                self.waypoints = waypoints_in

        def info(self):
                print "#### TRAJECTORY CLASS ######"
                print "LENGTH: ", self.get_length()

                [f0,df0] = np.around(self.evaluate_at(0),2)
                [f1,df1] = np.around(self.evaluate_at(1),2)
                print "START      : ", f0, " dir: ", df0
                print "GOAL       : ", f1, " dir: ", df1
                print "WAYPOINTS  : "
                print "             [   T   ]  [  WAYPOINT  ]"
                for tt in np.linspace(0,1,10):
                        [f,df]=self.evaluate_at(tt)
                        print "             [ ",np.around(tt,1)," ] ",np.around(f,decimals=2)

        def get_dimension(self):
                [F,dF] = self.evaluate_at(0)
                return F.shape[0]

        def waypoint_to_force(self, env, W):
                Ndims = W.shape[0]
                pt = np.array(((W[0],W[1],-0.1,0.001)))
                F = np.zeros((Ndims))
                F[0:3] = env.GetForceAtX(pt)
                return F

        def get_forces_at_waypoints(self, W, env):
                Ndim = W.shape[0]
                if W.ndim>1:
                        Nwaypoints = W.shape[1]
                else:
                        Nwaypoints = 1
                F = np.zeros((Ndim,Nwaypoints))
                if Nwaypoints>1:
                        for i in range(0,Nwaypoints):
                                F[:,i] = self.waypoint_to_force(env, W[:,i])
                else:
                        F[:,0] = self.waypoint_to_force(env,W)
                return F

        def normalize_velocity_vector(self, dW):
                Ndim = dW.shape[1]
                for i in range(0,Ndim):
                        dW[:,i] = dW[:,i]/np.linalg.norm(dW[:,i])
                return dW

        def plot_speed_profile(self, env):

                L = self.get_length()
                dt = 0.05
                Nwaypoints = int(L/dt)

                [W,dW] = self.get_waypoints(N=Nwaypoints)
                F = self.get_forces_at_waypoints(W, env)
                dW = self.normalize_velocity_vector(dW)
                ### we need F, W, dW(normalized)

                #V = np.zeros((Nwaypoints-1))
                #for i in range(0,Nwaypoints-1):
                        #V[i] = self.get_minimum_feasible_velocity(dt, W[:,i],W[:,i+1],dW[:,i],F[:,i])

                epsilon = 0.01
                T = np.linspace(0.0,1.0,num=Nwaypoints-1)
                V = self.get_minimum_feasible_velocity_vector(dt,epsilon, W,dW,F)

                plot(T,V,'or',linewidth=3,markersize=2)        
                title('Speed Profile Trajectory', fontsize=self.FONT_SIZE)
                xlabel('$\theta$', fontsize=self.FONT_SIZE)
                ylabel('$s$', fontsize=self.FONT_SIZE)
                plt.show()

        def plot_speed_profile(self, P):
                Nwaypoints = P.shape[0]
                S = np.linspace(0,1,Nwaypoints)

                fig=figure(facecolor='white')
                plot(S,P,'-r',linewidth=5,markersize=8)        
                plt.tick_params(labelsize=self.FONT_SIZE)
                title('Speed Profile Trajectory', fontsize=self.FONT_SIZE)
                xlabel('$s$', fontsize=2*self.FONT_SIZE)
                ylabel('$\dot s$', fontsize=2*self.FONT_SIZE)

                plt.show()


        def getControlMatrix(self, W):
                Ndim = W.shape[0]
                if W.ndim>1:
                        Nwaypoints = W.shape[1]
                else:
                        Nwaypoints = 1
                assert(Ndim==4)
                Kdim = 3
                R = np.zeros((Ndim,Kdim,Nwaypoints))
                for i in range(0,Nwaypoints):
                        if Nwaypoints>1:
                                t = W[3,i]
                        else:
                                t = W[3]
                        R[0,:,i] = np.array((cos(t),-sin(t),0.0))
                        R[1,:,i] = np.array((sin(t),cos(t),0.0))
                        R[2,:,i] = np.array((0.0,0.0,0.0))
                        R[3,:,i] = np.array((0.0,0.0,1.0))

                amin = np.array((-5,-5,-5))
                amax = np.array((5,5,5))
                return [R,amin,amax]

        def getWaypointDim(self, W):
                Ndim = W.shape[0]
                if W.ndim>1:
                        Nwaypoints = W.shape[1]
                else:
                        Nwaypoints = 1
                return [Ndim, Nwaypoints]

        def computeReachableSets(self, dt, env):
                L = self.get_length()
                Nwaypoints = int(L/0.01)
                [W,dW,ddW] = self.get_waypoints_second_order(N=Nwaypoints)

                [Ndim, Nwaypoints] = self.getWaypointDim(W)

                F = self.get_forces_at_waypoints(W, env)
                [R,amin,amax] = self.getControlMatrix(W)

                Rmin = np.zeros((Nwaypoints, Ndim))
                Rmax = np.zeros((Nwaypoints, Ndim))
                for i in range(0,Nwaypoints):
                        Rmin[i,:] = np.minimum(np.dot(R[:,:,i],amax),np.dot(R[:,:,i],amin))
                        Rmax[i,:] = np.maximum(np.dot(R[:,:,i],amax),np.dot(R[:,:,i],amin))

                dt2 = dt*dt*0.5
                Qoff = W[:,i] + dt*dW[:,i] + dt2*F[:,i]

                Lmin = Qoff + Rmin
                delta = np.zeros((Nwaypoints))
                #xvar = Variable(Ndim,Nwaypoints)
                #yvar = Variable(Ndim,Nwaypoints)
                constraints = []
                objfunc = 0.0
                for i in range(0,Nwaypoints-1):
                        xvar = Variable(Ndim)
                        yvar = Variable(Ndim)
                        #constraints.append( Qoff + Rmin[i,:] <= xvar )
                        #constraints.append( xvar <= Qoff + Rmax[i,:] ) ## Q
                        #constraints.append( W[:,i] + eta[i]*dW[:,i] == y[:,i] ) ## L
                        #constraints.append( yvar >= 0)
                        objfunc += norm(xvar)

                objective = Minimize(objfunc)
                prob = Problem(objective, constraints)

                print "solve minimal speed"
                result = prob.solve(solver=SCS)
                for i in range(0,Nwaypoints):
                        delta[i] = np.linalg.norm(x[:,i]-y[:,i])
                print delta
                plot(delta,'-r')
                plt.show()


                
        def computeDeltaProfile(self):
                pass

        def reparametrize(self, env):
                L = self.get_length()
                dt = 0.05
                Nwaypoints = int(L/dt)
                [W,dW,ddW] = self.get_waypoints_second_order(N=Nwaypoints)
                F = self.get_forces_at_waypoints(W, env)

                ### acceleration limits
                [R,amin,amax] = self.getControlMatrix(W)

                print F
                S = getSpeedProfileRManifold(F,R,amin,amax,W,dW,ddW)
                if S is None:
                        self.speed_profile = S
                        print "No parametrization solution -- sorry :-((("
                        return False
                return True

        def speed_profile_MVC(self, ravetraj):
                poly_traj = TOPPopenravepy.FromRaveTraj(self.robot, ravetraj)

                grav=[0,0,-9.8]
                ndof=self.robot.GetDOF()
                dof_lim=self.robot.GetDOFLimits()
                vel_lim=self.robot.GetDOFVelocityLimits()
                #robot.SetDOFLimits(-10*ones(ndof),10*ones(ndof)) # Overrides robot joint limits for TOPP computations
                #robot.SetDOFVelocityLimits(100*vel_lim) # Override robot velocity limits for TOPP computations
                
                ret = traj.RunComputeProfiles(0,0)
                if(ret == 1):
                        traj.ReparameterizeTrajectory()
                        return traj
                else:
                        print "Trajectory is not time-parameterizable"

        from util import Rz

        ### return the disturbance force after we counteract them with the best
        ### acceleration available (from our dt reachable set)
        def get_minimal_disturbance_forces(self, dt, W, F, u1min, u1max, u2min, u2max):
                dt2 = dt*dt/2
                Ndim = W.shape[0]
                Nsamples = W.shape[1]
                Fp = np.zeros((Ndim,Nsamples))

                for i in range(0,Nsamples):
                        theta = W[3,i]
                        Fp[0:3,i] = np.dot(Rz(theta),F[0:3,i])[0:3].flatten()

                ### question (1) : can we produce an acceleration ddW, such that it counteracts F?
                nearest_ddq = ForceCanBeCounteractedByAccelerationVector(dt, Fp, u1min, u1max, u2min, u2max)

                Ndim = Fp.shape[0]
                Nsamples = Fp.shape[1]
                dF = np.zeros((Ndim,Nsamples))

                for i in range(0,Nsamples):
                        for j in range(0,Ndim):
                                dF[j,i] = (Fp[j,i] - nearest_ddq[j,i])

                for i in range(0,Nsamples):
                        theta = W[3,i]
                        dF[0:3,i] = np.dot(Rz(theta).T,dF[0:3,i])[0:3].flatten()
                return dF

        def get_minimum_feasible_velocity_vector(self, dt, epsilon, W, dW, F):

                u1min = 0.0
                u1max = 2.0
                u2min = -2.0
                u2max = 2.0
                dt2 = dt*dt/2

                ## assume that W contains [x,y,z,theta]
                ## note: a force along x,y,z corresponds to a linear acc field, a force along theta would be a twisted acc field 

                ## project F onto the identity element, here the identity is at (0,0,0,0), with x pointing towards the 'right'
                Ndim = W.shape[0]
                Nsamples = W.shape[1]
                Fp = np.zeros((2,Nsamples))
                for i in range(0,Nsamples):
                        theta = W[3,i]
                        Fp[:,i] = np.dot(Rz(theta).T,F[:,i])[0:2].flatten()

                ### question (1) : can we produce an acceleration ddW, such that it counteracts F?
                nearest_ddq = ForceCanBeCounteractedByAccelerationVector(dt, Fp, u1min, u1max, u2min, u2max)

                ### question (2) : if ddW cannot counteract F, what is the minimal speed dW, such that we still follow tau?
                ## check that we are in an epsilon neighborhood of next waypoint!
                Ndim = Fp.shape[0]
                Nsamples = Fp.shape[1]
                dF = np.zeros((Ndim,Nsamples))

                for i in range(0,Nsamples):
                        for j in range(0,Ndim):
                                dF[j,i] = dt2*Fp[j,i] - nearest_ddq[j,i]

                vmin = GetMinimalSpeedToReachEpsilonNeighbordhoodVector(dt, epsilon, W, dW, dF)

                return vmin

        def get_minimum_feasible_velocity(self, dt, W0, W1, dW0, F):
                epsilon = 0.01

                ## assume that W contains [x,y,z,theta]
                ## note: a force along x,y,z corresponds to a linear acc field, a force along theta would be a twisted acc field 

                ## project F onto the identity element, here the identity is at (0,0,0,0), with x pointing towards the 'right'
                theta = W0[3]
                R = Rz(theta)
                Fp = np.dot(R.T,F)
                Fp = Fp[0:2].flatten()

                u1min = 0.0
                u1max = 2.0
                u2min = -2.0
                u2max = 2.0

                dt2 = dt*dt/2

                ### question (1) : can we produce an acceleration ddW, such that it counteracts F?
                [d,nearest_ddq] = ForceCanBeCounteractedByAcceleration(dt, Fp, u1min, u1max, u2min, u2max)

                if d <= 0.0001:
                        #print "force can be counteracted by system dynamics"
                        #print "no need for having a fixed velocity"
                        return 0.0
                
                ### question (2) : if ddW cannot counteract F, what is the minimal speed dW, such that we still follow tau?
                ## check that we are in an epsilon neighborhood of next waypoint!
                dF = np.zeros((2,1))
                dF[0] = dt2*Fp[0] - nearest_ddq[0]
                dF[1] = dt2*Fp[1] - nearest_ddq[1]


                vmin = GetMinimalSpeedToReachEpsilonNeighbordhood(dt, epsilon, W0, W1, dW0, dF)

                return vmin

        @abc.abstractmethod 
        def evaluate_at(self, t):
                pass

        @classmethod
        def from_waypoints(cls, W):
                pass

        @classmethod
        def from_ravetraj(cls, ravetraj):
                rave_traj = ravetraj
                N = ravetraj.GetNumWaypoints()
                W=[]
                for i in range(0,N):
                        w = np.array((ravetraj.GetWaypoint(i)[0],ravetraj.GetWaypoint(i)[1],ravetraj.GetWaypoint(i)[2],ravetraj.GetWaypoint(i)[3]))
                        W.append((w))
                W = np.array(W).T
                return cls.from_waypoints(W)

        def get_waypoints_second_order(self, N = 100):
                K = self.get_dimension()
                pts = np.zeros((K,N))
                dpts = np.zeros((K,N))
                ddpts = np.zeros((K,N))
                ctr = 0
                for t in np.linspace(0.0, 1.0, num=N):
                        [f0,df0,ddf0] = self.evaluate_at(t,der=2)
                        pts[:,ctr] = f0
                        dpts[:,ctr] = df0
                        ddpts[:,ctr] = ddf0
                        ctr = ctr+1
                return [pts,dpts,ddpts]

        def get_waypoints(self, N = 100):
                [pts,dpts,ddpts]=self.get_waypoints_second_order(N)
                return [pts,dpts]

        def get_length(self):
                dd = 0.0
                T = np.linspace(0.0, 1.0, num=100)
                for i in range(0,len(T)-1):
                        t = T[i]
                        tdt = T[i+1]
                        [ft,df0] = self.evaluate_at(t)
                        [ftdt,df0] = self.evaluate_at(tdt)
                        dd = dd + np.linalg.norm(ft-ftdt)
                return dd

        def forward_simulate_one_step(self, dt, oldspeed, W, dW, F, A):
                dt2=dt*dt/2.0
                x = W[0]
                y = W[1]
                z = W[2]
                theta = W[3]

                speedinc = A[0]*np.array((cos(theta),sin(theta),0,0)) + A[1]*np.array((0,0,0,1))
                Wcur = W + dW*oldspeed*dt + speedinc*dt2 + F*dt2
                dWcur = oldspeed*dW + F*dt + speedinc*dt
                speed = np.linalg.norm(dWcur)
                dWcur = dWcur/speed

                return [Wcur,dWcur,speed]

        def forward_simulate(self, acc_profile, env):
                Ndim = self.get_dimension()
                Nwaypoints=acc_profile.shape[0]
                W = np.zeros((Ndim,Nwaypoints))
                dW = np.zeros((Ndim,Nwaypoints))

                [f0,df0] = self.evaluate_at(0)
                W[:,0] = f0
                dW[:,0] = df0
                dt = 0.01

                S = np.linspace(0.0,1.0,num=Nwaypoints)
                oldspeed = 0.0

                for i in range(0,Nwaypoints-1):
                        scur = S[i]
                        F = self.waypoint_to_force(env, W[:,i])
                        A = acc_profile[i,:]
                        [W[:,i+1],dW[:,i+1],oldspeed]=self.forward_simulate_one_step(dt, oldspeed, W[:,i], dW[:,i], F, A)
                        #print W[:,i+1],dt,dW[:,i]

                self.handle = self.get_handle_draw_waypoints(env, W, dW)

        def get_handle_draw_waypoints(self, env, W, dW, ddW):
                Ndims = W.shape[0]
                Nwaypoints = W.shape[1]
                tmp_handle = []
                for i in range(0,Nwaypoints):
                        pt = np.array(((W[0,i],W[1,i],W[2,i])))
                        tmp_handle.append(env.env.plot3(points=pt,
                                           pointsize=self.ptsize,
                                           colors=self.trajectory_color,
                                           drawstyle=1))
                        dpt = np.array(((dW[0,i],dW[1,i],dW[2,i])))
                        ddpt = np.array(((ddW[0,i],ddW[1,i],ddW[2,i])))
                        dpt = self.dVECTOR_LENGTH*dpt/np.linalg.norm(dpt)
                        ddpt = self.dVECTOR_LENGTH*ddpt/np.linalg.norm(ddpt)

                        P = np.array(((pt[0],pt[1],pt[2]),
                                (pt[0]+dpt[0],pt[1]+dpt[1],pt[2]+dpt[2])))
                        h=env.env.drawlinestrip(points=P,linewidth=self.linsize,colors=np.array(((0.2,0.9,0.2,0.9))))
                        tmp_handle.append(h)

                        P = np.array(((pt[0],pt[1],pt[2]),
                                (pt[0]+ddpt[0],pt[1]+ddpt[1],pt[2]+ddpt[2])))
                        h=env.env.drawlinestrip(points=P,linewidth=self.linsize,colors=np.array(((0.9,0.9,0.9,0.9))))
                        tmp_handle.append(h)

                return tmp_handle

        def draw(self, env, keep_handle=True):
                Nwaypoints=200
                [W,dW,ddW] = self.get_waypoints_second_order(N=Nwaypoints)
                tmp_handle = self.get_handle_draw_waypoints(env, W, dW, ddW)

                if keep_handle:
                        self.handle = tmp_handle
                else:
                        return tmp_handle


        def draw_delete(self, env):
                self.handle = []
