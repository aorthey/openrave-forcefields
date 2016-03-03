import cvxpy as cvx
from cvxpy import Variable, Problem, Minimize, norm, SCS, OPTIMAL
from cvxpy import SCS, CVXOPT
from numpy import sqrt,sin,cos,pi

import abc
import time
import numpy as np
from util import Rz
from scipy.interpolate import interp1d,splev,splrep,splprep
from scipy.misc import derivative
from pylab import plot,title,xlabel,ylabel
import pylab as plt

from util_force import *

class Trajectory():
        __metaclass__ = abc.ABCMeta

        traj = []
        waypoints = []
        handle = []

        ptsize = 0.03
        linsize = 1.5
        FONT_SIZE = 10
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

        def get_forces_at_waypoints(self, W, env):
                Ndims = W.shape[0]
                Nwaypoints = W.shape[1]
                F = np.zeros((Ndims,Nwaypoints))
                for i in range(0,Nwaypoints):
                        pt = np.array(((W[0,i],W[1,i],-0.1,0.001)))
                        F[0:3,i] = env.GetForceAtX(pt)
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
                        Fp[0:3,i] = np.dot(Rz(theta).T,F[0:3,i])[0:3].flatten()

                ### question (1) : can we produce an acceleration ddW, such that it counteracts F?
                nearest_ddq = ForceCanBeCounteractedByAccelerationVector(dt, Fp, u1min, u1max, u2min, u2max)

                Ndim = Fp.shape[0]
                Nsamples = Fp.shape[1]
                dF = np.zeros((Ndim,Nsamples))

                for i in range(0,Nsamples):
                        for j in range(0,Ndim):
                                dF[j,i] = (Fp[j,i] - nearest_ddq[j,i])
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
                N = ravetraj.GetNumWaypoints()
                W=[]
                for i in range(0,N):
                        w = np.array((ravetraj.GetWaypoint(i)[0],ravetraj.GetWaypoint(i)[1],ravetraj.GetWaypoint(i)[2],ravetraj.GetWaypoint(i)[3]))
                        W.append((w))
                W = np.array(W).T
                return cls.from_waypoints(W)

        def get_waypoints(self, N = 100):
                K = self.get_dimension()
                pts = np.zeros((K,N))
                dpts = np.zeros((K,N))
                ctr = 0
                for t in np.linspace(0.0, 1.0, num=N):
                        [f0,df0] = self.evaluate_at(t)
                        pts[:,ctr] = f0
                        dpts[:,ctr] = df0
                        ctr = ctr+1
                return [pts,dpts]

        def get_length(self):
                dt = 0.05
                dd = 0.0
                for t in np.linspace(0.0, 1.0-dt, num=1000):
                        [ft,df0] = self.evaluate_at(t)
                        [ftdt,df0] = self.evaluate_at(t+dt)
                        dd = dd + np.linalg.norm(ft-ftdt)
                return dd

        def draw(self, env, keep_handle=True):
                Nwaypoints=200
                [W,dW] = self.get_waypoints(N=Nwaypoints)

                tmp_handle = []
                for i in range(0,Nwaypoints):
                        pt = np.array(((W[0,i],W[1,i],W[2,i])))
                        tmp_handle.append(env.env.plot3(points=pt,
                                           pointsize=self.ptsize,
                                           colors=self.trajectory_color,
                                           drawstyle=1))
                        dpt = np.array(((dW[0,i],dW[1,i],dW[2,i])))
                        dpt = self.dVECTOR_LENGTH*dpt/np.linalg.norm(dpt)

                        P = np.array(((pt[0],pt[1],pt[2]),
                                (pt[0]+dpt[0],pt[1]+dpt[1],pt[2]+dpt[2])))
                        h=env.env.drawlinestrip(points=P,linewidth=self.linsize,colors=np.array(((0.2,0.9,0.2,0.9))))
                        tmp_handle.append(h)


                if keep_handle:
                        self.handle = tmp_handle
                else:
                        return tmp_handle


        def draw_delete(self, env):
                self.handle = []
