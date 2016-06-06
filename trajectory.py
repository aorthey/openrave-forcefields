import cvxpy as cvx
from numpy import sqrt,sin,cos,pi
from scipy.interpolate import PPoly

import abc
import time
import numpy as np
from util import Rz
from scipy.interpolate import interp1d,splev,splrep,splprep
from scipy.misc import derivative
from pylab import plot,title,xlabel,ylabel,figure
from matplotlib.collections import LineCollection
import pylab as plt

from util_force import *
from util_mvc import *
from util_mvc_approx import *
from topp_interface import TOPPInterface

class Trajectory():
        __metaclass__ = abc.ABCMeta
        DEBUG = 0

        DISCRETIZATION_TIME_STEP = 0.02

        rave_traj = []
        traj = []
        waypoints = []
        handle = []
        #env_ptr = []

        SMOOTH_CONSTANT=0 ##TODO: do not change to >0 => problems with PPoly
        POLYNOMIAL_DEGREE=3
        MIN_NUMBER_WAYPOINTS = 5

        ##drawing parameters
        ptsize = 0.03
        critical_pt_size = 0.07

        show_tangent_vector = False
        show_orientation_vector = True

        lw_path = 10
        lw_tangent = 3
        lw_orientation = 3
        FONT_SIZE = 20
        dVECTOR_LENGTH = 0.5
        trajectory_orientation_color = np.array((0.9,0.9,0.9,0.3))
        trajectory_tangent_color = np.array((0.9,0.2,0.9,0.3))
        trajectory_color_deformation = np.array((0.9,0.2,0.9,0.9))

        tangent_zoffset = -0.02

        trajectory_color_feasible = np.array((0.2,0.9,0.2,0.9))
        trajectory_color_infeasible = np.array((0.9,0.2,0.2,0.9))

        trajectory_color = np.array((0.9,0.2,0.2,0.9))
        critical_pt_color = np.array((0.9,0.2,0.2,0.9))

        bspline = []
        trajectorystring = []
        Ndim = 0

        ### waypoints: real matrix of dimensionality Ndim x Nwaypoints
        def __init__(self, waypoints_in):
                self.new_from_waypoints(waypoints_in)

        def new_from_waypoints(self, waypoints_in):
                self.waypoints = self.prettify(waypoints_in)
                self.Ndim = self.waypoints.shape[0]
                [B,T,D] = self.computeTrajectoryStringForTOPP(self.waypoints)
                self.bspline = B
                self.trajectorystring = T
                self.durationVector = D

        def evaluate_at(self, t, der=1):
                f = np.zeros((self.Ndim))
                df = np.zeros((self.Ndim))
                if der>1:
                        ddf = np.zeros((self.Ndim))
                for i in range(0,self.Ndim):
                        f[i] = splev(t,self.bspline[i])
                        df[i] = splev(t,self.bspline[i],der=1)
                        if der>1:
                                ddf[i] = splev(t,self.bspline[i],der=2)

                df[2] = 0
                if der>1:
                        ddf[2] = 0
                        return [f,df,ddf]
                else:
                        return [f,df]

        #def ForwardSimulate(self, dt = 0.01):

        def prettify(self, W):
                if W.ndim <= 1:
                        print "cannot create trajectory with only one waypoint"
                        sys.exit(1)
                W = self.removeDuplicateWaypoints(W)
                Nwaypoints = W.shape[1]
                if Nwaypoints <= 3:
                        W = self.addMinimalWaypoints(W)
                        Nwaypoints = W.shape[1]
                return W
        def removeDuplicateWaypoints(self, W):
                ACCURACY_DUPLICATE = 1e-5
                Nwaypoints = W.shape[1]
                if Nwaypoints == 1:
                        print "cannot create trajectory with only one waypoint"
                        sys.exit(1)
                if Nwaypoints == 2:
                        if np.linalg.norm(W[:,0]-W[:,1])<ACCURACY_DUPLICATE:
                                print "only two identical waypoints. abort"
                                sys.exit(1)
                        return W

                assert(Nwaypoints>2)
                if np.linalg.norm(W[:,0]-W[:,1])<ACCURACY_DUPLICATE:
                        W = np.delete(W,0,axis=1)
                if np.linalg.norm(W[:,-2]-W[:,-1])<ACCURACY_DUPLICATE:
                        W = np.delete(W,-1,axis=1)

                return W

        def addMinimalWaypoints(self, W):
                Nwaypoints = W.shape[1]
                if Nwaypoints >= self.MIN_NUMBER_WAYPOINTS:
                        return W
                Wtmp = W
                for i in range(0,Nwaypoints-1):
                        Wnew = W[:,i]+0.5*(W[:,i+1]-W[:,i])
                        Wtmp = np.insert(Wtmp, 2*i+1, values=Wnew, axis=1)
                W = Wtmp
                return self.addMinimalWaypoints(W)

        ### SYSTEM DYNAMICS
        def getControlMatrix(self, W):
                AM = 5

                ### car/sailboat
                amin = np.array((-AM,-AM,-3))
                amax = np.array((AM,AM,3))

                ### bacteriophage
                #amin = np.array((0,0,-AM))
                #amax = np.array((AM,0,AM))

                [Ndim,Nwaypoints] = self.getWaypointDim(W)
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

                return [R,amin,amax]


        def computeTrajectorySubstringForTOPP(self, Win, Nc):
                W = Win[:,0:Nc]
                return self.computeTrajectoryStringForTOPP(W, DEBUG=1)

        def computeTrajectoryStringForTOPP(self, W, DEBUG=0):
                Ndim = W.shape[0]
                Nwaypoints = W.shape[1]

                ###############################################################
                ##### get number of intervals between breakpoints
                tvec = np.linspace(0,1,W.shape[1])
                trajectory = splrep(tvec,W[0,:],k=self.POLYNOMIAL_DEGREE,s=self.SMOOTH_CONSTANT)
                poly= PPoly.from_spline(trajectory)
                [B,idx]=np.unique(poly.x,return_index=True)
                Ninterval = B.shape[0]-1
                ###############################################################
                P = np.zeros((Ninterval, Ndim, 4))

                durationVector = np.zeros((Ninterval))
                for j in range(1,B.shape[0]):
                        d = B[j]-B[j-1]
                        durationVector[j-1] = d
                bspline = []
                for i in range(0,Ndim):

                        tvec = np.linspace(0,1,Nwaypoints)
                        trajectory = splrep(tvec,W[i,:],k=self.POLYNOMIAL_DEGREE,s=self.SMOOTH_CONSTANT)
                        bspline.append(trajectory)

                        poly= PPoly.from_spline(trajectory)
                        dpoly = poly.derivative(1)
                        ddpoly = poly.derivative(2)
                        [B,idx]=np.unique(poly.x,return_index=True)
                        coeff = poly.c[:,idx]

                        try:
                                P[:,i,0] = coeff[3,:-1]
                                P[:,i,1] = coeff[2,:-1]
                                P[:,i,2] = coeff[1,:-1]
                                P[:,i,3] = coeff[0,:-1]
                        except Exception as e:
                                print Ninterval
                                print "TOPP EXCEPTION: ",e
                                print B
                                print coeff
                                sys.exit(0)

                        #### TODO: remove z-coordinates for now
                        if i == 2:
                                P[:,i,3]=0
                                P[:,i,2]=0
                                P[:,i,1]=0

                for i in range(0,durationVector.shape[0]):
                        duration = durationVector[i]
                        if i==0:
                                trajectorystring = str(duration)
                        else:
                                trajectorystring += "\n" + str(duration)
                        trajectorystring += "\n" + str(Ndim)

                        for j in range(Ndim):
                                trajectorystring += "\n"
                                trajectorystring += string.join([str(P[i,j,0]),str(P[i,j,1]),str(P[i,j,2]),str(P[i,j,3])])

                return [bspline, trajectorystring, durationVector]


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


        def plot_reachableset(self, env):
                thetavec = [-pi/2,-pi/4,0,pi/4,pi/2]
                #thetavec = [0]
                for theta in thetavec:
                        p = np.array((0,0,0,theta))

                        #dp = np.array((1,0.1,0,0))
                        force = np.array((0,-6.5,0,0))
                        dp = np.array((1,0.1,0,0))
                        #force = np.array((0,-2.5,0,0))
                        smax = 0.1
                        s = 0.1

                        [R,amin,amax] = self.getControlMatrix(p)

                        from reachable_set import ReachableSet
                        self.reach = ReachableSet( p, s, dp, force, R[:,:,0], amin, amax)
                        self.reach.Plot()


        def get_dimension(self):
                [F,dF] = self.evaluate_at(0)
                return F.shape[0]

        def waypoint_to_force(self, env, W):
                Ndims = W.shape[0]
                pt = np.array(((W[0],W[1],-0.1,0.001)))
                F = np.zeros((Ndims))
                F[0:3] = env.GetForceAtX(pt)
                r = 0.5

                theta = W[3]
                rcom = r * np.dot(Rz(theta),ex)
                torque = np.cross(rcom,F[0:2])

                F[3] = np.sign(torque[2])*np.linalg.norm(torque)
                return F
        
        def IsInCollision(self, env, W):
                [Ndim,Nwaypoints] = self.getWaypointDim(W)

                for i in range(0,Nwaypoints):
                        if env.CheckCollisionAtX(W[:,i]):
                                return True
                return False

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

        def getWaypointDim(self, W):
                Ndim = W.shape[0]
                if W.ndim>1:
                        Nwaypoints = W.shape[1]
                else:
                        Nwaypoints = 1
                return [Ndim, Nwaypoints]

        def IsReparametrizable(self, env):
                S = self.reparametrize(env, ploting=False)
                if S is None:
                        return False
                else:
                        return True

        def PlotParametrization(self, env):
                [W,dW,ddW] = self.get_waypoints_second_order()

                [Ndim, Nwaypoints] = self.getWaypointDim(W)
                F = self.get_forces_at_waypoints(W, env)
                [R,amin,amax] = self.getControlMatrix(W)

                self.topp = TOPPInterface(self, self.durationVector, self.trajectorystring, -F,R,amin,amax,W,dW)
                if self.topp.ReparameterizeTrajectory():
                        self.topp.PlotTrajectory(env)
                        print self.topp
                else:
                        print "Trajectory has no ReParameterization"

        def getVelocityIntervalWithoutForceField(self, env, W, dW, ddW):
                [Ndim, Nwaypoints] = self.getWaypointDim(W)
                Fzero = np.zeros((Ndim, Nwaypoints))
                [R,amin,amax] = self.getControlMatrix(W)

                self.topp = TOPPInterface(self, self.durationVector, self.trajectorystring, -Fzero,R,amin,amax,W,dW)
                if self.topp.ReparameterizeTrajectory():
                        return self.topp.traj0
                else:
                        print "WARNING: without force field, TOPP couldn't find a valid \
                        velocity profile. The system is likely not STLC"
                        return None
                        #sys.exit(1)


        def getCriticalPointFromWaypoints(self, env, W, dW, ddW, oldNc = 0):

                [Ndim, Nwaypoints] = self.getWaypointDim(W)
                F = self.get_forces_at_waypoints(W, env)
                [R,amin,amax] = self.getControlMatrix(W)

                self.topp = TOPPInterface(self, self.durationVector, self.trajectorystring, -F,R,amin,amax,W,dW)
                Nc = self.topp.getCriticalPoint()

                if Nc < 0:
                        print "return oldNc=",oldNc
                        Nc = oldNc

                ### AVP on the subtrajectory between 0 and Nc
                #if Nc > 0:
                #        [bspline, trajstr, durationVector] = self.computeTrajectorySubstringForTOPP(W, Nc)
                #        [semin, semax] = self.topp.getSpeedIntervalAtCriticalPoint(Nc, durationVector, trajstr)
                #else:
                #        semin = 0.0
                #        semax = 0.0
                #print "MAX SPEED AT CRITICAL PT",Nc,"/",Nwaypoints," is:",semin,semax

                return Nc

        def getCriticalPoint(self, env):
                L = self.get_length()
                ds = 0.01
                [W,dW,ddW] = self.get_waypoints_second_order()

                N = self.getCriticalPointFromWaypoints(env, W, dW, ddW)
                if N is not None:
                        print "CRITICAL POINT:",N,"/",W.shape[1]
                        print W[:,N]
                else:
                        print "No critical point found"

        def reparametrize(self, env, ploting=False):
                L = self.get_length()
                [W,dW,ddW] = self.get_waypoints_second_order()
                F = self.get_forces_at_waypoints(W, env)

                ### acceleration limits
                [R,amin,amax] = self.getControlMatrix(W)

                S = getSpeedProfileRManifold(-F,R,amin,amax,W,dW,ddW,ploting)
                #if S is None:
                        #self.speed_profile = S
                        #print "No parametrization solution -- sorry :-((("
                        #return False
                return S

        @classmethod
        def from_ravetraj(cls, ravetraj):
                N = ravetraj.GetNumWaypoints()
                W=[]
                for i in range(0,N):
                        w = np.array((ravetraj.GetWaypoint(i)[0],ravetraj.GetWaypoint(i)[1],ravetraj.GetWaypoint(i)[2],ravetraj.GetWaypoint(i)[3]))
                        W.append((w))
                W = np.array(W).T
                return cls(W)

        def get_waypoints_second_order(self, N=None):
                ###############################################################
                ### obtain waypoints along the trajectory at constant spacing of
                ### DISCRETIZATION_TIME_STEP
                ###############################################################
                if N is None:
                        dt = self.DISCRETIZATION_TIME_STEP
                        L = self.get_length()
                        N = int(L/dt)
                        #print N,"WAYPOINTS <<<<<<<<<<<"
                ###############################################################
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

        def get_waypoints(self, N = None):
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
                if Nwaypoints > 0:
                        with env.env:
                                Wext = np.zeros((3, 2*Nwaypoints))
                                Wtheta = np.zeros((3, 3*Nwaypoints))

                                for i in range(0,Nwaypoints):
                                        pt = np.array(((W[0,i],W[1,i],W[2,i]+self.tangent_zoffset)))
                                        dpt = np.array(((dW[0,i],dW[1,i],dW[2,i]+self.tangent_zoffset)))
                                        dpt = self.dVECTOR_LENGTH*dpt/np.linalg.norm(dpt)

                                        Wext[:,2*i] = pt
                                        Wext[:,2*i+1] = np.array( (pt[0]+dpt[0],pt[1]+dpt[1],pt[2]+dpt[2]) )

                                        #### orientation of system
                                        theta = W[3,i]
                                        etheta = self.dVECTOR_LENGTH*np.dot(Rz(theta),ex)
                                        etheta[2] += self.tangent_zoffset
                                        Wtheta[:,3*i] = pt
                                        Wtheta[:,3*i+1] = np.array( (pt[0]+etheta[0],pt[1]+etheta[1],pt[2]+etheta[2]) )
                                        Wtheta[:,3*i+2] = pt

                                #Wext = np.array(zip(W[0:3,:],Wdir)).flatten().reshape(Ndim,2*Nwaypoints)
                                if self.show_tangent_vector:
                                        h=env.env.drawlinestrip(points=Wext.T,linewidth=self.lw_tangent,colors=self.trajectory_tangent_color)
                                        tmp_handle.append(h)

                                if self.show_orientation_vector:
                                        h=env.env.drawlinestrip(points=Wtheta.T,linewidth=self.lw_orientation,colors=self.trajectory_orientation_color)
                                        tmp_handle.append(h)

                                h=env.env.drawlinestrip(points=W[0:3,:].T,linewidth=self.lw_path,colors=self.trajectory_color)
                                tmp_handle.append(h)


                return tmp_handle

        def draw(self, env, keep_handle=True, critical_pt = None):
                [W,dW,ddW] = self.get_waypoints_second_order()
                t1 = time.time()
                if critical_pt == None:
                        N = self.getCriticalPointFromWaypoints(env, W, dW, ddW)
                else:
                        N = critical_pt

                t2 = time.time()
                tmp_handle = []

                self.trajectory_color = self.trajectory_color_feasible
                tmp_handle.append(self.get_handle_draw_waypoints(env, W[:,0:N], dW[:,0:N], ddW[:,0:N]))
                self.trajectory_color = self.trajectory_color_infeasible
                tmp_handle.append(self.get_handle_draw_waypoints(env, W[:,N:], dW[:,N:], ddW[:,N:]))
                t3 = time.time()
                if self.DEBUG:
                        print "draw loop: ",t3-t1," topp:",t2-t1," draw:",t3-t2

                if keep_handle:
                        self.handle = tmp_handle
                else:
                        return tmp_handle

        def draw_delete(self):
                self.handle = []
