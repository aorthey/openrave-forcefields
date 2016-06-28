import cvxpy as cvx
from numpy import sqrt,sin,cos,pi
from scipy.interpolate import PPoly

import abc
import time
import numpy as np
from util import Rz, PrintNumpy
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

        DISCRETIZATION_TIME_STEP = 0.01
        SMOOTH_CONSTANT = 0 ##TODO: do not change to >0 => problems with PPoly
        POLYNOMIAL_DEGREE = 1
        MIN_NUMBER_WAYPOINTS = 5

        rave_traj = []
        traj = []
        waypoints = []
        handle = []
        #env_ptr = []

        ##drawing parameters
        ptsize = 0.03
        critical_pt_size = 0.08

        show_tangent_vector = True
        show_orientation_vector = True

        lw_path = 10
        lw_tangent = 3
        lw_orientation = 3
        FONT_SIZE = 20
        dVECTOR_LENGTH = 0.4
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
                self.bspline = self.computeSplineFromWaypoints(self.waypoints)
                [self.waypoints,dW,ddW] = self.get_waypoints_second_order()
                [T,D] = self.computeTrajectoryStringForTOPP(self.waypoints,dW)
                self.trajectorystring = T
                self.durationVector = D

        def save(self, filename):
                np.save(filename, self.waypoints)

        def load(self, filename):
                W = np.load(filename+'.npy')
                return W

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
                ACCURACY_DUPLICATE = self.DISCRETIZATION_TIME_STEP/1e3
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
                        print "deleting waypoint"
                        W = np.delete(W,0,axis=1)
                if np.linalg.norm(W[:,-2]-W[:,-1])<ACCURACY_DUPLICATE:
                        print "deleting waypoint"
                        W = np.delete(W,-1,axis=1)

                return W

        def addMinimalWaypoints(self, W):
                Nwaypoints = W.shape[1]
                if Nwaypoints >= self.MIN_NUMBER_WAYPOINTS:
                        return W
                Wtmp = W
                for i in range(0,Nwaypoints-1):
                        print "adding waypoint"
                        Wnew = W[:,i]+0.5*(W[:,i+1]-W[:,i])
                        Wtmp = np.insert(Wtmp, 2*i+1, values=Wnew, axis=1)
                W = Wtmp
                return self.addMinimalWaypoints(W)

        ### SYSTEM DYNAMICS
        def getControlMatrix(self, W):
                import parameters_dynamical_system as params
                #AM = 1

                amin = params.amin
                amax = params.amax
                [Ndim,Nwaypoints] = self.getWaypointDim(W)
                R = params.ControlPerWaypoint(W, Ndim, Nwaypoints)
                return [R,amin,amax]

        def computeTrajectoryStringForTOPP_deprecated(self, W, DEBUG=0):
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
                ###############################################################

                Kcoeff = self.POLYNOMIAL_DEGREE+1
                P = np.zeros((Ninterval, Ndim, Kcoeff))

                durationVector = np.zeros((Ninterval))
                for j in range(1,B.shape[0]):
                        d = B[j]-B[j-1]
                        durationVector[j-1] = d
                bspline = []

                for i in range(0,Ndim):


                        tvec = np.linspace(0,1,Nwaypoints)
                        WP = W[i,:]

                        #tvec = np.hstack((-0.5,tvec,1.5))
                        #WP = np.hstack((W[i,0],W[i,:],W[i,-1]))

                        trajectory = splrep(tvec,WP,k=self.POLYNOMIAL_DEGREE,s=self.SMOOTH_CONSTANT)
                        bspline.append(trajectory)

                        poly= PPoly.from_spline(trajectory)
                        dpoly = poly.derivative(1)
                        ddpoly = poly.derivative(2)
                        [B,idx]=np.unique(poly.x,return_index=True)
                        coeff = poly.c[:,idx]

                        #print poly.x
                        #print W[i,:]
                        #print splev(0.0,trajectory)
                        #print splev(1.0,trajectory)
                        #print splev(0.0,trajectory,der=1)
                        #print splev(1.0,trajectory,der=1)
                        #print coeff
                        #sys.exit(0)

                        try:
                                for j in range(0,Kcoeff):
                                        P[:,i,j] = coeff[Kcoeff-j-1,:-1]
                                #P[:,i,0] = coeff[3,:-1]
                                #P[:,i,1] = coeff[2,:-1]
                                #P[:,i,2] = coeff[1,:-1]
                                #P[:,i,3] = coeff[0,:-1]
                        except Exception as e:
                                print "TOPP EXCEPTION: ",e
                                print "Kcoeff:",Kcoeff,"Ninterval:",Ninterval
                                print "P.shape",P.shape
                                print i,"/",Ndim
                                print B
                                print coeff
                                sys.exit(0)

                        #### TODO: remove z-coordinates for now
                        if i == 2:
                                for j in range(1,Kcoeff):
                                        P[:,i,j]=0
                                #P[:,i,3]=0
                                #P[:,i,2]=0
                                #P[:,i,1]=0

                for i in range(0,durationVector.shape[0]):
                        duration = durationVector[i]
                        if i==0:
                                trajectorystring = str(duration)
                        else:
                                trajectorystring += "\n" + str(duration)
                        trajectorystring += "\n" + str(Ndim)

                        for j in range(Ndim):
                                trajectorystring += "\n"
                                trajectorystring += string.join(map(str,P[i,j,:]))
                                #trajectorystring += string.join([str(P[i,j,0]),str(P[i,j,1]),str(P[i,j,2]),str(P[i,j,3])])
                                #trajectorystring += string.join([str(P[i,j,0]),str(P[i,j,1])])

                #print trajectorystring
                #sys.exit(0)

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
                thetavec = [0]
                for theta in thetavec:
                        p = np.array((0,1e-4,0,theta))

                        ##sailboat
                        force = np.array((2.5,-3.5,0,5.0))
                        dp = np.array((1,0.0,0,0.2))
                        ##car
                        force = np.array((0.5,-2.5,0,5.0))
                        dp = np.array((1,0.1,0,0.2))
                        #force = np.array((0,0,0,0))
                        s = 0.2

                        #p =np.array( [-2.608897298112661, -0.2047366548132623, 0.1000000008651387, -2.987014267911084] )
                        #dp =np.array( [-16.487757154509474, -4.837209926510647, 0.0, 4.899656054155726] )
                        #force =np.array( [0.0, 3.0, 0.0, -1.4821147977289577] )
                        #s= 3.86342505425

                        p =np.array( [-2.644855145592143, 1.1657167012260294, 0.10000000000235312, -3.7242126663188597] )
                        dp =np.array( [-1.1866896124736256e-20, 1.910882585951823e-18, 1.4346537045461183e-27, -9.702129945467387e-19] )
                        force =np.array( [0.0, 0.0, 0.0, 0.0] )
                        ds= 0.01
                        s= 0.341259539422



                        [R,amin,amax] = self.getControlMatrix(p)

                        #from reachable_set import ReachableSet
                        #self.reach = ReachableSet( p, s, dp, force, R[:,:,0], amin, amax)
                        #self.reach.Plot()

                        from reachable_set3d import ReachableSet3D
                        self.reach = ReachableSet3D( self.DISCRETIZATION_TIME_STEP, p, s, dp, force, R[:,:,0], amin, amax)
                        self.reach.Plot()
                        self.reach.PlotShow()
                        self.reach.PlotSave("images/reachable_set_projection2.png")


        def get_dimension(self):
                [F,dF] = self.evaluate_at(0)
                return F.shape[0]

        def waypoint_to_force(self, env, W):
                Ndims = W.shape[0]
                pt = np.array(((W[0],W[1],-0.1,W[3])))
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

        def GetFirstCollisionPointIdx(self, env, W):
                [Ndim,Nwaypoints] = self.getWaypointDim(W)

                for i in range(0,Nwaypoints):
                        if env.CheckCollisionAtX(W[:,i]):
                                return i
                return None

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
                #Fzero = np.zeros((Ndim, Nwaypoints))
                Fzero = 0.0*self.get_forces_at_waypoints(W, env)
                [R,amin,amax] = self.getControlMatrix(W)

                self.topp = TOPPInterface(self, self.durationVector, self.trajectorystring, Fzero,R,amin,amax,W,dW)
                if self.topp.ReparameterizeTrajectory():
                        return self.topp.traj0
                else:
                        print "WARNING: without force field, TOPP couldn't find a valid \
                        velocity profile. Path not continuous or system not STLC"
                        print W.shape
                        print Fzero
                        print W
                        print self.trajectorystring
                        self.info()
                        sys.exit(1)
                        return None


        def getCriticalPointFromWaypoints(self, env, W, dW, ddW, oldNc = 0):

                [Ndim, Nwaypoints] = self.getWaypointDim(W)
                F = self.get_forces_at_waypoints(W, env)
                [R,amin,amax] = self.getControlMatrix(W)

                self.topp = TOPPInterface(self, self.durationVector, self.trajectorystring, -F,R,amin,amax,W,dW)
                Nc = self.topp.getCriticalPoint()-1

                if Nc < 0:
                        print "return oldNc=",oldNc
                        Nc = oldNc
                return Nc

        def GetSpeedIntervalAtCriticalPoint(self, env, Win, dWin, Nc_in):
#
#                for Nc in range(3,Nc_in):
#                        W = Win[:,0:Nc]
#                        dW = dWin[:,0:Nc]
#
#                        [Ndim, Nwaypoints] = self.getWaypointDim(W)
#                        F = self.get_forces_at_waypoints(W, env)
#                        [R,amin,amax] = self.getControlMatrix(W)
#
#                        [trajsubstr, durationVector] = self.computeTrajectorySubstringForTOPP(Win, dWin, Nc)
#
#                        self.topp = TOPPInterface(self, durationVector, trajsubstr, -F,R,amin,amax,W,dW)
#
#                        #print durationVector,trajsubstr
#                        #print W,dW
#                        [semin,semax] = self.topp.getSpeedIntervalAtCriticalPoint(
#                                        Nc,
#                                        durationVector,
#                                        trajsubstr)
#                        print "TOPP velocity at CP",Nc,":",semin,semax
#                sys.exit(0)

                W = Win[:,0:Nc_in]
                dW = dWin[:,0:Nc_in]

                [Ndim, Nwaypoints] = self.getWaypointDim(W)
                F = self.get_forces_at_waypoints(W, env)
                [R,amin,amax] = self.getControlMatrix(W)

                [trajsubstr, durationVector] = self.computeTrajectorySubstringForTOPP(Win, dWin, Nc_in)

                self.topp = TOPPInterface(self, durationVector, trajsubstr, -F,R,amin,amax,W,dW)

                #print durationVector,trajsubstr
                #print W,dW
                [semin,semax] = self.topp.getSpeedIntervalAtCriticalPoint(
                                Nc_in,
                                durationVector,
                                trajsubstr)
                print "TOPP velocity at CP",Nc_in,":",semin,semax
                ### AVP on the subtrajectory between 0 and Nc
                #if Nc > 0:
                #        [bspline, trajstr, durationVector] = self.computeTrajectorySubstringForTOPP(W, Nc)
                #        [semin, semax] = self.topp.getSpeedIntervalAtCriticalPoint(Nc, durationVector, trajstr)
                #else:
                #        semin = 0.0
                #        semax = 0.0
                #print "MAX SPEED AT CRITICAL PT",Nc,"/",Nwaypoints," is:",semin,semax
                return [semin,semax]

        def computeTrajectorySubstringForTOPP(self, Win, dWin, Nc):
                W = Win[:,0:Nc]
                dW = dWin[:,0:Nc]
                return self.computeTrajectoryStringForTOPP(W,dW)


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

        @classmethod
        def from_file(cls, filename):
                W = np.load(filename+'.npy')
                return cls(W)

        def get_number_waypoints(self):
                dt = self.DISCRETIZATION_TIME_STEP
                L = self.get_length()
                N = int(L/dt)
                return N

        def get_waypoints_second_order(self, N=None):
                ###############################################################
                ### obtain waypoints along the trajectory at constant spacing of
                ### DISCRETIZATION_TIME_STEP
                ###############################################################
                if N is None:
                        N = self.get_number_waypoints()
                ###############################################################
                K = self.get_dimension()
                pts = np.zeros((K,N))
                dpts = np.zeros((K,N))
                ddpts = np.zeros((K,N))
                ctr = 0
                for t in np.linspace(0.0, 1.0, num=N):
                        [f0,df0] = self.evaluate_at(t,der=1)
                        pts[:,ctr] = f0
                        dpts[:,ctr] = df0
                        #ddpts[:,ctr] = ddf0
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

        def execute(self, env, robot, tsleep=0.01, stepping=False):
                tstep = 0.01
                xt = self.topp.traj0

                with env.env:
                        robot.GetLinks()[0].SetStatic(True)
                        env.env.StopSimulation() 

                t = 0.0
                tstep = 0.01
                robot.SetDOFValues(xt.Eval(t))
                env.MakeRobotVisible()

                while t < xt.duration:
                        q = xt.Eval(t)
                        dq = xt.Evald(t)
                        ddq = xt.Evaldd(t)

                        qn = q + tstep*dq + 0.5*tstep*tstep*ddq
                        robot.SetDOFValues(qn)

                        env.env.StepSimulation(tstep)
                        time.sleep(tsleep)
                        t += tstep
                        if stepping:
                                raw_input('Press Key to Step. Time: '+str(t)+'/'+str(xt.duration))

        def SimpleInterpolate(self,q0,q1,qd0,qd1,T):
                a=((qd1-qd0)*T-2*(q1-q0-qd0*T))/T**3
                b=(3*(q1-q0-qd0*T)-(qd1-qd0)*T)/T**2
                c=qd0
                d=q0
                return [d,c,b,a]

        def computeSplineFromWaypoints(self,W):
                Ndim = W.shape[0]
                Nwaypoints = W.shape[1]
                bspline=[]
                for i in range(0,Ndim):

                        tvec = np.linspace(0,1,Nwaypoints)
                        WP = W[i,:]

                        #dws = (W[i,1]-W[i,0])
                        #tws = tvec[1]
                        #dwe = (W[i,-1]-W[i,-2])
                        #twe = 1.0-tvec[-2]

                        #M = 10
                        #for k in range(1,M):
                        #        tvec = np.hstack((-k*tws,tvec,k*twe+1.0))
                        #        WP = np.hstack((W[i,0]-k*dws,WP,W[i,-1]+k*dwe))

                        trajectory = splrep(tvec,WP,k=self.POLYNOMIAL_DEGREE,s=self.SMOOTH_CONSTANT)
                        bspline.append(trajectory)
                return bspline

        def GetIntervalIndex(self, W):
                Ndim = W.shape[0]
                Nwaypoints = W.shape[1]
                idx = []
                epsilon = 1e-1
                i = 0
                idx.append(0)
                while i < Nwaypoints-1:
                        k = 1
                        d = 0.0
                        while d < epsilon:
                                if i+k > Nwaypoints-1:
                                        k = Nwaypoints-i
                                        break
                                #print "i",i,"k",k,"d",d,Nwaypoints
                                q0 = W[:,i]
                                q1 = W[:,i+k]
                                d = np.linalg.norm(q0-q1)
                                k+=1
                        k=k-1
                        i = i + k
                        idx.append(i)
                return np.array(idx)

        def computeTrajectoryStringForTOPP(self, W, dW, DEBUG=0):
                #idx = self.GetIntervalIndex(W)
                #W = W[:,idx]
                #dW = dW[:,idx]

                Ndim = W.shape[0]
                Nwaypoints = W.shape[1]
                Kcoeff = 2
                Ninterval = Nwaypoints-1
                P = np.zeros((Ninterval, Ndim, Kcoeff))

                #poly= PPoly.from_spline(self.bspline[0])
                #Ninterval = poly.c.shape[1]
                ##print "Ninterval:",Ninterval,"Nwaypoints:",Nwaypoints
                durationVector = np.zeros((Ninterval))

                #for i in range(0,Ndim):
                #        T = 1.0/float(Ninterval)
                #        durationVector[i] = T
                #        poly= PPoly.from_spline(self.bspline[i])
                #        dpoly = poly.derivative(1)
                #        ddpoly = poly.derivative(2)
                #        coeff = poly.c
                #        c = coeff[0]
                #        b = coeff[1]
                #        a = coeff[2]

                #        P[:,i,0] = a
                #        P[:,i,1] = b
                #        P[:,i,2] = c

                #        qq = a[0] + T*b[0] + T*T*c[0]
                #        print a,b,c
                #        print qq

                #print coeff.shape


                qd0 = np.zeros((Ndim))
                qd0 += dW[:,0]
                #print W[0,:]
                #print dW[0,:]
                for j in range(0,Ninterval):
                        q0 = W[:,j]
                        q1 = W[:,j+1]
                        qd0 = dW[:,j]
                        #qd0 += dW[:,j]
                        #qd0 += dW[:,j]
                        qd1 = dW[:,j+1]
                        #print q0,q1,qd0,qd1
                        #durationVector[j] = np.linalg.norm(q0-q1)
                        #[a,b,c,d]=self.SimpleInterpolate(q0,q1,qd0,qd1,T)
                        a = q0

                        dq = np.linalg.norm(q1-q0)
                        durationVector[j] = dq
                        T = durationVector[j]
                        qnd = (q1-q0)/dq
                        b = dq/T*qnd
                        #c = (qd1-qd0)/(T)
                        #c = (q1-q0-qd0*T)/(T*T)

                        #print qd0
                        #qd0 = qd0 + (qd1-qd0) + 2*c*T

                        #qd0 = qd0 + 2*c*T
                        #c = (qd1-qd0)/(2*T)
                        #P[j,i,0]=0
                        P[j,:,0] = a
                        P[j,:,1] = b
                        #P[j,:,2] = c
                        #P[j,:,3] = d
                        #P[j,2,1]=0
                        #P[j,2,2]=0
                        #P[j,2,3]=0

                        #qspline0 = a
                        #qspline1 = a+T*b+T*T*c+T*T*T*d
                        #print a,b,c,d,T,qspline1,q1

                #self.CheckPolynomial(W,P,durationVector)
                for i in range(0,durationVector.shape[0]):
                        duration = durationVector[i]
                        if i==0:
                                trajectorystring = str(duration)
                        else:
                                trajectorystring += "\n" + str(duration)
                        trajectorystring += "\n" + str(Ndim)

                        for j in range(Ndim):
                                trajectorystring += "\n"
                                trajectorystring += string.join(map(str,P[i,j,:]))

                return [trajectorystring, durationVector]

        def EvalPoly(self, P,i,t):
                Ndim = P.shape[1]
                Kcoeff = P.shape[2]
                x = np.zeros((Ndim))
                dx = np.zeros((Ndim))
                for k in range(0,Kcoeff):
                        x += P[i,:,k]*(t**k)
                        dx += k*P[i,:,k]*(t**(max(k-1,0)))
                return [x,dx]
        def CheckPolynomial(self, W, P, D):
                [Ninterval, Ndim, Kcoeff] = P.shape

                X = []
                dX = []
                tvec = []
                M = D.shape[0]
                Wplt = []
                WTplt =[]
                #M = 10
                tall = 0.0
                Dall = 0.0
                for i in range(0,M):
                        t=0.0
                        tstep = D[i]/1000

                        WTplt.append(Dall)
                        Wplt.append(P[i,:,1])
                        Dall+=D[i]
                        while t <= D[i]:
                                tvec.append(tall)

                                [x,dx] = self.EvalPoly(P,i,t)
                                #x = np.zeros((Ndim))
                                #dx = np.zeros((Ndim))
                                #for k in range(0,Kcoeff):
                                #        x += P[i,:,k]*(t**k)
                                #        dx += k*P[i,:,k]*(t**(max(k-1,0)))

                                #x = P[i,:,0] + P[i,:,1]*t + P[i,:,2]*t*t + P[i,:,3]*t*t*t
                                #dx = P[i,:,1] + 2*P[i,:,2]*t + 3*P[i,:,3]*t*t
                                X.append(x)
                                dX.append(dx)
                                t+=tstep
                                tall += tstep

                        [q1spl,qd1spl] = self.EvalPoly(P,i,D[i])
                        q1 = W[:,i+1]

                        if np.linalg.norm(q1-q1spl)>1e-5:
                                print "qnext mismatch"
                                print "q1",q1
                                print "q1spl",q1spl
                                sys.exit(0)

                        q0 = W[:,i]
                        [q0spl,qd0spl] = self.EvalPoly(P,i,0)

                        if np.linalg.norm(q0-q0spl)>1e-5:
                                print "q0 mismatch"
                                print "q0",q0
                                print "q0spl",q0spl
                                sys.exit(0)

                tvec=np.array(tvec)
                X=np.array(X).T
                dX=np.array(dX).T
                WTplt=np.array(WTplt).T
                Wplt=np.array(Wplt).T
                print WTplt.shape
                print Wplt.shape
                #plt.plot(tvec,dX[0,:],'-k',linewidth=3)
                plt.plot(X[0,:],X[1,:],'-r',linewidth=3)
                #plt.plot(tvec,X[0,:],'-r',linewidth=3)

                #plt.plot(WTplt,Wplt[0,:],'-ok',linewidth=3)
                #plt.plot(tvec,X[1,:],'-k',linewidth=3)
                #plt.plot(tvec,X[2,:],'-b',linewidth=3)
                #plt.plot(tvec,X[3,:],'-g',linewidth=3)
                #plt.plot(tvec,dX[0,:],'--r',linewidth=3)
                #plt.plot(tvec,dX[1,:],'--k',linewidth=3)
                #plt.plot(tvec,dX[2,:],'--b',linewidth=3)
                #plt.plot(tvec,dX[3,:],'--g',linewidth=3)
                #plt.plot(W[0,0:M+1],W[1,0:M+1],'ok',markersize=6)
                #plt.plot(P[0:M,0,0],P[0:M,1,0],'ob',markersize=6)
                plt.show()

        def InsertMissingWaypoints(self,Wnext,max_dist):
                ictr=0
                for i in range(0,Wnext.shape[1]-1):
                        d = np.linalg.norm(Wnext[:,i]-Wnext[:,i+1])
                        if d > max_dist:
                                ## insert some points
                                dW = Wnext[:,i+1]-Wnext[:,i]

                                dstep = max_dist
                                dall = np.linalg.norm(dW)
                                Nnew = int(dall/dstep) - 1

                                dcur = 0.0
                                for k in range(1,Nnew):
                                        Wstep = (Wnext[:,i+k]-Wnext[:,i])
                                        Wstep /= np.linalg.norm(Wstep)
                                        Wnew = Wnext[:,i] + k*dstep*Wstep
                                        Wnext = np.insert(Wnext, i+k, Wnew, axis=1)
                                        ictr+=1
                print "added",ictr
                return Wnext

        def RepairTrajectory(self,W,max_dist):
                ictr=0
                i=0
                while i < W.shape[1]-1:

                        d = np.linalg.norm(W[:,i+1]-W[:,i])

                        istep = 1
                        if d > max_dist:

                                dstep = max_dist
                                Nnew = int(d/dstep)

                                dcur = 0.0
                                #print i,d,"/",max_dist,"d/dstep",d/dstep,Nnew
                                for k in range(1,Nnew):
                                        Wstep = (W[:,i+k]-W[:,i])
                                        Wstep /= np.linalg.norm(Wstep)
                                        Wnew = W[:,i] + k*dstep*Wstep
                                        W= np.insert(W, i+k, Wnew, axis=1)
                                        ictr+=1 
                                istep = Nnew

                        i = i+istep
                print "new points:",ictr
                return W

