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
import copy

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

        show_tangent_vector = False
        show_orientation_vector = False

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
        Ndim = 0

        ### waypoints: real matrix of dimensionality Ndim x Nwaypoints
        def __init__(self, waypoints_in):
                self.new_from_waypoints(waypoints_in)

        def new_from_waypoints(self, waypoints_in):
                self.waypoints = self.prettify(waypoints_in)
                self.Ndim = self.waypoints.shape[0]
                self.Nwaypoints = self.waypoints.shape[1]
                self.bspline = self.computeSplineFromWaypoints(self.waypoints)
                [self.waypoints,dW,ddW] = self.get_waypoints_second_order()

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
                ACCURACY_DUPLICATE = self.DISCRETIZATION_TIME_STEP/1e10
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

        def info(self):
                print "#### TRAJECTORY CLASS ######"
                print "LENGTH: ", self.get_length()
                print "DISCRETIZATION:",self.DISCRETIZATION_TIME_STEP
                print "WAYPOINTS: ", self.waypoints.shape[1]

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

                        ##sailboat
                        force = np.array((2.5,-3.5,0,5.0))
                        dp = np.array((1,0.0,0,0.2))
                        ##car

                        ds= 0.01

                        ### extra large dp
                        p =np.array( [-2.608897298112661, -0.2047366548132623, 0.1000000008651387, -2.987014267911084] )
                        dp =np.array( [-16.487757154509474, -4.837209926510647, 0.0, 4.899656054155726] )
                        force =np.array( [0.0, 3.0, 0.0, -1.4821147977289577] )
                        speed= 3.86342505425
                        ### no force, no dp
                        p =np.array( [-2.644855145592143, 1.1657167012260294, 0.10000000000235312, -3.7242126663188597] )
                        dp =np.array( [-1.1866896124736256e-20, 1.910882585951823e-18, 1.4346537045461183e-27, -9.702129945467387e-19] )
                        force =np.array( [0.0, 0.0, 0.0, 0.0] )
                        speed= 0.341259539422
                        ### no force, no dp
                        p =np.array( [-2.644855145592143, 1.1657167012260294, 0.10000000000235312, -3.7242126663188597] )
                        dp =np.array( [0.0, 0.0, 0.0, 0.0] )
                        force =np.array( [0.0, 0.0, 0.0, 0.0] )
                        speed= 0.0

                        ### close to large dp
                        p =np.array( [-4.113833032000303, -0.0026206376782298216, 0.0999999999999477, -3.139440349128482] )
                        dp =np.array( [-982166.3702744454, -1372.8293479195447, -8.057068259496254e-08, 1124.7453741523088] )
                        force =np.array( [0.0, 0.0, 0.0, 0.0] )
                        speed= 1.07398453113

                        p = np.array( [ -5.05e+00, -4.08e-03, 1.00e-01, -3.14e+00] )
                        dp = np.array([ -2.97e+08, -4.32e+05, 5.76e-06, 3.48e+05] )
                        F = np.array( [ 0.       ,      0.  , 0.      , 0.])
                        speed = 1.07398453113 
                        #ds 1.19162456455 
                        #dt 2.19793777714e-11
                        p =np.array( [-2.5974191119628514, -0.0004997370983040311, 0.0999999999999989, -3.141179458965877] )
                        dp =np.array( [-1.2469468647720259, -0.002999767139199421, 0.0, 0.0017880891982438346] )
                        force =np.array( [0.0, 0.0, 0.0, 0.0] )
                        speed= 1.51000061035
                        #ds: 0.00305489294111

                        ### nice vis
                        p = np.array((0,1e-4,0,0))
                        force = np.array((0.5,-3.5,0,1.0))
                        dp = np.array((1,0.1,0,0.2))
                        speed = 0.2

                        ##### error
                        ###p =np.array( [-6.012902723369949, -0.005123021667037449, 0.1, -3.1373188699666956] )
                        ###force =np.array( [0.0, 0.0, 0.0, 0.0] )
                        ###dp =np.array( [-344.35813853264824, -0.46681954515138263, 0.0, 0.389387904541571] )
                        ###speed= 1.19162642333

                        #from reachable_set import ReachableSet
                        #self.reach = ReachableSet( p, s, dp, force, R[:,:,0], amin, amax)
                        #self.reach.Plot()
                        [R,amin,amax] = self.getControlMatrix(p)

                        from reachable_set3d import ReachableSet3D
                        self.reach = ReachableSet3D(env)
                        dnp = dp/np.linalg.norm(dp)
                        dpnormal = -0.001*(force - np.dot(force,dnp)*dnp)

                        self.reach.PlotSingleSet(self.DISCRETIZATION_TIME_STEP, 
                                        p, dp, speed, force, 
                                        R[:,:,0], amin, amax)
                        #self.reach.PlotTotalSet(self.DISCRETIZATION_TIME_STEP, 
                        #                p, dp, speed, force, 
                        #                R[:,:,0], amin, amax)
                        #self.reach.PlotMoveAgainstWrenchField(self.DISCRETIZATION_TIME_STEP, 
                        #                p, dpnormal, dp, speed, force, 
                        #                R[:,:,0], amin, amax)
                        #self.reach.PlotStretch(self.DISCRETIZATION_TIME_STEP, 
                        #                p, dp, speed, force, 
                        #                R[:,:,0], amin, amax)

        def get_dimension(self):
                [F,dF] = self.evaluate_at(0)
                return F.shape[0]

        def waypoint_to_force(self, env, W):
                Ndims = W.shape[0]
                pt = np.array(((W[0],W[1],-0.1,W[3])))
                F = np.zeros((Ndims))
                F[0:3] = env.GetForceAtX(pt)
                r = 0.3

                theta = W[3]
                rcom = r * np.dot(Rz(theta),ex)
                torque = np.cross(rcom,F[0:2])

                ### TORQUE application
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

        def get_tangent_forces_at_waypoints(self, W, dW, env):
                F = self.get_forces_at_waypoints(W, env)
                if F.ndim>1:
                        Nwaypoints = F.shape[1]
                else:
                        Nwaypoints = 1

                FT = copy.copy(F)
                for i in range(0,Nwaypoints):
                        T = dW[:,i]/np.linalg.norm(dW[:,i])
                        FT[:,i] = np.dot(F[:,i],T)*T

                return FT

        def get_normal_forces_at_waypoints(self, W, dW, env):
                F = self.get_forces_at_waypoints(W, env)
                if F.ndim>1:
                        Nwaypoints = F.shape[1]
                else:
                        Nwaypoints = 1

                FN = copy.copy(F)
                for i in range(0,Nwaypoints):
                        T = dW[:,i]/np.linalg.norm(dW[:,i])
                        FN[:,i] = F[:,i] - np.dot(F[:,i],T)*T

                return FN

        def get_forces_at_waypoints(self, W, env):
                Ndim = W.shape[0]
                if W.ndim>1:
                        Nwaypoints = W.shape[1]
                else:
                        F = self.waypoint_to_force(env,W)
                        return F

                F = np.zeros((Ndim,Nwaypoints))
                for i in range(0,Nwaypoints):
                        F[:,i] = self.waypoint_to_force(env, W[:,i])
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
                self.topp = TOPPInterface(W, env)
                if self.topp.ReparameterizeTrajectory():
                        self.topp.PlotTrajectory(env)
                        print self.topp
                else:
                        print "Trajectory has no ReParameterization"


        def getCriticalPointFromWaypoints(self, env, W, oldNc = 0):
                self.topp = TOPPInterface(W, env)
                Nc = self.topp.getCriticalPoint()-1
                if Nc < 0:
                        print "return oldNc=",oldNc
                        Nc = oldNc
                return Nc

        def getCriticalPoint(self, env):
                #[W,dW,ddW] = self.get_waypoints_second_order()
                N = self.getCriticalPointFromWaypoints(env, self.waypoints)
                if N is None:
                        print "No critical point found"
                        sys.exit(1)
                return N

        def getVelocityIntervalWithoutForceField(self, env, W):

                self.topp = TOPPInterface(W, env, zeroForce=True)

                if self.topp.ReparameterizeTrajectory():
                        return self.topp.traj0
                else:
                        #### without a force field any path should be executable
                        #### at near zero speed => catch it if not
                        print "WARNING1: without force field, TOPP couldn't find a valid \
                        velocity profile. Path not continuous or system not STLC"
                        print "discrtimestep=",self.topp.discrtimestep

                        self.topp.SaveToFile('clc2')

                        CP = self.topp.getCriticalPoint()
                        traj = self.topp.traj0
                        q = np.array([traj.Eval(t) for t in np.linspace(0,traj.duration,1e5)]).T
                        plot(q[0,:],q[1,:],'-r',linewidth=2)
                        plot(W[0,:],W[1,:],'ok',markersize=2)
                        plot(W[0,CP],W[1,CP],'og',markersize=10)
                        plt.show()
                        sys.exit(1)

        def getTOPPTrajectoryWithoutForceField(self, env, W):
                self.topp = TOPPInterface(W, env, zeroForce=True)
                if self.topp.ReparameterizeTrajectory():
                        return self.topp
                else:
                        self.info()
                        print "WARNING2: without force field, TOPP couldn't find a valid \
                        velocity profile. Path not continuous or system not STLC"
                        print W.shape
                        print W
                        sys.exit(1)

        #def getTOPPTrajectory(self, env, W, dW, F):

        #        self.topp = TOPPInterface(W, env)
        #        if self.topp.ReparameterizeTrajectory():
        #                return self.topp
        #        else:
        #                print W.shape
        #                print W
        #                self.info()
        #                sys.exit(1)
        #                return None


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
                        #[Ndim,Nwaypoints] = self.getWaypointDim(self.waypoints)
                        #N = Nwaypoints
                ###############################################################
                K = self.get_dimension()
                pts = np.zeros((K,N))
                dpts = np.zeros((K,N))
                ddpts = np.zeros((K,N))
                ctr = 0
                for t in np.linspace(0.0, self.get_length(), num=N):
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
                return self.length
                #dd = 0.0
                #T = np.linspace(0.0, 1.0, num=self.waypoints.shape[1])
                #for i in range(0,len(T)-1):
                #        t = T[i]
                #        tdt = T[i+1]
                #        [ft,df0] = self.evaluate_at(t)
                #        [ftdt,df0] = self.evaluate_at(tdt)
                #        dd = dd + np.linalg.norm(ft-ftdt)
                #return dd

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

        def draw(self, env, keep_handle=True, critical_pt = None, DrawRobot=False):

                [W,dW,ddW] = self.get_waypoints_second_order()
                t1 = time.time()
                if critical_pt == None:
                        N = self.getCriticalPointFromWaypoints(env, W)
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

        def draw_robot_along_path(self, env, robot, N=10):
                xt = self.topp.traj0
                with env.env:
                        robot.GetLinks()[0].SetStatic(True)
                        env.env.StopSimulation() 

                env.MakeRobotVisible()

                ictr=0
                print "before:",env.env.GetRobots()
                R = env.env.GetRobot("clone"+str(ictr))
                while R is not None:
                        env.env.Remove(R)
                        ictr+=1
                        R = env.env.GetRobot("clone"+str(ictr))
                print "after deleting:",env.env.GetRobots()

                with env.env:
                        ictr=0
                        for t in np.linspace(0,xt.duration,N):
                                print "t",t,"/",xt.duration
                                robot_t = RaveCreateRobot(env.env,robot.GetXMLId())
                                robot_t.Clone(robot,0)
                                robot_t.SetName("clone"+str(ictr))
                                env.env.AddRobot(robot_t,True)
                                #env.ChangeTransparencyRobot(robot_t, 1.0)
                                robot_t.SetDOFValues(xt.Eval(t))
                                ictr+=1
                        #robot_t.SetDOFValues(xt.Eval(t))
                print "after adding:",env.env.GetRobots()

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

                q = xt.Eval(0)
                dq = xt.Evald(0)
                while t < xt.duration:
                        #q = xt.Eval(t)
                        #dq = xt.Evald(t)
                        ddq = xt.Evaldd(t)

                        #qn = q + tstep*dq + 0.5*tstep*tstep*ddq
                        #robot.SetDOFValues(qn)
                        q = q + tstep*dq + 0.5*tstep*tstep*ddq
                        dq = dq + tstep*ddq
                        robot.SetDOFValues(q)

                        env.env.StepSimulation(tstep)
                        time.sleep(tsleep)
                        t += tstep
                        if stepping:
                                raw_input('Press Key to Step. Time: '+str(t)+'/'+str(xt.duration))

        def computeSplineFromWaypoints(self,W):
                Ndim = W.shape[0]
                Nwaypoints = W.shape[1]
                bspline=[]

                tvec = np.linspace(0,1,Nwaypoints)
                for i in range(Nwaypoints-1):
                            tvec[i+1] = tvec[i]+linalg.norm(W[:,i]-W[:,i+1])

                self.length = tvec[i+1]

                for i in range(0,Ndim):
                        WP = W[i,:]
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

