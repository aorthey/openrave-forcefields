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

        DISCRETIZATION_TIME_STEP = 0.005

        rave_traj = []
        traj = []
        waypoints = []
        handle = []

        #env_ptr = []
        SMOOTH_CONSTANT=0
        POLYNOMIAL_DEGREE=3
        MIN_NUMBER_WAYPOINTS = 5

        ##drawing parameters
        ptsize = 0.03
        critical_pt_size = 0.07
        linsize = 1.5
        FONT_SIZE = 20
        dVECTOR_LENGTH = 0.5
        trajectory_tangent_color = np.array((0.9,0.2,0.9,0.9))
        trajectory_color_deformation = np.array((0.9,0.2,0.9,0.9))

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

        def prettify(self, W):
                if W.ndim <= 1:
                        print "cannot create trajectory with only one waypoint"
                        sys.exit(1)
                Nwaypoints = W.shape[1]
                if Nwaypoints <= 3:
                        W = self.addMinimalWaypoints(W)
                        Nwaypoints = W.shape[1]
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
                #amin = np.array((-5,-5,-5))
                #amax = np.array((5,5,5))
                AM = 5
                amin = np.array((-AM,-AM,-AM))
                amax = np.array((AM,AM,AM))
                #amin = np.array((-1,-1,-1))
                #amax = np.array((1,1,1))

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


        def computeTrajectoryStringForTOPP(self, W):
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

                        P[:,i,0] = coeff[3,:-1]
                        P[:,i,1] = coeff[2,:-1]
                        P[:,i,2] = coeff[1,:-1]
                        P[:,i,3] = coeff[0,:-1]

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

        def get_dimension(self):
                [F,dF] = self.evaluate_at(0)
                return F.shape[0]

        def waypoint_to_force(self, env, W):
                Ndims = W.shape[0]
                pt = np.array(((W[0],W[1],-0.1,0.001)))
                F = np.zeros((Ndims))
                F[0:3] = env.GetForceAtX(pt)
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

        def computeDelta(self, W, dW, Qoff, Rmin, Rmax):
                Ndim = Rmin.shape[0]
                eta = cvx.Variable(1)
                xvar = cvx.Variable(Ndim)
                yvar = cvx.Variable(Ndim)

                constraints = []
                constraints.append( Qoff + Rmin <= xvar )
                constraints.append( xvar <= Qoff + Rmax ) ## Q
                constraints.append( W + eta*dW == yvar ) ## L
                constraints.append( eta > 0)
                objfunc = cvx.norm(xvar-yvar)

                objective = cvx.Minimize(objfunc)
                prob = cvx.Problem(objective, constraints)
                result = prob.solve(solver=cvx.SCS, eps=1e-5)
                #result = prob.solve(solver=SCS, use_indirect=True, eps=1e-2, verbose=True)
                #prob.solve(verbose=True, abstol_inacc=1e-2,reltol_inacc=1e-2,max_iters= 300, reltol=1e-2)

                delta = np.linalg.norm(xvar.value-yvar.value)

                Qnearest = np.array(xvar.value).flatten()
                Lnearest = np.array(yvar.value).flatten()

                dL = dW/np.linalg.norm(dW)
                dQ = (Qnearest-W)/np.linalg.norm(Qnearest-W)
                angle = acos(np.dot(dL,dQ))

                return [delta, angle, Lnearest, Qnearest]

        def drawReachableSet(self, env, dt, speed, W, dW, amin, amax, Lnearest, Qnearest, R, F, angle, delta):
                dt2 = dt*dt*0.5
                from scipy.spatial import ConvexHull
                import itertools
                Ndim = W.shape[0]

                witv = dt
                plot([W[0],W[0]+witv*dW[0]],[W[1],W[1]+witv*dW[1]],'-k',linewidth=3)
                plot([W[0]+witv*dW[0],W[0]+(witv+witv/3.0)*dW[0]],[W[1]+witv*dW[1],W[1]+(witv+witv/3.0)*dW[1]],'--k',linewidth=3)
                plot(W[0],W[1],'or',markersize=10)
                plot([Lnearest[0],Qnearest[0]],[Lnearest[1],Qnearest[1]],'--og',markersize=10,linewidth=4)
                plot([Lnearest[0]],[Lnearest[1]],'ok',markersize=20)
                Mlq = Lnearest+0.5*(Qnearest-Lnearest)
                strdelta = ' $\\delta$='+str(np.around(delta,decimals=2))
                plt.text(Mlq[0], Mlq[1], strdelta, fontsize=22)
                dQ = Qnearest-W
                plot([W[0],W[0]+dQ[0]],[W[1],W[1]+dQ[1]],'--g',linewidth=3)
                Mlq = W+0.2*(Qnearest-W)
                strangle = ' $\\theta$='+str(np.around(angle,decimals=2))
                plt.text(Mlq[0], Mlq[1], strangle, fontsize=22)


                ctr = 0
                Adim = amin.shape[0]
                Mf = 2**Adim
                Qt = np.zeros((Mf,Ndim))
                for val in itertools.product([0,1], repeat=Adim):
                        ## if 0, then choose amin, if 1, then choose amax
                        a = val*amax + abs(1-np.array(val))*amin
                        Qt[ctr,:] = W + speed*dt*dW + dt2*F.flatten() + dt2*np.dot(R,a) 

                        ctr=ctr+1

                #plt.fill_betweenx(Qt[:,0], Qt[:,1])
                plot(Qt[:,0],Qt[:,1],'ob',markersize=5)
                hull = ConvexHull(Qt[:,0:2])
                plt.fill(Qt[hull.vertices,0], Qt[hull.vertices,1], facecolor='blue', alpha=0.2, lw=2)
                plt.show()

        def visualizeReachableSetsAtT(self, env, t0, dt = 0.1):
                speed = 0.0
                dt2 = dt*dt*0.5
                [W,dW,ddW] = self.evaluate_at(t0,der=2)
                print "REACHABLE SET AT POSITION:",W
                [Ndim, Nwaypoints] = self.getWaypointDim(W)
                F = self.get_forces_at_waypoints(W, env)
                #F = np.zeros((F.shape))
                [R,amin,amax] = self.getControlMatrix(W)
                Rmin = np.zeros((Ndim))
                Rmax = np.zeros((Ndim))
                Rmin = dt2*np.minimum(np.dot(R[:,:,0],amax),np.dot(R[:,:,0],amin))
                Rmax = dt2*np.maximum(np.dot(R[:,:,0],amax),np.dot(R[:,:,0],amin))

                Qoff = W + speed*dt*dW + dt2*F.flatten()
                [delta0, angle0, Lnearest, Qnearest] = self.computeDelta(W, dW, Qoff, Rmin, Rmax)
                self.drawReachableSet(env, dt, speed, W, dW, amin, amax, Lnearest, Qnearest, R[:,:,0], F, angle0, delta0)
                
                M = 50
                lambda1 = np.linspace(-0.02,0.02,M)
                lambda2 = np.linspace(0.0,1,M)

                delta = np.zeros((M,M))
                angle = np.zeros((M,M))
                ctri = 0
                for l1 in lambda1:
                        ctrj = 0
                        for l2 in lambda2:
                                speed = sqrt(l2*amax[0]*2)
                                trajNormal = np.array((-W[1],W[0],0.0,0.0))
                                trajNormal = trajNormal/np.linalg.norm(trajNormal)

                                Fproj = np.dot(F.flatten().T, trajNormal)
                                Qoff = W + speed*dt*dW + dt2*F.flatten() - l1*Fproj*trajNormal
                                [delta[ctri,ctrj], angle[ctri,ctrj], Lnearest, Qnearest] = self.computeDelta(W, dW, Qoff, Rmin, Rmax)
                                ctrj = ctrj+1
                        ctri = ctri+1

                from mpl_toolkits.mplot3d import Axes3D
                from matplotlib import cm
                from matplotlib.ticker import LinearLocator, FormatStrFormatter
                import matplotlib.pyplot as plt

                fig = plt.figure()
                ax = fig.gca(projection='3d')
                X, Y = np.meshgrid(lambda2, lambda1)
                surf = ax.plot_surface(X, Y, angle, rstride=1, cstride=1, cmap=cm.coolwarm,
                                       linewidth=0, antialiased=False)

                ax.set_xlabel('$\\lambda_2$',fontsize=22)
                ax.set_ylabel('$\\lambda_1$',fontsize=22)
                ax.set_zlabel('$\\theta$',fontsize=22)

                dmax = np.amax(delta)
                ax.plot([0],[0],[angle0],'ob',markersize=10)
                #ax.set_zlim(-1.01, 1.01)

                #fig.colorbar(surf, shrink=0.5, aspect=5)

                plt.show()
                #print delta

        def computeReachableSets(self, dt, env):
                [W,dW,ddW] = self.get_waypoints_second_order()

                [Ndim, Nwaypoints] = self.getWaypointDim(W)

                F = self.get_forces_at_waypoints(W, env)
                [R,amin,amax] = self.getControlMatrix(W)

                Rmin = np.zeros((Nwaypoints, Ndim))
                Rmax = np.zeros((Nwaypoints, Ndim))
                for i in range(0,Nwaypoints):
                        Rmin[i,:] = np.minimum(np.dot(R[:,:,i],amax),np.dot(R[:,:,i],amin))
                        Rmax[i,:] = np.maximum(np.dot(R[:,:,i],amax),np.dot(R[:,:,i],amin))

                dt2 = dt*dt*0.5

                delta = np.zeros((Nwaypoints))
                eta = cvx.Variable(Nwaypoints)
                xvar = cvx.Variable(Ndim,Nwaypoints)
                yvar = cvx.Variable(Ndim,Nwaypoints)
                constraints = []
                objfunc = 0.0
                for i in range(0,Nwaypoints):
                        Qoff = W[:,i] + dt*dW[:,i] + dt2*F[:,i]
                        constraints.append( Qoff + Rmin[i,:] <= xvar[:,i] )
                        constraints.append( xvar[:,i] <= Qoff + Rmax[i,:] ) ## Q
                        constraints.append( W[:,i] + eta[i]*dW[:,i] == yvar[:,i] ) ## L
                        constraints.append( eta[i] >= 0)
                        objfunc += cvx.norm(xvar[:,i]-yvar[:,i])

                objective = cvx.Minimize(objfunc)
                prob = cvx.Problem(objective, constraints)

                print "solve minimal speed"
                result = prob.solve(solver=cvx.SCS)
                for i in range(0,Nwaypoints):
                        delta[i] = np.linalg.norm(xvar[:,i].value-yvar[:,i].value)

                plt.clf()
                T = np.linspace(0,1.0,Nwaypoints)
                subplot(2,1,1)
                title('Delta and Force Profiles', fontsize=self.FONT_SIZE)
                plot(T,delta,'-r',linewidth=4)
                ylabel('Delta $\delta$', fontsize=self.FONT_SIZE)
                subplot(2,1,2)
                plot(T,F[0,:],'-b',linewidth=4,label="$F_x$")
                plot(T,F[1,:],'-g',linewidth=4,label="$F_y$")
                xlabel('Path parameter $s$',fontsize=self.FONT_SIZE)
                ylabel('Force $F$', fontsize=self.FONT_SIZE)
                legend = plt.legend(loc='upper center', shadow=True, fontsize=self.FONT_SIZE)
                plt.show()

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

        def getCriticalPointFromWaypoints(self, env, W, dW, ddW, oldNc = 0):

                [Ndim, Nwaypoints] = self.getWaypointDim(W)
                F = self.get_forces_at_waypoints(W, env)
                [R,amin,amax] = self.getControlMatrix(W)

                self.topp = TOPPInterface(self, self.durationVector, self.trajectorystring, -F,R,amin,amax,W,dW)
                Nc = self.topp.getCriticalPoint()
                if Nc < 0:
                        print "return oldNc=",oldNc
                        return oldNc
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
                rave_traj = ravetraj
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
                with env.env:
                        for i in range(0,Nwaypoints):
                                pt = np.array(((W[0,i],W[1,i],W[2,i])))
                                tmp_handle.append(env.env.plot3(points=pt,
                                                   pointsize=self.ptsize,
                                                   colors=self.trajectory_color,
                                                   drawstyle=1))
                                dpt = np.array(((dW[0,i],dW[1,i],dW[2,i])))
                                ddpt = np.array(((ddW[0,i],ddW[1,i],ddW[2,i])))
                                dpt = self.dVECTOR_LENGTH*dpt/np.linalg.norm(dpt)
                                #ddpt = self.dVECTOR_LENGTH*ddpt/np.linalg.norm(ddpt)

                                P = np.array(((pt[0],pt[1],pt[2]),
                                        (pt[0]+dpt[0],pt[1]+dpt[1],pt[2]+dpt[2])))
                                h=env.env.drawlinestrip(points=P,linewidth=self.linsize,colors=self.trajectory_tangent_color)
                                tmp_handle.append(h)
                                #P = np.array(((pt[0],pt[1],pt[2]),
                                #        (pt[0]+ddpt[0],pt[1]+ddpt[1],pt[2]+ddpt[2])))
                                #h=env.env.drawlinestrip(points=P,linewidth=self.linsize,colors=np.array(((0.9,0.9,0.9,0.9))))
                                #tmp_handle.append(h)

                return tmp_handle

        def draw(self, env, keep_handle=True):
                [W,dW,ddW] = self.get_waypoints_second_order()
                N = self.getCriticalPointFromWaypoints(env, W, dW, ddW)
                tmp_handle = []

                self.trajectory_color = self.trajectory_color_feasible
                tmp_handle.append(self.get_handle_draw_waypoints(env, W[:,0:N], dW[:,0:N], ddW[:,0:N]))
                self.trajectory_color = self.trajectory_color_infeasible
                tmp_handle.append(self.get_handle_draw_waypoints(env, W[:,N:], dW[:,N:], ddW[:,N:]))

                if keep_handle:
                        self.handle = tmp_handle
                else:
                        return tmp_handle

        def draw_delete(self):
                self.handle = []
