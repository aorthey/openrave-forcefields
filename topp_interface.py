import string,time
from pylab import *
import numpy as np
from openravepy import *
import TOPP
from TOPP import Utilities
from TOPP import TOPPbindings
from TOPP import TOPPpy
from TOPP import Trajectory
from TOPP import TOPPopenravepy

class TOPPInterface():
        #DURATION_DISCRETIZATION = 0.0001
        #DURATION_DISCRETIZATION = 1
        DURATION_DISCRETIZATION = 0.1

        TRAJECTORY_ACCURACY_REQUIRED = 1e-1
        traj0 = []
        trajstr = []
        durationVector = []
        length = -1
        Ndim = -1
        Nwaypoints = -1
        critical_point = -1
        critical_point_value = -1.0
        topp_inst = []
        semin = -1
        semax = -1
        F_ = []
        R_ = []
        amin_ = []
        amax_ = []
        W_ = []
        dW_ = []
        ddW_ = []

        def initializeFromSpecifications(self, F, R, amin, amax, W, dW, ddW):
                self.Ndim = W.shape[0]
                self.Nwaypoints = W.shape[1]

                [self.trajstr,self.durationVector] = self.GetTrajectoryString(W, dW, ddW)
                self.traj0 = Trajectory.PiecewisePolynomialTrajectory.FromString(self.trajstr)
                self.length = np.sum(self.durationVector)

                dendpoint = np.linalg.norm(self.traj0.Eval(self.length)-W[:,-1])

                if dendpoint > self.TRAJECTORY_ACCURACY_REQUIRED:
                        print self.length
                        print "###############"
                        print "FINAL POINT on piecewise C^2 traj:",self.traj0.Eval(self.length)
                        print "FINAL WAYPOINT                   :",W[:,-1]
                        print "###############"
                        sys.exit(1)

                [a,b,c] = self.GetABC(F, R, amin, amax)

                vmax = 1000*np.ones(self.Ndim)
                durationQ = self.traj0.duration/(self.Nwaypoints-1)
                self.topp_inst = TOPP.QuadraticConstraints(self.traj0, durationQ, vmax, list(a), list(b), list(c))

        def __init__(self, F, R, amin, amax, W, dW, ddW):
                self.F_ = F
                self.R_ = R
                self.amin_ = amin
                self.amax_ = amax
                self.W_ = W
                self.dW_ = dW
                self.ddW_ = ddW
                self.initializeFromSpecifications(F, R, amin, amax, W, dW, ddW)

        def getSpeedIntervalAtCriticalPoint(self, N):
                assert(N > -1)
                Fc = self.F_[:,0:N]
                Rc = self.R_[:,:,0:N]
                Wc = self.W_[:,0:N]
                dWc = self.dW_[:,0:N]
                ddWc = self.ddW_[:,0:N]
                self.initializeFromSpecifications(Fc, Rc, self.amin_, self.amax_, Wc, dWc, ddWc)
                x = self.topp_inst.solver
                [self.semin,self.semax] = self.topp_inst.AVP(0, 0)
                return [self.semin,self.semax]

        def PlotTrajectory(self):

                #################################
                #Npts = 600
                #tvect = np.linspace(0,self.traj0.duration, Npts)
                #qdvect = np.array([self.traj0.Evald(t) for t in tvect])
                #np.set_printoptions(suppress=True)
                #np.set_printoptions(precision=2)
                #print self.traj0.Eval(0.1900)
                #print self.traj0.Eval(0.1900)
                #print self.traj0.Eval(0.1900)
                #print self.traj0.Eval(0.1900)
                #print tvect.shape
                #print qdvect.shape
                #print np.vstack((qdvect[0:10,0:2].T,tvect[0:10])).T
                #Npts = 66
                #tvect = np.linspace(0,self.traj0.duration, Npts)
                #qdvect = np.array([self.traj0.Evald(t) for t in tvect])
                #print self.traj0.Eval(0.19)
                #print self.traj0.Eval(0.19)
                #print self.traj0.Eval(0.19)
                #print self.traj0.Eval(0.19)
                #
                dt = float(self.DURATION_DISCRETIZATION)/2.0
                Npts = int(self.traj0.duration/dt)
                print Npts
                tvect = np.linspace(0,self.traj0.duration, Npts)
                qvect = np.array([self.traj0.Eval(t) for t in tvect])
                qdvect = np.array([self.traj0.Evald(t) for t in tvect])
                qddvect = np.array([self.traj0.Evaldd(t) for t in tvect])

                #################################
                twvect = np.linspace(0,self.traj0.duration, self.Nwaypoints)

                lw = 3
                fs = 22
                lloc = 'lower right'
                ANCHOR = (0.5,1)
                
                subplot(3,1,1)

                plot(twvect, self.W_[0,:], '-k', linewidth = lw)
                plot(twvect, self.W_[1,:], '-k', linewidth = lw)
                plot(twvect, self.W_[2,:], '-k', linewidth = lw)
                plot(twvect, self.W_[3,:], '-k', linewidth = lw)
                plot(tvect, qvect[:,0], linewidth = lw, label = "$x$")
                plot(tvect, qvect[:,1], linewidth = lw, label = "$y$")
                plot(tvect, qvect[:,2], linewidth = lw, label = "$z$")
                plot(tvect, qvect[:,3], linewidth = lw, label = "$\\theta$")
                #plot(tvect, qdvect, f, linewidth=2)
                #plot(tvect, qddvect, f, linewidth=2)
                title('Position/Velocity/Acceleration Path', fontsize=fs)
                ylabel('Position', fontsize=fs)
                legend = plt.legend(loc=lloc, shadow=True, fontsize=fs)

                subplot(3,1,2)
                ylabel('Velocity', fontsize=fs)
                plot(tvect, qdvect[:,0], linewidth = lw, label = "$\dot x$")
                plot(tvect, qdvect[:,1], linewidth = lw, label = "$\dot y$")
                plot(tvect, qdvect[:,2], linewidth = lw, label = "$\dot z$")
                plot(tvect, qdvect[:,3], linewidth = lw, label = "$\dot \\theta$")
                legend = plt.legend(loc=lloc, shadow=True, fontsize=fs)

                subplot(3,1,3)
                ylabel('Acceleration', fontsize=fs)
                plot(tvect, qddvect[:,0], linewidth = lw, label = "$\ddot{x}$")
                plot(tvect, qddvect[:,1], linewidth = lw, label = "$\ddot{y}$")
                plot(tvect, qddvect[:,2], linewidth = lw, label = "$\ddot{z}$")
                plot(tvect, qddvect[:,3], linewidth = lw, label = "$\ddot{\\theta}$")
                xlabel('Time $t$',fontsize=fs)
                #box = ax.get_position()
                #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
                #legend = plt.legend(loc=lloc, shadow=True, fontsize=fs, bbox_to_anchor=ANCHOR)

                plt.show()

        def ReparameterizeTrajectory(self):
                x = self.topp_inst.solver
                ret = x.RunComputeProfiles(0.0,0.0)
                if ret == 1:
                        x.ReparameterizeTrajectory()
                        x.WriteResultTrajectory()
                        self.traj0 = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
                        return True
                else:
                        return False

        def getCriticalPoint(self):
                x = self.topp_inst.solver
                #x.integrationtimestep = 0.001
                #x.reparamtimestep = 0.001
                #x.extrareps = 10

                self.critical_point = self.Nwaypoints
                try:
                        ret = x.RunComputeProfiles(0.0,0.0)
                        if ret == 4:
                                #[semin,semax] = topp_inst.AVP(0, 0)
                                #print "TOPP critical pt:",semin,semax
                                self.critical_point = x.GetCriticalPoint()
                                self.critical_point_value = x.GetCriticalPointValue()
                                #print "TOPP critical pt:",self.critical_point,self.critical_point_value
                                return self.critical_point
                        if ret == 1:
                                print "TOPP: success"
                                #x.ReparameterizeTrajectory()
                        
                                #ret_param = x.ReparameterizeTrajectory()
                                #print "reparametrized"
                                #x.WriteProfilesList()
                                #x.WriteSwitchPointsList()
                                #profileslist = TOPPpy.ProfilesFromString(x.resprofilesliststring)
                                #switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
                                #TOPPpy.PlotProfiles(profileslist,switchpointslist,4)
                                #x.WriteResultTrajectory()
                                #self.traj0 = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
                                #print "reparametrized"
                                #TOPPpy.PlotKinematics(self.traj0,traj1)

                                return self.critical_point
                        if ret== 0:
                                M = 10
                                np.set_printoptions(precision=3)
                                dt = self.DURATION_DISCRETIZATION
                                tvect = np.linspace(0,self.traj0.duration, self.traj0.duration/dt)
                                qvect = np.array([self.traj0.Eval(t) for t in tvect])
                                qdvect = np.array([self.traj0.Evald(t) for t in tvect])
                                qddvect = np.array([self.traj0.Evaldd(t) for t in tvect])
                                print "TOPP: unspecified error"
                                print "W=",repr(self.W_[0:2,0:M])
                                print "q=",repr(qvect[0:M,0:2].T)
                                print "dW=",repr(self.dW_[0:2,0:M])
                                print "dq=",repr(qdvect[0:M,0:2].T)
                                print "ddW=",repr(self.ddW_[0:2,0:M])
                                print "ddq=",repr(qddvect[0:M,0:2].T)
                                self.PlotTrajectory()
                                sys.exit(0)
                        else:
                                print "TOPP: ",ret
                                sys.exit(0)

                except Exception as e:
                        print "TOPP EXCEPTION: ",e
                        #print self.durationVector
                        print self.traj0.duration
                        sys.exit(0)
                        #return -1

        def GetABC(self, F, R, amin, amax):
                ### compute a,b,c
                q = np.zeros((self.Ndim,self.Nwaypoints))
                qs = np.zeros((self.Ndim,self.Nwaypoints))
                qss = np.zeros((self.Ndim,self.Nwaypoints))
                for i in range(0,self.Nwaypoints):
                        duration = np.sum(self.durationVector[0:i])
                        q[:,i] = self.traj0.Eval(duration)
                        qs[:,i] = self.traj0.Evald(duration)
                        qss[:,i] = self.traj0.Evaldd(duration)

                I = np.identity(self.Ndim)
                G = np.vstack((I,-I))

                a = np.zeros((self.Nwaypoints, 2*self.Ndim))
                b = np.zeros((self.Nwaypoints, 2*self.Ndim))
                c = np.zeros((self.Nwaypoints, 2*self.Ndim))

                for i in range(0,self.Nwaypoints):
                        Rmax = np.maximum(np.dot(R[:,:,i],amin),np.dot(R[:,:,i],amax))
                        Rmin = np.minimum(np.dot(R[:,:,i],amin),np.dot(R[:,:,i],amax))
                        H1 = F[:,i] - Rmax
                        H2 = -F[:,i] + Rmin
                        for j in range(self.Ndim):
                                if H2[j] > -H1[j]:
                                        print H2[j],"<= q[",j,"]<=",-H1[j]
                                        sys.exit(1)

                        c[i,:] = np.hstack((H1,H2)).flatten()

                for i in range(0,self.Nwaypoints):
                        a[i,:] = np.dot(G,qs[:,i]).flatten()
                        b[i,:] = np.dot(G,qss[:,i]).flatten()
                return [a,b,c]

        def GetTrajectoryString(self, W, dW, ddW):
                ### get trajectory string for any R subspaces
                Ndim = W.shape[0]
                Nwaypoints = W.shape[1]
                durationVector = np.zeros((Nwaypoints-1))

                DEBUG = 0
                for i in range(0,Nwaypoints-1):
                        #ds = np.linalg.norm(W[:,i+1]-W[:,i])
                        #dv = np.linalg.norm(dW[:,i])
                        #duration = ds/dv
                        duration = self.DURATION_DISCRETIZATION
                        durationVector[i] = duration
                        if i==0:
                                trajectorystring = str(duration)
                        else:
                                trajectorystring += "\n" + str(duration)

                        trajectorystring += "\n" + str(Ndim)

                        A = np.zeros(Ndim)
                        B = np.zeros(Ndim)
                        C = np.zeros(Ndim)
                        D = np.zeros(Ndim)
                        E = np.zeros(Ndim)
                        F = np.zeros(Ndim)

                        T = duration
                        TT = duration*duration
                        TTT = duration*duration*duration

                        p0 = W[:,i]
                        dp0 = dW[:,i]
                        ddp0 = ddW[:,i]
                        p1 = W[:,i+1]
                        dp1 = dW[:,i+1]
                        ddp1 = ddW[:,i+1]


                        U1 = (p1-p0-dp0*T-0.5*ddp0*TT)/TTT
                        U2 = (dp1-dp0-ddp0*T)/TT
                        U3 = (ddp1-ddp0)/T

                        #A=p0
                        #B=dp0
                        #C=0.5*ddp0
                        #D = (10*U1 - 4*U2 + 0.5*U3)
                        #E = (-15*U1 + 7*U2 - U3)/T
                        #F = (6*U1 - 3*U2 + 0.5*U3)/TT

                        #A = p0
                        #B = dp0
                        #C = (3*(p1-p0-dp0*T)-(dp1-dp0)*T)/TT
                        #D = ( (dp1-dp0)*T - 2*(p1-p0-dp0*T))/TTT

                        A = p0
                        B = (p1-p0)/T
                        C = 0.5*(dp1-dp0)/TT

                        def g(dt,A,B,C,D,E,F):
                                return A + B*dt + C*dt**2 + D*dt**3 + E*dt**4 + F*dt**5
                        def dg(dt,A,B,C,D,E,F):
                                return B + 2*C*dt + 3*D*dt**2 + 4*E*dt**3 + 5*F*dt**4
                        def ddg(dt,A,B,C,D,E,F):
                                return 2*C + 6*D*dt + 12*E*dt**2 + 20*F*dt**3

                        dd = np.linalg.norm(W[:,i+1]-g(T,A,B,C,D,E,F))
                        if dd > self.TRAJECTORY_ACCURACY_REQUIRED:
                                np.set_printoptions(precision=4)
                                print "##############################"
                                print "ERROR: Waypoints and C2-approximation do not match up"
                                print "##############################"
                                #print W[:,i],g(0,A,B,C,D,E,F)
                                print "W    : ",W[:,i+1]
                                print "C2   : ",g(duration,A,B,C,D,E,F)
                                print "##############################"
                                #print dW[:,i],dg(0,A,B,C,D,E,F)
                                print "dW   : ",dW[:,i+1]
                                print "dC2  : ",dg(duration,A,B,C,D,E,F)
                                print "##############################"
                                #print ddW[:,i],ddg(0,A,B,C,D,E,F)
                                print "ddW  : ",ddW[:,i+1]
                                print "ddC2 : ",ddg(duration,A,B,C,D,E,F)
                                tvec = np.linspace(0,duration,100)
                                G = np.array(map(lambda t: g(t,A,B,C,D,E,F), tvec)).T
                                dG = np.array(map(lambda t: dg(t,A,B,C,D,E,F), tvec)).T
                                ddG = np.array(map(lambda t: ddg(t,A,B,C,D,E,F), tvec)).T
                                subplot(3,1,1)
                                print G.shape
                                ms = 10
                                plot(0,W[0,i],'ok',markersize=ms)
                                plot(duration,W[0,i+1],'ok',markersize=ms)
                                plot(tvec,G[0,:],'-r',linewidth=2)

                                plot(0,W[1,i],'ok',markersize=ms)
                                plot(duration,W[1,i+1],'ok',markersize=ms)
                                plot(tvec,G[1,:],'-b',linewidth=2)
                                subplot(3,1,2)
                                plot(0,dW[0,i],'ok',markersize=ms)
                                plot(duration,dW[0,i+1],'ok',markersize=ms)
                                plot(tvec,dG[0,:],'-r',linewidth=2)
                                subplot(3,1,3)
                                plot(0,ddW[0,i],'ok',markersize=ms)
                                plot(duration,ddW[0,i+1],'ok',markersize=ms)
                                plot(tvec,ddG[0,:],'-r',linewidth=2)
                                plt.show()
                                sys.exit(0)
                                sys.exit(0)


                        #tvec = np.linspace(0,duration,100)
                        #G = np.array(map(lambda t: g(t,A,B,C,D,E,F), tvec)).T
                        #dG = np.array(map(lambda t: dg(t,A,B,C,D,E,F), tvec)).T
                        #ddG = np.array(map(lambda t: ddg(t,A,B,C,D,E,F), tvec)).T
                        #subplot(3,1,1)
                        #print G.shape
                        #ms = 10
                        #plot(0,W[0,i],'ok',markersize=ms)
                        #plot(duration,W[0,i+1],'ok',markersize=ms)
                        #plot(tvec,G[0,:],'-r',linewidth=2)
                        #subplot(3,1,2)
                        #plot(0,dW[0,i],'ok',markersize=ms)
                        #plot(duration,dW[0,i+1],'ok',markersize=ms)
                        #plot(tvec,dG[0,:],'-r',linewidth=2)
                        #subplot(3,1,3)
                        #plot(0,ddW[0,i],'ok',markersize=ms)
                        #plot(duration,ddW[0,i+1],'ok',markersize=ms)
                        #plot(tvec,ddG[0,:],'-r',linewidth=2)
                        #plt.show()
                        #sys.exit(0)

                        for j in range(Ndim):
                                trajectorystring += "\n"
                                trajectorystring += string.join([str(A[j]),str(B[j]),str(C[j]),str(D[j]),str(E[j]),str(F[j])])

                return [trajectorystring, durationVector]

