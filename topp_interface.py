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
        DURATION_DISCRETIZATION = 0.001

        TRAJECTORY_ACCURACY_REQUIRED = 1e-5
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
        path_length = -1
        F_ = []
        R_ = []
        amin_ = []
        amax_ = []
        W_ = []
        dW_ = []
        ddW_ = []
        trajectoryclass_ = []

        def initializeFromSpecifications(self, durationVector_in, trajectorystring, F, R, amin, amax, W, dW):
                self.Ndim = W.shape[0]
                self.Nwaypoints = W.shape[1]

                self.trajstr = trajectorystring
                self.durationVector = durationVector_in

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

        def __init__(self, trajectoryclass, durationVector, trajectorystring, F, R, amin, amax, W, dW):
                self.trajectoryclass_ = trajectoryclass
                self.path_length = trajectoryclass.get_length()
                self.F_ = F
                self.R_ = R
                self.amin_ = amin
                self.amax_ = amax
                self.W_ = W
                self.dW_ = dW
                self.initializeFromSpecifications(durationVector, trajectorystring, F, R, amin, amax, W, dW)

        def getSpeedIntervalAtCriticalPoint(self, N, Subtraj_dvec, Subtraj_str):
                if (N < 0) or (N > self.Nwaypoints):
                        print "Critical Point",N,"is out of bounds of trajectory"
                        sys.exit(1)
                if N == 0:
                        return [0,0]

                Fc = self.F_[:,0:N]
                Rc = self.R_[:,:,0:N]
                Wc = self.W_[:,0:N]
                dWc = self.dW_[:,0:N]
                amin = self.amin_
                amax = self.amax_
                self.initializeFromSpecifications(Subtraj_dvec, Subtraj_str, Fc, Rc, amin, amax, Wc, dWc)
                x = self.topp_inst.solver
                ### AVP(sbmin, sbmax)
                #[semin,semax] = self.topp_inst.AVP(0.0, 0.0)
                return_code = x.RunVIP(0.0, 0.0)
                if return_code != 1:
                        print "TOPP Error:", return_code
                        print "waypoint: ",N,"/",self.Nwaypoints
                        self.critical_point = x.GetCriticalPoint()
                        self.critical_point_value = x.GetCriticalPointValue()
                        print "CRITICAL POINT:",self.critical_point,self.critical_point_value
                        sys.exit(1)

                semin = x.sdendmin
                semax = x.sdendmax
                return [semin, semax]

        def PlotPrettifiedAxes(self, ax, fs):
                plt.axvspan(0, 1.0, facecolor='k', alpha=0.1)
                plt.axvspan(1.0, self.traj0.duration, facecolor='g', alpha=0.1)
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=fs)

        def PlotVerticalLineOverSubplots(self, x, ax1, ax2, ax3, ax4):
                ax1.axvline(x=x,ymin=-1.2,ymax=1  ,c="k",ls='--',linewidth=3,zorder=0, clip_on=False)
                ax2.axvline(x=x,ymin=-1.2,ymax=1  ,c="k",ls='--',linewidth=3,zorder=0, clip_on=False)
                ax3.axvline(x=x,ymin=-1.2,ymax=1  ,c="k",ls='--',linewidth=3,zorder=0, clip_on=False)
                ax4.axvline(x=x,ymin=-0.2,ymax=1.2,c="k",ls='--',linewidth=3,zorder=0, clip_on=False)

        def PlotTrajectory(self, env=None):
                lw = 4
                fs = 22
                limit_ls = ':'
                limit_lw = 4
                color_x_coordinate = 'b'
                color_y_coordinate = 'r'
                color_fx_coordinate = (0.7,0.7,0.9)
                color_fy_coordinate = (0.9,0.7,0.7)
                color_z_coordinate = 'g'
                color_t_coordinate = 'm'
                color_a1_coordinate = (1.0,0.0,0.0)
                color_a2_coordinate = (1.0,0.5,0.0)
                color_a3_coordinate = (1.0,0.8,0.6)
                offset_a_coordinate = 0.3

                #################################
                        #dt = self.DISCRETIZATION_TIME_STEP
                        #L = self.get_length()
                        #N = int(L/dt)
                        #Npts = int(self.path_length/dt)
                dt = float(self.DURATION_DISCRETIZATION)
                Npts = int(self.traj0.duration/dt)
                tvect = np.linspace(0,self.traj0.duration, Npts)
                qvect = np.array([self.traj0.Eval(t) for t in tvect])
                qdvect = np.array([self.traj0.Evald(t) for t in tvect])
                qddvect = np.array([self.traj0.Evaldd(t) for t in tvect])

                #################################
                path = self.trajectoryclass_
                Adim = 3
                a = np.zeros((Adim, Npts))
                if env is not None:
                        [R,amin,amax] = path.getControlMatrix(qvect.T)
                        F = path.get_forces_at_waypoints(qvect.T, env)
                        for i in range(0,Npts):
                                Ri = R[:,:,i]
                                Fi = F[:,i]
                                qdd = qddvect[i,:]
                                a[:,i] = np.dot(Ri.T,qdd.T-Fi)

                #################################
                twvect = np.linspace(0,np.sum(self.durationVector), self.Nwaypoints)
                tavect = np.linspace(0,self.traj0.duration, self.Nwaypoints)
                
                fig=figure(facecolor='white')
                ax1 = subplot(4,1,1)
                plot(twvect, self.W_[0,:], '--', color = color_x_coordinate, linewidth = lw)
                plot(twvect, self.W_[1,:], '--', color = color_y_coordinate, linewidth = lw)
                plot(twvect, self.W_[2,:], '--', color = color_z_coordinate, linewidth = lw)
                plot(twvect, self.W_[3,:], '--', color = color_t_coordinate, linewidth = lw)
                plot(tvect, qvect[:,0], color = color_x_coordinate, linewidth = lw, label = "$x$")
                plot(tvect, qvect[:,1], color = color_y_coordinate, linewidth = lw, label = "$y$")
                plot(tvect, qvect[:,2], color = color_z_coordinate, linewidth = lw, label = "$z$")
                plot(tvect, qvect[:,3], color = color_t_coordinate, linewidth = lw, label = "$\\theta$")
                #plot(tvect, qdvect, f, linewidth=2)
                #plot(tvect, qddvect, f, linewidth=2)
                title('TOPP-Profile Trajectory', fontsize=fs)
                ylabel('Position $(m)$', fontsize=fs)
                ax1.set_xticklabels(())
                self.PlotPrettifiedAxes(ax1, fs)

                ax2 = subplot(4,1,2)
                ylabel('Velocity $\\left(\\frac{m}{s}\\right)$', fontsize=fs)
                plot(twvect, self.dW_[0,:], '--', color = color_x_coordinate, linewidth = lw)
                plot(twvect, self.dW_[1,:], '--', color = color_y_coordinate, linewidth = lw)
                plot(twvect, self.dW_[2,:], '--', color = color_z_coordinate, linewidth = lw)
                plot(twvect, self.dW_[3,:], '--', color = color_t_coordinate, linewidth = lw)
                plot(tvect, qdvect[:,0], color = color_x_coordinate, linewidth = lw, label = "$\dot x$")
                plot(tvect, qdvect[:,1], color = color_y_coordinate, linewidth = lw, label = "$\dot y$")
                plot(tvect, qdvect[:,2], color = color_z_coordinate, linewidth = lw, label = "$\dot z$")
                plot(tvect, qdvect[:,3], color = color_t_coordinate, linewidth = lw, label = "$\dot \\theta$")
                ax2.set_xticklabels(())
                self.PlotPrettifiedAxes(ax2, fs)

                ax3 = subplot(4,1,3)
                ylabel('Acceleration $\\left(\\frac{m^2}{s^2}\\right)$', fontsize=fs)
                plot(tvect, qddvect[:,0], color = color_x_coordinate, linewidth = lw, label = "$\ddot{x}$")
                plot(tvect, qddvect[:,1], color = color_y_coordinate, linewidth = lw, label = "$\ddot{y}$")
                plot(tvect, qddvect[:,2], color = color_z_coordinate, linewidth = lw, label = "$\ddot{z}$")
                plot(tvect, qddvect[:,3], color = color_t_coordinate, linewidth = lw, label = "$\ddot{\\theta}$")
                plot(tvect, F[0,:], color = color_fx_coordinate, linewidth = lw, label = "$F_{x}$")
                plot(tvect, F[1,:], color = color_fy_coordinate, linewidth = lw, label = "$F_{y}$")
                ax3.set_xticklabels(())
                self.PlotPrettifiedAxes(ax3, fs)

                if env is not None:
                        ax4 = subplot(4,1,4)
                        ylabel('Control $\\left(\\frac{m^2}{s^2}\\right)$', fontsize=fs)
                        plot(tvect, a[0,:], color = color_a1_coordinate, linewidth = lw, label = "${a_1}(Thruster)$")
                        plot(tvect, a[1,:], color = color_a2_coordinate, linewidth = lw, label = "${a_2}(Lie Bracket)$")
                        plot(tvect, a[2,:], color = color_a3_coordinate, linewidth = lw, label = "${a_3}(Steer)$")

                        plot(tvect, np.repeat(amin[2]-offset_a_coordinate,tvect.size), lw = limit_lw, ls = limit_ls, color = color_a3_coordinate)
                        plot(tvect, np.repeat(amax[2]+offset_a_coordinate,tvect.size), lw = limit_lw, ls = limit_ls, color = color_a3_coordinate)
                        plot(tvect, np.repeat(amin[1]-offset_a_coordinate,tvect.size), lw = limit_lw, ls = limit_ls, color = color_a2_coordinate)
                        plot(tvect, np.repeat(amax[1]+offset_a_coordinate,tvect.size), lw = limit_lw, ls = limit_ls, color = color_a2_coordinate)
                        plot(tvect, np.repeat(amin[0]-offset_a_coordinate,tvect.size), lw = limit_lw, ls = limit_ls, color = color_a1_coordinate)
                        plot(tvect, np.repeat(amax[0]+offset_a_coordinate,tvect.size), lw = limit_lw, ls = limit_ls, color = color_a1_coordinate)
                        xlabel('Time ($s$)',fontsize=fs)
                        self.PlotPrettifiedAxes(ax4, fs)
                        self.PlotVerticalLineOverSubplots(self.traj0.duration, ax1, ax2, ax3, ax4)
                        self.PlotVerticalLineOverSubplots(1.0, ax1, ax2, ax3, ax4)
                else:
                        self.PlotVerticalLineOverSubplots(self.traj0.duration, ax1, ax2, ax3)
                        self.PlotVerticalLineOverSubplots(1.0, ax1, ax2, ax3)

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
                                self.critical_point = x.GetCriticalPoint()
                                self.critical_point_value = x.GetCriticalPointValue()
                                #print "TOPP critical pt:",self.critical_point,self.critical_point_value
                                return self.critical_point
                        if ret == 1: ##TOPP_OK
                                #print "TOPP: success"
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
                                return self.Nwaypoints
                        if ret== 0:
                                #M = 10
                                #np.set_printoptions(precision=3)
                                #tvect = np.linspace(0,self.traj0.duration, self.traj0.duration/dt)
                                #qvect = np.array([self.traj0.Eval(t) for t in tvect])
                                #qdvect = np.array([self.traj0.Evald(t) for t in tvect])
                                #qddvect = np.array([self.traj0.Evaldd(t) for t in tvect])
                                print "TOPP: unspecified error"
                                #print "W=",repr(self.W_[0:2,0:M])
                                #print "q=",repr(qvect[0:M,0:2].T)
                                #print "dW=",repr(self.dW_[0:2,0:M])
                                #print "dq=",repr(qdvect[0:M,0:2].T)
                                #self.critical_point = x.GetCriticalPoint()
                                #self.critical_point_value = x.GetCriticalPointValue()
                                #print self.critical_point,self.critical_point_value
                                return -1
                                #sys.exit(0)
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

                        #for j in range(self.Ndim):
                                #print H2[j],"<= q[",j,"]<=",-H1[j]

                        c[i,:] = np.hstack((H1,H2)).flatten()
                #sys.exit(0)

                for i in range(0,self.Nwaypoints):
                        a[i,:] = np.dot(G,qs[:,i]).flatten()
                        b[i,:] = np.dot(G,qss[:,i]).flatten()
                return [a,b,c]
