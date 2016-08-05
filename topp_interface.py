import string,time
import os
from pylab import *
import numpy as np
from openravepy import *
from util import *
import parameters_dynamical_system as params
import TOPP
from TOPP import Utilities
from TOPP import TOPPbindings
from TOPP import TOPPpy
from TOPP import Trajectory
from TOPP import TOPPopenravepy

### TOPP ERROR CODES
#       TOPP_UNSPEC = 0
#       TOPP_OK = 1
#       TOPP_CANNOT_PREPROCESS = 2
#       TOPP_SHORT_TRAJ = 3
#       TOPP_MVC_HIT_ZERO = 4
#       TOPP_CLC_ERROR = 5
#       TOPP_SDBEGMIN_TOO_HIGH = 6
#       TOPP_SDENDMIN_TOO_HIGH = 7
#       TOPP_FWD_HIT_ZERO = 8
#       TOPP_BWD_HIT_ZERO = 9
#       TOPP_FWD_FAIL = 10
#       TOPP_BWD_FAIL = 11
#       
#       MESSAGES = {
#           TOPP_UNSPEC: "unspecified error",
#           TOPP_OK: "everything OK",
#           TOPP_CANNOT_PREPROCESS: "cannot preprocess trajectory",
#           TOPP_SHORT_TRAJ: "trajectory too short",
#           TOPP_MVC_HIT_ZERO: "MVC hit the sd=0 axis",
#           TOPP_CLC_ERROR: "some CLC error",
#           TOPP_SDBEGMIN_TOO_HIGH: "sdbegmin is too high",
#           TOPP_SDENDMIN_TOO_HIGH: "sdendmin is too high",
#           TOPP_FWD_HIT_ZERO: "forward integration hit the sd=0 axis",
#           TOPP_BWD_HIT_ZERO: "backward integration hit the sd=0 axis",
#           TOPP_FWD_FAIL: "forward integration failed",
#           TOPP_BWD_FAIL: "backward integration failed"
#       }

class TOPPInterface():
        fs_title = 45
        fs_labels = 30
        fs_ticks = 24
        fs_legend = 31

        discrtimestep = 5*1e-3
        TRAJECTORY_ACCURACY_REQUIRED = 1e-10
        traj0 = []
        env_ptr = None
        traj1 = []
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
        zeroForce_ = False

        def initializeFromSpecifications(self, W):
                self.filename = 'images/topp_'+str(params.FILENAME)+'.png'

                self.Ndim = W.shape[0]
                self.Nwaypoints = W.shape[1]
                self.traj0 = Utilities.InterpolateViapoints(W)
                self.length = self.traj0.duration

                dendpoint = np.linalg.norm(self.traj0.Eval(self.length)-W[:,-1])
                if dendpoint > self.TRAJECTORY_ACCURACY_REQUIRED:
                        print "###############"
                        print "TOPP INTERFACE ERROR"
                        print "path duration:",self.length
                        print "###############"
                        print "FINAL POINT on piecewise C^2 traj:",self.traj0.Eval(self.length)
                        print "FINAL WAYPOINT                   :",W[:,-1]
                        print "distance between points:",dendpoint
                        print "###############"
                        sys.exit(1)

                [a,b,c] = self.GetABC(self.traj0, self.discrtimestep)
                vmax = 1e5*np.ones(self.Ndim)
                self.topp_inst = TOPP.QuadraticConstraints(self.traj0, self.discrtimestep, vmax, list(a), list(b), list(c))
                self.waypoints = W

        def __init__(self, W, env, zeroForce=False):

                self.zeroForce_ = zeroForce
                self.env_ptr = env
                self.initializeFromSpecifications(W)

        def GetABC(self, traj_in, discrtimestep):
                ndiscrsteps = int((traj_in.duration + 1e-10) / discrtimestep) + 1
                Adim = params.amin.shape[0]

                q = np.zeros((self.Ndim,ndiscrsteps))
                qs = np.zeros((self.Ndim,ndiscrsteps))
                qss = np.zeros((self.Ndim,ndiscrsteps))
                a = np.zeros((ndiscrsteps, 2*Adim))
                b = np.zeros((ndiscrsteps, 2*Adim))
                c = np.zeros((ndiscrsteps, 2*Adim))


                for i in range(0,ndiscrsteps):
                        duration = i*discrtimestep

                        q[:,i] = traj_in.Eval(duration)
                        qs[:,i] = traj_in.Evald(duration)
                        qss[:,i] = traj_in.Evaldd(duration)

                        ### G*qdd + h <= 0

                        if self.zeroForce_:
                                Fi = np.zeros((q[:,i].shape))
                        else:
                                Fi = params.waypoint_to_force(self.env_ptr, q[:,i])

                        Ri = params.GetControlMatrixAtWaypoint(q[:,i])
                        [G,h] = params.GetControlConstraintMatricesFromControl(Ri, Fi)
                        a[i,:] = np.dot(G,qs[:,i]).flatten()
                        b[i,:] = np.dot(G,qss[:,i]).flatten()
                        c[i,:] = h

                ### make them accesible to the outside world
                self.a = a
                self.b = b
                self.c = c

                self.q = q
                self.qs = qs
                self.qss = qss

                return [a,b,c]


        def ReparameterizeTrajectory(self):
                x = self.topp_inst.solver
                #x.integrationtimestep = self.discrtimestep
                #x.reparamtimestep = 0.0001
                ret = x.RunComputeProfiles(0.0,0.0)
                if ret == 1:
                        x.ReparameterizeTrajectory()
                        x.WriteResultTrajectory()
                        try:
                                self.traj1 = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
                        except Exception as e:
                                print "Exception in ReparameterizeTrajectory:",e
                                print "trajectorystring is :",x.restrajectorystring
                                print "discrtimestep=",self.discrtimestep
                                #print "trajectorystring=\"\"\"",self.traj0,"\"\"\""
                                #PrintNumpy("a=",self.a)
                                #PrintNumpy("b=",self.b)
                                #PrintNumpy("c=",self.c)
                                sys.exit(1)
                        return True
                else:
                        return False

        def CriticalTimeToWaypoint(self, t):
                kctr = 0

                if t>self.traj0.duration:
                        return self.Nwaypoints

                dt = self.traj0.chunkcumulateddurationslist[kctr]
                while dt < t:
                        kctr+=1
                        dt = self.traj0.chunkcumulateddurationslist[kctr]
                return kctr
                
        def getCriticalPoint(self):
                x = self.topp_inst.solver
                x.integrationtimestep = self.discrtimestep
                #x.reparamtimestep = 0.01
                #x.extrareps = 10

                self.critical_point = self.Nwaypoints
                print self.Nwaypoints
                try:
                        ret = x.RunComputeProfiles(0.0,0.0)
                        if ret == 4:
                                CP_dt = x.GetCriticalPoint()
                                print "CPdt",CP_dt
                                t_CP = CP_dt * x.integrationtimestep
                                self.critical_point = self.CriticalTimeToWaypoint(t_CP)
                                print "TOPP: MVC hit zero | CP:",
                                print self.critical_point,"/",self.Nwaypoints,
                                print "(t=",t_CP,"/",self.traj0.duration,")"
                                #print self.getSpeedIntervalAtPoint(4)
                                #self.ReparameterizeTrajectory()
                                #sys.exit(0)
                                return self.critical_point
                        if ret == 1: 
                                print "TOPP: success"
                                return self.Nwaypoints
                        if ret == 0:
                                CP_dt = x.GetCriticalPoint()
                                t_CP = CP_dt * x.integrationtimestep
                                print CP_dt,t_CP
                                self.critical_point = self.CriticalTimeToWaypoint(t_CP)
                                print "TOPP: unspecified error | CP:",
                                print self.critical_point,"/",self.Nwaypoints,
                                print "(t=",t_CP,"/",self.traj0.duration,")"
                                return self.critical_point
                        else:
                                print "TOPP: error code ",ret,"(not handled yet)"
                                sys.exit(0)

                except Exception as e:
                        print "TOPP EXCEPTION: ",e
                        print self.traj0.duration
                        sys.exit(0)

        def getSpeedIntervalAtPoint(self, N):
                if (N < 0) or (N > self.Nwaypoints):
                        print "Critical Point",N,"is out of bounds of trajectory"
                        sys.exit(1)
                if N == 0:
                        return [0,0]

                Wn = self.waypoints[:,0:N]
                trajN = Utilities.InterpolateViapoints(Wn)
                [a,b,c] = self.GetABC(trajN, self.discrtimestep)
                vmax = 1e5*np.ones(self.Ndim)

                topp_inst = TOPP.QuadraticConstraints(trajN, self.discrtimestep, vmax, list(a), list(b), list(c))
                x = topp_inst.solver
                ret = x.RunVIP(0.0, 0.0)

                if ret != 1:
                        semin = 0.0
                        semax = 0.0
                else:
                        semin = x.sdendmin
                        semax = x.sdendmax
                        if semin<=1e-5:
                                ret2 = x.RunComputeProfiles(0.0, semin)
                                if ret2 != 1:
                                        #print N,ret,ret2,semin,semax
                                        #np.savetxt('topp/a2',a)
                                        #np.savetxt('topp/b2',b)
                                        #np.savetxt('topp/c2',c)

                                        #tstr = "trajectorystring=\"\"\" %s \"\"\""%(str(trajN))
                                        #with open("topp/traj2", "w") as fh:
                                        #        fh.write("%s" % str(trajN))
                                        self.SaveToFile('vip')

                                        sys.exit(0)


                return [semin, semax]

        def SaveToFile(self, postfix):
                prefix='topp/'
                np.savetxt(prefix+postfix+'_a',self.a)
                np.savetxt(prefix+postfix+'_b',self.b)
                np.savetxt(prefix+postfix+'_c',self.c)
                np.savetxt(prefix+postfix+'_W',self.waypoints)

                #tstr = "trajectorystring=\"\"\" %s \"\"\""%(str(self.traj0))

                with open(prefix+postfix+"_traj", "w") as fh:
                        fh.write("%s" % str(self.traj0))

                print "Saved TOPP to files",prefix+postfix+"*"

        @classmethod
        def LoadFromFile(cls, prepostfix, env):
                self.a = np.loadtxt(prepostfix+'_a')
                self.b = np.loadtxt(prepostfix+'_b')
                self.c = np.loadtxt(prepostfix+'_c')
                self.waypoints = np.loadtxt(prepostfix+'_W')
                with open(prepostfix+"_traj", "r") as fh:
                        self.traj0 = "%s" % fh.read()
                return cls(self.waypoints, env)



        def PlotPrettifiedAxes(self, ax):
                plt.axvspan(0, 1.0, facecolor='k', alpha=0.1)
                plt.axvspan(1.0, self.traj1.duration, facecolor='g', alpha=0.1)
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=self.fs_legend)
                ax.tick_params(axis='both', which='major', pad=15)
                ax.tick_params(axis='both', which='minor', pad=15)
                ax.yaxis.get_offset_text().set_size(self.fs_labels)
                ax.xaxis.get_offset_text().set_size(self.fs_labels)

                for tick in ax.xaxis.get_major_ticks():
                        tick.label.set_fontsize(self.fs_ticks) 
                for tick in ax.yaxis.get_major_ticks():
                        tick.label.set_fontsize(self.fs_ticks) 
                ax.xaxis.labelpad = 20
                ax.yaxis.labelpad = 20

        def PlotVerticalLineOverSubplots(self, x, ax1, ax2, ax3, ax4):
                ax1.axvline(x=x,ymin=-1.2,ymax=1  ,c="k",ls='--',linewidth=3,zorder=0, clip_on=False)
                ax2.axvline(x=x,ymin=-1.2,ymax=1  ,c="k",ls='--',linewidth=3,zorder=0, clip_on=False)
                ax3.axvline(x=x,ymin=-1.2,ymax=1  ,c="k",ls='--',linewidth=3,zorder=0, clip_on=False)
                ax4.axvline(x=x,ymin=-0.2,ymax=1.2,c="k",ls='--',linewidth=3,zorder=0, clip_on=False)

        def PlotTrajectory(self, env=None):
                lw = 4
                fs = self.fs_legend
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
                offset_a_coordinate = 0.1

                dt = float(self.discrtimestep)
                while dt > self.traj1.duration:
                        dt/=2.0
                Npts = int(self.traj1.duration/dt)
                print "dt",dt,"duration:",self.traj1.duration

                tvect = np.linspace(0,self.traj1.duration, Npts)
                qvect = np.array([self.traj1.Eval(t) for t in tvect])
                qdvect = np.array([self.traj1.Evald(t) for t in tvect])
                qddvect = np.array([self.traj1.Evaldd(t) for t in tvect])

                Noripts = int(self.traj0.duration/dt)
                tvect_ori= np.linspace(0,self.traj0.duration, Noripts)
                qori_vect = np.array([self.traj0.Eval(t) for t in tvect_ori])
                qori_dvect = np.array([self.traj0.Evald(t) for t in tvect_ori])
                qori_ddvect = np.array([self.traj0.Evaldd(t) for t in tvect_ori])

                #################################
                Adim = params.amin.shape[0]
                a = np.zeros((Adim, Npts))
                if env is not None:
                        [R,amin,amax] = path.getControlMatrix(qvect.T)
                        amin = params.amin
                        amax = params.amax
                        Ndim = qvect.shape[0]
                        Nwaypoints = qvect.shape[1]
                        assert(Ndim==4)
                        R = params.ControlPerWaypoint(qvect.T, Ndim, Nwaypoints)
                        F = params.waypoints_to_force(env, qvect.T)
                        #F = path.get_forces_at_waypoints(qvect.T, env)
                        for i in range(0,Npts):
                                Ri = R[:,:,i]
                                Fi = F[:,i]
                                qdd = qddvect[i,:]
                                Rinv = np.linalg.pinv(Ri)
                                I = np.identity(self.Ndim)
                                a[:,i] = np.dot(Rinv,np.dot(I,qdd)-Fi)
                                for k in range(0,Adim):
                                        if a[k,i] > amax[k]:
                                                a[k,i] = amax[k]
                                        if a[k,i] < amin[k]:
                                                a[k,i] = amin[k]
                                #if (a[:,i]>amax).any():
                                #        print "i",i,"/",Npts
                                #        print "a[:,",i,",]=",a[:,i],">",amax
                                #        for k in range(0,Adim):
                                #                print amin[k],"<=",a[k,i],"<=",amax[k]
                                #        print "R",Ri,Rinv
                                #        print "W",W
                                #        print "Fi",Fi
                                #        print "I",I
                                #        print "qdd",qdd
                                #        sys.exit(0)
                #################################

                twvect = np.linspace(0,np.sum(self.durationVector), self.Nwaypoints)
                tavect = np.linspace(0,self.traj1.duration, self.Nwaypoints)
                
                fig=figure(facecolor='white')

                ax1 = subplot(4,1,1)

                #plot(twvect, self.W_[0,:], '--', color = color_x_coordinate, linewidth = lw)
                #plot(twvect, self.W_[1,:], '--', color = color_y_coordinate, linewidth = lw)
                #plot(twvect, self.W_[3,:], '--', color = color_t_coordinate, linewidth = lw)

                plot(tvect, qvect[:,0], color = color_x_coordinate, linewidth = lw, label = "$x$")
                plot(tvect, qvect[:,1], color = color_y_coordinate, linewidth = lw, label = "$y$")
                #plot(tvect, qvect[:,2], color = color_z_coordinate, linewidth = lw, label = "$z$")
                plot(tvect, qvect[:,3], color = color_t_coordinate, linewidth = lw, label = "$\\theta$")
                #plot(tvect, qdvect, f, linewidth=2)
                #plot(tvect, qddvect, f, linewidth=2)
                #title('TOPP-Profile '+str(params.system_name), fontsize=self.fs_title, y=1.1)
                ylabel('Position \n$\\left(m\\right)$,$\\left(rad\\right)$', fontsize=self.fs_labels)
                ax1.set_xticklabels(())
                self.PlotPrettifiedAxes(ax1)

                ax2 = subplot(4,1,2)
                ylabel('Velocity\n$\\left(\\frac{m}{s}\\right)$, $\\left(\\frac{rad}{s}\\right)$', fontsize=self.fs_labels)
                #plot(twvect, 0*self.dW_[0,:], '-', color=black)
                #plot(twvect, self.dW_[0,:], '--', color = color_x_coordinate, linewidth = lw)
                #plot(twvect, self.dW_[1,:], '--', color = color_y_coordinate, linewidth = lw)
                #plot(twvect, self.dW_[3,:], '--', color = color_t_coordinate, linewidth = lw)

                plot(tvect, qdvect[:,0], color = color_x_coordinate, linewidth = lw, label = "$\dot x$")
                plot(tvect, qdvect[:,1], color = color_y_coordinate, linewidth = lw, label = "$\dot y$")
                #plot(tvect, qdvect[:,2], color = color_z_coordinate, linewidth = lw, label = "$\dot z$")
                plot(tvect, qdvect[:,3], color = color_t_coordinate, linewidth = lw, label = "$\dot \\theta$")
                ax2.set_xticklabels(())
                self.PlotPrettifiedAxes(ax2)

                ax3 = subplot(4,1,3)
                ylabel('Acceleration\n$\\left(\\frac{m}{s^2}\\right)$,$\\left(\\frac{rad}{s^2}\\right)$', fontsize=self.fs_labels)
                #plot(twvect, 0*self.dW_[0,:], '-', color=black)
                #plot(tvect_ori, qori_ddvect[:,0], '--', color = color_x_coordinate, linewidth = lw)
                #plot(tvect_ori, qori_ddvect[:,1], '--', color = color_y_coordinate, linewidth = lw)
                #plot(tvect_ori, qori_ddvect[:,3], '--', color = color_t_coordinate, linewidth = lw)

                #plot(twvect, self.dW_[1,:], '--', color = color_y_coordinate, linewidth = lw)
                #plot(twvect, self.dW_[3,:], '--', color = color_t_coordinate, linewidth = lw)

                plot(tvect, qddvect[:,0], color = color_x_coordinate, linewidth = lw, label = "$\ddot{x}$")
                plot(tvect, qddvect[:,1], color = color_y_coordinate, linewidth = lw, label = "$\ddot{y}$")
                #plot(tvect, qddvect[:,2], color = color_z_coordinate, linewidth = lw, label = "$\ddot{z}$")
                plot(tvect, qddvect[:,3], color = color_t_coordinate, linewidth = lw, label = "$\ddot{\\theta}$")
                #plot(tvect, F[0,:], color = color_fx_coordinate, linewidth = lw, label = "$F_{x}$")
                #plot(tvect, F[1,:], color = color_fy_coordinate, linewidth = lw, label = "$F_{y}$")
                ax3.set_xticklabels(())
                self.PlotPrettifiedAxes(ax3)

                if env is not None:
                        ax4 = subplot(4,1,4)
                        ylabel('Control\n$\\left(\\frac{m}{s^2}\\right)$,$\\left(\\frac{rad}{s^2}\\right)$', fontsize=self.fs_labels)
                        #plot(twvect, 0*self.dW_[0,:], '-', color=black)
                        plot(tvect, a[0,:], color = color_a1_coordinate, linewidth = lw, label = "${a_1}(Thruster)$")
                        plot(tvect, a[1,:], color = color_a2_coordinate, linewidth = lw, label = "${a_2}(Lie Bracket)$")
                        plot(tvect, a[2,:], color = color_a3_coordinate, linewidth = lw, label = "${a_3}(Steer)$")

                        plot(tvect, np.repeat(2*pi*amin[2]-offset_a_coordinate,tvect.size), lw = limit_lw, ls = limit_ls, color = color_a3_coordinate)
                        plot(tvect, np.repeat(2*pi*amax[2]+offset_a_coordinate,tvect.size), lw = limit_lw, ls = limit_ls, color = color_a3_coordinate)

                        plot(tvect, np.repeat(amin[1]-offset_a_coordinate,tvect.size), lw = limit_lw, ls = limit_ls, color = color_a2_coordinate)
                        plot(tvect, np.repeat(amax[1]+offset_a_coordinate,tvect.size), lw = limit_lw, ls = limit_ls, color = color_a2_coordinate)

                        plot(tvect, np.repeat(amin[0]-offset_a_coordinate,tvect.size), lw = limit_lw, ls = limit_ls, color = color_a1_coordinate)
                        plot(tvect, np.repeat(amax[0]+offset_a_coordinate,tvect.size), lw = limit_lw, ls = limit_ls, color = color_a1_coordinate)
                        xlabel('Time ($s$)',fontsize=self.fs_labels)
                        self.PlotPrettifiedAxes(ax4)
                        self.PlotVerticalLineOverSubplots(self.traj1.duration, ax1, ax2, ax3, ax4)
                        self.PlotVerticalLineOverSubplots(1.0, ax1, ax2, ax3, ax4)
                else:
                        self.PlotVerticalLineOverSubplots(self.traj1.duration, ax1, ax2, ax3)
                        self.PlotVerticalLineOverSubplots(1.0, ax1, ax2, ax3)

                #plt.gca().tight_layout()
                plt.gcf().subplots_adjust(bottom=0.15,right=0.8)
                plt.show()

                fig.savefig(self.filename,format='svg', dpi=1200)

                svgname = self.filename
                pdfname = os.path.splitext(self.filename)[0]+'.pdf'

                syscmd = 'mogrify -format pdf -trim '+svgname
                os.system(syscmd)
                syscmd = 'pdfcrop '+pdfname+' '+pdfname
                os.system(syscmd)
                os.system('cp /home/`whoami`/git/openrave/sandbox/WPI/images/*.pdf /home/`whoami`/git/papers/images/simulation/')


