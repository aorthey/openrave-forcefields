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

class TOPPInterface():
        fs_title = 45
        fs_labels = 35
        fs_ticks = 24
        fs_legend = 31

        #DURATION_DISCRETIZATION = 0.0001
        #DURATION_DISCRETIZATION = 1
        DURATION_DISCRETIZATION = 0.01

        TRAJECTORY_ACCURACY_REQUIRED = 1e-1
        traj0 = []
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
        F_ = []
        R_ = []
        amin_ = []
        amax_ = []
        W_ = []
        dW_ = []
        ddW_ = []
        trajectoryclass_ = []

        def initializeFromSpecifications(self, durationVector_in, trajectorystring, F, R, amin, amax, W, dW, dt =None):
                self.filename = 'images/topp_'+str(params.FILENAME)+'.png'
                self.Ndim = W.shape[0]
                self.Nwaypoints = W.shape[1]

                [self.trajstr, self.durationVector] = self.computeTrajectoryStringFromWaypoints(F,R,amin,amax,W,dW,dt)
                #self.trajstr = trajectorystring
                #self.durationVector = durationVector_in

                self.traj0 = Trajectory.PiecewisePolynomialTrajectory.FromString(self.trajstr)
                self.length = np.sum(self.durationVector)
                dendpoint = np.linalg.norm(self.traj0.Eval(self.length)-W[:,-1])

                #if dendpoint > self.TRAJECTORY_ACCURACY_REQUIRED:
                #        print "###############"
                #        print "TOPP INTERFACE ERROR"
                #        print "path duration:",self.length
                #        print "###############"
                #        print "FINAL POINT on piecewise C^2 traj:",self.traj0.Eval(self.length)
                #        print "FINAL WAYPOINT                   :",W[:,-1]
                #        print "distance between points:",dendpoint
                #        print "###############"
                #        #sys.exit(1)

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
                #[trajectorystring, durationVector] = self.computeTrajectoryStringFromWaypoints(F,R,amin,amax,W,dW)
                self.initializeFromSpecifications(durationVector, trajectorystring, F, R, amin, amax, W, dW)


        def getSpeedIntervalAtCriticalPoint(self, N, Subtraj_dvec, Subtraj_str, dt=None):
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
                self.initializeFromSpecifications(Subtraj_dvec, Subtraj_str, Fc, Rc, amin, amax, Wc, dWc,dt)
                x = self.topp_inst.solver
                ### AVP(sbmin, sbmax)
                #[semin,semax] = self.topp_inst.AVP(0.0, 0.0)
                return_code = x.RunVIP(0.0, 0.0)
                if return_code != 1:
                        #print "GET SPEED INTERVAL"
                        #print "TOPP Error:", return_code
                        #print "waypoint: ",N,"/",self.Nwaypoints
                        self.critical_point = x.GetCriticalPoint()
                        self.critical_point_value = x.GetCriticalPointValue()
                        #print "CRITICAL POINT:",self.critical_point,self.critical_point_value
                        #plt.plot(Wc[0,:],Wc[1,:],'-or',markersize=5,linewidth=3)
                        #plt.show()
                        #sys.exit(1)
                        semin = 0.0
                        semax = 0.0
                else:
                        semin = x.sdendmin
                        semax = x.sdendmax
                #print Subtraj_dvec
                #print Wc
                #print dWc
                return [semin, semax]

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

                #################################
                        #dt = self.DISCRETIZATION_TIME_STEP
                        #L = self.get_length()
                        #N = int(L/dt)
                        #Npts = int(self.path_length/dt)
                dt = float(self.DURATION_DISCRETIZATION)
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
                path = self.trajectoryclass_
                Adim = params.amin.shape[0]
                a = np.zeros((Adim, Npts))
                if env is not None:
                        [R,amin,amax] = path.getControlMatrix(qvect.T)
                        F = path.get_forces_at_waypoints(qvect.T, env)
                        for i in range(0,Npts):
                                Ri = R[:,:,i]
                                Fi = F[:,i]
                                qdd = qddvect[i,:]
                                Rinv = np.linalg.pinv(Ri)
                                I = np.identity(self.Ndim)
                                a[:,i] = np.dot(Rinv,np.dot(I,qdd)-Fi)
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
                plot(twvect, self.W_[0,:], '--', color = color_x_coordinate, linewidth = lw)
                plot(twvect, self.W_[1,:], '--', color = color_y_coordinate, linewidth = lw)
                #plot(twvect, self.W_[2,:], '--', color = color_z_coordinate, linewidth = lw)
                plot(twvect, self.W_[3,:], '--', color = color_t_coordinate, linewidth = lw)
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
                plot(twvect, 0*self.dW_[0,:], '-', color=black)
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
                plot(twvect, 0*self.dW_[0,:], '-', color=black)
                plot(tvect_ori, qori_ddvect[:,0], '--', color = color_x_coordinate, linewidth = lw)
                plot(tvect_ori, qori_ddvect[:,1], '--', color = color_y_coordinate, linewidth = lw)
                plot(tvect_ori, qori_ddvect[:,3], '--', color = color_t_coordinate, linewidth = lw)
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
                        plot(twvect, 0*self.dW_[0,:], '-', color=black)
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

        def ReparameterizeTrajectory(self):
                x = self.topp_inst.solver
                #x.integrationtimestep = 0.001
                #x.integrationtimestep = 0.001
                #x.reparamtimestep = 0.01
                ret = x.RunComputeProfiles(0.0,0.0)
                if ret == 1:
                        x.ReparameterizeTrajectory()
                        x.WriteResultTrajectory()
                        self.traj1 = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
                        return True
                else:
                        return False

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

        def getCriticalPoint(self):
                x = self.topp_inst.solver
                #x.integrationtimestep = 0.001
                #x.reparamtimestep = 0.01
                #x.extrareps = 10

                self.critical_point = self.Nwaypoints
                try:
                        ret = x.RunComputeProfiles(0.0,0.0)
                        if ret == 4:
                                ### MVC Hit Zero
                                self.critical_point = x.GetCriticalPoint()
                                self.critical_point_value = x.GetCriticalPointValue()
                                print "TOPP critical pt:",self.critical_point,self.critical_point_value
                                return self.critical_point
                        if ret == 1: ##TOPP_OK
                                print "TOPP: success"
                                return self.Nwaypoints
                        if ret== 0:
                                self.critical_point = x.GetCriticalPoint()
                                self.critical_point_value = x.GetCriticalPointValue()
                                print "TOPP: unspecified error"
                                print "TOPP critical pt:",self.critical_point,self.critical_point_value
                                return self.critical_point
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
                Adim = amin.shape[0]
                q = np.zeros((self.Ndim,self.Nwaypoints))
                qs = np.zeros((self.Ndim,self.Nwaypoints))
                qss = np.zeros((self.Ndim,self.Nwaypoints))
                for i in range(0,self.Nwaypoints):
                        duration = np.sum(self.durationVector[0:i])
                        q[:,i] = self.traj0.Eval(duration)
                        qs[:,i] = self.traj0.Evald(duration)
                        qss[:,i] = self.traj0.Evaldd(duration)

                a = np.zeros((self.Nwaypoints, 2*Adim))
                b = np.zeros((self.Nwaypoints, 2*Adim))
                c = np.zeros((self.Nwaypoints, 2*Adim))

                for i in range(0,self.Nwaypoints):
                        ### G*qdd + h <= 0
                        [G,h] = params.GetControlConstraintMatricesFromControl(R[:,:,i], F[:,i])
                        a[i,:] = np.dot(G,qs[:,i]).flatten()
                        b[i,:] = np.dot(G,qss[:,i]).flatten()
                        c[i,:] = h
                        #print a[i,:],b[i,:],c[i,:]
                

                return [a,b,c]

        def GetNearestPolynomial(self, p):
                pass
        
        def computeTrajectoryStringFromWaypoints(self, F, R, amin, amax, W, dW, dt=None):
                DEBUG = True
                Ndim = W.shape[0]
                Nwaypoints = W.shape[1]
                Kcoeff = 4
                Ninterval = Nwaypoints-1
                P = np.zeros((Ninterval, Ndim, Kcoeff))
                durationVector = np.zeros((Ninterval))

                q0 = W[:,0]
                qd0 = dW[:,0]
                for j in range(0,Ninterval):
                        #q0 = W[:,j]
                        qd0 = dW[:,j]

                        qd0n = qd0/np.linalg.norm(qd0)
                        #qd0n = qd0

                        ### time info
                        if dt is None:
                                ds = np.dot(W[:,j+1]-q0,qd0n)
                        else:
                                ds = dt*np.linalg.norm(qd0)

                        qdd0 = W[:,j+1] - (q0 + ds*qd0n)
                        qdd0 = qdd0/(ds*ds)

                        a = q0
                        b = qd0n
                        c = qdd0
                        d = 0

                        P[j,:,0] = a
                        P[j,:,1] = b
                        P[j,:,2] = c
                        P[j,:,3] = d

                        P[j,2,1]=0
                        P[j,2,2]=0
                        P[j,2,3]=0

                        durationVector[j] = ds

                        [q0,qd0] = self.EvalPoly(P,j,ds)

                        #q1 = a + ds*qd0 + 0.5*ds*ds*qdd0
                        if np.linalg.norm(q0-W[:,j+1])>1e-10:
                                #print "d(q0,q1)",np.linalg.norm(W[:,j]-W[:,j+1])
                                #print "q0",W[:,j]
                                #print "q1",W[:,j+1]
                                print "mistmatch at",j,j+1,"dist",np.linalg.norm(q0-W[:,j+1])
                                sys.exit(0)
                        #qd0 = dqnext

                #sys.exit(0)

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
                        x += P[i,:,k]*np.power(t,k)
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
                        tstep = D[i]/100

                        WTplt.append(Dall)
                        Wplt.append(P[i,:,1])
                        Dall+=D[i]
                        while t <= D[i]:
                                tvec.append(tall)

                                [x,dx] = self.EvalPoly(P,i,t)
                                X.append(x)
                                dX.append(dx)
                                t+=tstep
                                tall += tstep

                        [q1spl,qd1spl] = self.EvalPoly(P,i,D[i])
                        q1 = W[:,i+1]

                        if np.linalg.norm(q1-q1spl)>1e-10:
                                print "qnext mismatch at",i,"/",M
                                print "q1",q1
                                print "q1spl",q1spl
                                sys.exit(0)

                        q0 = W[:,i]
                        [q0spl,qd0spl] = self.EvalPoly(P,i,0)

                        if np.linalg.norm(q0-q0spl)>1e-10:
                                print "q0 mismatch"
                                print "q0",q0
                                print "q0spl",q0spl
                                sys.exit(0)

                tvec=np.array(tvec)
                X=np.array(X).T
                dX=np.array(dX).T
                WTplt=np.array(WTplt).T
                Wplt=np.array(Wplt).T
                #print WTplt.shape
                #print Wplt.shape
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
