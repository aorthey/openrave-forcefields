import abc
import sys
import time
import numpy as np
from trajectory import *
import copy
import parameters_dynamical_system as params

DEFORM_SUCCESS = 0
DEFORM_OK = 1
DEFORM_COLLISION = 2
DEFORM_NOPROGRESS = 3
DEFORM_MAX_ITER_REACHED = 4

class Deformation():
        DEBUG = 0
        __metaclass__ = abc.ABCMeta

        ## number of steps of a linear homotopy deformation between two given
        ## trajectories
        DEFORMATION_STEPS = 50
        COLLISION_ENABLED = True
        PAUSE_BEFORE_DEFORM = False

        env = []
        traj_ori = []
        traj_deformed = []
        traj_current = []

        traj_display = []
        handle = []
        forcehandle = []

        critical_pt = None

        def __init__(self, trajectory, environment):
                self.env = environment

                ## traj_ori : original trajectory before any calls
                ## traj_current : the current deformation
                ## traj_deformed : the working copy of current for the next iteration
                ## traj_display : the trajectory which we display at the moment
                self.traj_ori = trajectory 
                self.traj_current = copy.copy(trajectory)
                self.traj_deformed = copy.copy(trajectory)
                self.traj_display = copy.copy(trajectory)
                self.critical_pt = 0

        def deform(self, N_iter = 1):
                computeNewCriticalPoint = True
                for i in range(0,N_iter):
                        res = self.deform_onestep(computeNewCriticalPoint)
                        self.traj_current = self.traj_deformed
                        self.draw_deformation()
                        if res == DEFORM_SUCCESS:
                                print "DEFORM_SUCCESS"
                                return True
                        elif res == DEFORM_OK:
                                print "DEFORM_OK",
                        elif res == DEFORM_COLLISION:
                                print "DEFORM_COLLISION",
                                return False
                        elif res== DEFORM_NOPROGRESS:
                                print "DEFORM_NOPROGRESS",
                                return False
                        else:
                                print "DEFORM_UNKNOWN"
                                print "(Deformation returned unknown error code)"
                                sys.exit(0)

                        print "(Iteration:",i,")"

                print "DEFORM_MAX_ITER_REACHED"
                return False

        def plot_topp(self):
                self.traj_deformed.PlotParametrization(self.env)

        def execute(self, robot, tsleep=0.01, stepping=False):
                tstep = 0.01
                if self.traj_deformed.topp.traj0 is not None:
                #if deform_success:
                        #td.traj_deformed.PlotParametrization(env)

                        xt = self.traj_deformed.topp.traj0

                        with self.env.env:
                                robot.GetLinks()[0].SetStatic(True)
                                self.env.env.StopSimulation() 

                        t = 0.0
                        tstep = 0.01
                        robot.SetDOFValues(xt.Eval(t))
                        self.env.MakeRobotVisible()

                        while t < xt.duration:
                                q = xt.Eval(t)
                                dq = xt.Evald(t)
                                ddq = xt.Evaldd(t)

                                qn = q + tstep*dq + 0.5*tstep*tstep*ddq
                                robot.SetDOFValues(qn)

                                self.env.env.StepSimulation(tstep)
                                time.sleep(tsleep)
                                t += tstep
                                if stepping:
                                        raw_input('Press Key to Step. Time: '+str(t)+'/'+str(xt.duration))

                        robot.WaitForController(0)
                        print "xt = td.traj_current.topp.traj0"
                        raw_input('Enter any key to quit. ')
                else:
                        print "deformation not successful"
                        raw_input('Enter any key to quit. ')

        def GetForcesAtWaypoints(self, W):
                Ndim = W.shape[0]
                if W.ndim > 1:
                        Nwaypoints = W.shape[1]
                else:
                        Nwaypoints = 1
                        pt = np.array(((W[0],W[1],-0.1,0.001)))
                        F = np.zeros((Ndim))
                        F[0:3] = self.env.GetForceAtX(pt)
                        return F

                F = np.zeros((Ndim,Nwaypoints))
                for i in range(0,Nwaypoints):
                        pt = np.array(((W[0,i],W[1,i],-0.1,0.001)))
                        F[0:3,i] = self.env.GetForceAtX(pt)
                return F

        def draw_forces_at_waypoints(self, W, F):
                for i in range(0,W.shape[1]):
                        df = np.linalg.norm(F[:,i])
                        if df > 0.3:
                                nF = 0.3*F[:,i]/df
                        else:
                                nF = F[:,i]
                        W1 = np.array((W[0,i],W[1,i],W[2,i]))
                        W2 = np.array((W[0,i]+nF[0],W[1,i]+nF[1],W[2,i]+nF[2]))
                        P = np.array(((W1),(W2)))
                        ls = self.traj_current.linsize
                        h=self.env.env.drawlinestrip(points=P,linewidth=ls,colors=np.array(((0.8,0.2,0.2,0.9))))
                        self.forcehandle.append(h)

        def InterpolateWaypoints(self, W0, W1, k):
                ### k: between [0,1]
                Ndim = W0.shape[0]
                Nwaypoints = W0.shape[1]
                Wk = np.zeros((Ndim,Nwaypoints))
                Wk[0:3,:] = (1-k)*W0[0:3,:] + k*W1[0:3,:]

                for i in range(0,Nwaypoints):
                        t0 = W0[3,i]
                        t1 = W1[3,i]

                        d1 = abs(t1-t0)
                        d2 = 2*pi-abs(t1-t0)

                        if d1 > d2:
                                ### interpolation needs to cross boundary
                                ### between charts!
                                if t1 > t0:
                                        ## going into minus direction
                                        if k*d2 < abs(t0+pi):
                                                ## before chart boundary
                                                tnew = t0 - k*d2
                                        else:
                                                ## after chart boundary
                                                tnew = t1 + (d2-k*d2)
                                else:
                                        ### t0 > t1
                                        ## going into plus direction
                                        if k*d2 < abs(t0-pi):
                                                ## before chart boundary
                                                tnew = t0 + k*d2
                                        else:
                                                ## after chart boundary
                                                tnew = t1 - (d2-k*d2)

                                Wk[3,i] = tnew
                                #print "k:","t0:",t0,"t1:",t1,"d2:",d2,"tnew:",tnew
                        else:
                        #        print ">",d1,t1,t0
                                Wk[3,i] = (1-k)*W0[3,i] + k*W1[3,i]
                #sys.exit(0)
                return Wk



        def draw_deformation(self):
                ## visualize linear homotopy between two trajectories

                [W0,dW] = self.traj_display.get_waypoints() 
                Nwaypoints = W0.shape[1]
                [W1,dW] = self.traj_current.get_waypoints(Nwaypoints) 
                self.traj_current.new_from_waypoints(W1)

                self.forcehandle = []
                if np.linalg.norm(W0-W1)<1e-10:
                        print "No deformation"
                        self.handle = self.traj_current.draw(self.env, keep_handle=False)
                        return

                t1 = time.time()

                tdraw = 0.0
                for i in range(0,self.DEFORMATION_STEPS):
                        k = float(i)/float(self.DEFORMATION_STEPS)

                        Wk = (1-k)*W0 + k*W1
                        self.traj_current.new_from_waypoints(Wk)

                        ti1 = time.time()
                        self.handle = self.traj_current.draw(self.env, keep_handle=False, critical_pt = self.critical_pt)
                        ti2 = time.time()

                        if self.PAUSE_BEFORE_DEFORM and i==0:
                                raw_input('Press <ENTER> to start deforming.')
                        tdraw += ti2-ti1

                self.handle = self.traj_current.draw(self.env, keep_handle=False) 
                t2 = time.time()
                if self.DEBUG:
                        print "DEFORMATION TIME(s): ",t2-t1," (drawing:",tdraw,")"

                self.traj_display = copy.copy(self.traj_current)


        def draw_trajectory_original(self):
                self.traj_ori.draw(self.env)

        def draw_trajectory_deformed(self):
                self.traj_current.draw(self.env)

        def extractInfoFromTrajectory(self, traj):
                eta = 1.0

                L = traj.get_length()

                ### traj.DISCRETIZATION_TIME_STEP controls ds
                [Wori,dWori,ddWori] = traj.get_waypoints_second_order()
                [Ndim, Nwaypoints] = traj.getWaypointDim(Wori)
                F = traj.get_forces_at_waypoints(Wori, self.env)
                [R,amin,amax] = traj.getControlMatrix(Wori)
                ### get forces in normal direction to trajectory
                FN = self.getForceNormalComponent(F, dWori)

                self.critical_pt = traj.getCriticalPointFromWaypoints(self.env, Wori, dWori, ddWori, 0)

                #### compute min/max velocity profile from path without forces
                #### (if available). otherwise use [0,0]
                dpmin = np.zeros((1,Nwaypoints))
                dpmax = np.zeros((1,Nwaypoints))

                self.traj_velprofile = traj.getVelocityIntervalWithoutForceField(self.env, Wori, dWori, ddWori)

                if self.traj_velprofile is not None:
                        Tend = self.traj_velprofile.duration
                        Tstart = 0.0
                        Tstep = Tend/1e4

                        for i in range(0,Nwaypoints):
                                Tcur =Tstart
                                p = Wori[:,i]
                                q = self.traj_velprofile.Eval(Tcur)

                                dold = 1e5
                                dnew = np.linalg.norm(p-q)
                                while dnew < dold:
                                        dold = dnew
                                        Tcur += Tstep
                                        q = self.traj_velprofile.Eval(Tcur)
                                        dnew = np.linalg.norm(p-q)

                                dq = self.traj_velprofile.Evald(Tcur)
                                dpmax[:,i] = np.linalg.norm(dq)
                                Tstart = Tcur


                DeformInfo = {}
                DeformInfo['Ndim'] = Ndim
                DeformInfo['Nwaypoints'] = Nwaypoints
                DeformInfo['smoothing_factor'] = self.smoothing_factor
                DeformInfo['critical_pt'] = self.critical_pt
                DeformInfo['traj'] = traj
                DeformInfo['Wori'] = Wori
                DeformInfo['dWori'] = dWori
                DeformInfo['ddWori'] = ddWori
                DeformInfo['F'] = F
                DeformInfo['FN'] = FN
                DeformInfo['R'] = R
                DeformInfo['amin'] = amin
                DeformInfo['amax'] = amax
                DeformInfo['dpmin'] = dpmin
                DeformInfo['dpmax'] = dpmax
                DeformInfo['eta'] = eta
                DeformInfo['env'] = self.env
                return DeformInfo

        def IsProjectable(self, DeformInfo):

                env = DeformInfo['env']
                traj = DeformInfo['traj']
                Wori = DeformInfo['Wori']
                dWori = DeformInfo['dWori']
                ddWori = DeformInfo['ddWori']
                Ndim = DeformInfo['Ndim']
                Nwaypoints = DeformInfo['Nwaypoints']
                critical_pt = DeformInfo['critical_pt']

                FT = traj.get_tangent_forces_at_waypoints(Wori, dWori, env)
                topp = traj.getTOPPTrajectoryWithTangentForceField(env, Wori, dWori)
                #topp = traj.getTOPPTrajectoryWithoutForceField(env, Wori, dWori)

                #topp.PlotTrajectory(env)
                ### execute blindly!?
                xt = self.traj_deformed.topp.traj1

                Q0 = []
                QD0 = []
                Qori = []
                Qdori = []
                D = 0

                t = 0.0
                dt = 0.1

                q = xt.Eval(0)
                dq = xt.Evald(0)
                while t < xt.duration:
                        #dq = xt.Evald(t)
                        #print np.linalg.norm(xt.Evald(t))
                        #print np.linalg.norm(dq),np.linalg.norm(dq[0]),np.linalg.norm(dq[1])
                        #print np.linalg.norm(dq[0])
                        ### get component of force orthogonal to path

                        F = self.GetForcesAtWaypoints(q)
                        dqn = xt.Evald(t)/np.linalg.norm(xt.Evald(t))
                        FN = F - np.dot(F,dqn)*dqn

                        ddq = xt.Evaldd(t)
                        #GetNearestControlPoint(q, dq, pnext, F, dt, speed, Acontrol, bcontrol, returnControl=False):
                        [qnew, ddq_adjust] = params.GetBestControlPathInvariant( q, dq, ddq, xt.Eval(t+dt), F, dt)

                        Qori.append(xt.Eval(t))
                        Qdori.append(xt.Evald(t))
                        Q0.append(q)
                        QD0.append(dq)

                        q0 = copy.copy(q)

                        dt2 = 0.5*dt*dt

                        #q = q + dt*dq + dt2*ddq + dt2*FN
                        q = q + dt*dq + dt2*ddq_adjust

                        ### compute here the best control to reject normal
                        ### component of force
                        #dq = dq + dt*ddq + dt*FN
                        dq = dq + dt*ddq_adjust

                        D+=np.linalg.norm(q-q0)
                        t += dt
                        print "t",t,"/",xt.duration,"speed",np.linalg.norm(dq)

                print "overall dist",D
                ### do not move waypoint derivatives
                Q0 = np.array(Q0).T
                QD0 = np.array(QD0).T
                Qori = np.array(Qori).T
                Qdori = np.array(Qdori).T

                plt.plot(Qori[0,:],Qori[1,:],'-ok',linewidth=3,markersize=5)
                #plt.plot(Wori[0,:],Wori[1,:],'og',markersize=12)
                plt.plot(Q0[0,:],Q0[1,:],'or',linewidth=3)

                #for i in range(0,QD0.shape[1]):
                        #tangent = 0.1*QD0[0:2,i]/np.linalg.norm(QD0[0:2,i])
                        #plt.plot([Q0[0,i],Q0[0,i]+tangent[0]],[Q0[1,i],Q0[1,i]+tangent[1]],'-m')
                #for i in range(0,QD0.shape[1]):
                        #tangent = 0.05*Qdori[0:2,i]/np.linalg.norm(Qdori[0:2,i])
                        #plt.plot([Qori[0,i],Qori[0,i]+tangent[0]],[Qori[1,i],Qori[1,i]+tangent[1]],'-m',linewidth=2)

                #for i in range(0,Wori.shape[1]):
                        #tangent = 0.05*dWori[0:2,i]/np.linalg.norm(dWori[0:2,i])
                        #plt.plot([Wori[0,i],Wori[0,i]+tangent[0]],[Wori[1,i],Wori[1,i]+tangent[1]],'-g',linewidth=3)

                ### create new trajectory string for that
                plt.show()
                #sys.exit(0)


                print "waypoints:",Q0.shape[1]
                N = Q0.shape[1]
                Q0 = Q0[:,0:N]
                QD0 = QD0[:,0:N]
                #Qori = Qori[:,0:N]
                t2 = Trajectory(Q0)
                [smin,smax] = t2.GetSpeedIntervalAtCriticalPoint(env, Q0, QD0, N-10)
                print "speed interval",smin,smax

                #t2.PlotParametrization(env)
                #topp = t2.getTOPPTrajectoryWithoutForceField(env, Q0, QD0)
                #topp.PlotTrajectory(env)

                #t2.PlotParametrization2(env,Q0,QD0)

                sys.exit(0)

        def IsDynamicallyFeasible(self, DeformInfo):

                traj = DeformInfo['traj']
                Wori = DeformInfo['Wori']
                dWori = DeformInfo['dWori']
                ddWori = DeformInfo['ddWori']
                Ndim = DeformInfo['Ndim']
                Nwaypoints = DeformInfo['Nwaypoints']
                critical_pt = DeformInfo['critical_pt']
                self.critical_pt = traj.getCriticalPointFromWaypoints(self.env, Wori, dWori, ddWori, critical_pt)

                print "###########################################"
                print "CRITICAL WAYPOINT: ",self.critical_pt,"/",Nwaypoints

                if self.critical_pt >= Nwaypoints-1:
                        print "No deformation necessary => Trajectory dynamically feasible"
                        print traj.getCriticalPointFromWaypoints(self.env, Wori, dWori, ddWori, self.critical_pt)
                        print "###########################################"
                        return True
                else:
                        self.Nc_handle = self.env.env.plot3(points=Wori[0:3,critical_pt],
                                        pointsize=0.02,
                                        colors=np.array(((1.0,0.2,0.2))),
                                        drawstyle=1)
                print "###########################################"
                return False

        def getForceNormalComponent(self, F, dWori):
                Nwaypoints = dWori.shape[1]
                FN = np.zeros((F.shape))
                for i in range(0,Nwaypoints):
                        dWn = dWori[:,i]/np.linalg.norm(dWori[:,i])
                        FN[:,i] = F[:,i]-np.dot(F[:,i],dWn)*dWn
                return FN

        def SafetyCheckUpdate(self, dU):
                epsilon=1e-10
                dstart = np.linalg.norm(dU[:,0])
                dend = np.linalg.norm(dU[:,-1])
                if dstart > epsilon:
                        print "ERROR: start point deformation -> not allowed"
                        print dU[:,0],dstart
                        sys.exit(1)
                if dend > epsilon:
                        print "ERROR: end point deformation -> not allowed"
                        print dU[:,-1],dend
                        sys.exit(1)
