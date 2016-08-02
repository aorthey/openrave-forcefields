import sys
import time
import numpy as np
from trajectory import *
from trajectory import Trajectory
import copy
import parameters_dynamical_system as params

class ProjectorSimple():

        def __init__(self):
                pass

        def IsProjectable(self, DeformInfo):
                #return self.IsProjectableTOPP(DeformInfo)
                return self.IsProjectableReachableSet(DeformInfo)

        projectedWaypoints = None

        def GetProjectableWaypoints(self):
                if self.projectedWaypoints is None:
                        print "not projectable, cannot return waypoints"
                        sys.exit(0)
                return self.projectedWaypoints

        def GetTimeAtCriticalPoint( self, xt, Wori, CP):
                Wc = Wori[:,CP]

                t =0.0
                dt = 0.01
                dbest = 1e5
                tc = -1
                while t < xt.duration:
                        W = xt.Eval(t)
                        d = np.linalg.norm(W-Wc)
                        if d<dbest:
                                dbest = d
                                tc = t
                        t +=dt
                print "CP t=",tc,"/",xt.duration

                return tc


        def IsProjectableReachableSet(self, DeformInfo):
                env = DeformInfo['env']
                traj = DeformInfo['traj']
                Wori = DeformInfo['Wori']
                dWori = DeformInfo['dWori']
                ddWori = DeformInfo['ddWori']
                Ndim = DeformInfo['Ndim']
                Nwaypoints = DeformInfo['Nwaypoints']
                critical_pt = DeformInfo['critical_pt']
                topp = traj.getTOPPTrajectoryWithoutForceField(env, Wori)

                ### strategy: execute blindly but reject forces
                xt = topp.traj1
                tc = self.GetTimeAtCriticalPoint( xt, Wori, critical_pt)

                dt = Trajectory.DISCRETIZATION_TIME_STEP
                Q0 = []
                QD0 = []
                Qori = []
                Qdori = []
                D = 0
                t = 0

                ictr=0

                while t < tc:
                        Qori.append(xt.Eval(t))
                        Qdori.append(xt.Evald(t))
                        Q0.append(xt.Eval(t))
                        QD0.append(xt.Evald(t))
                        t+=dt

                q = xt.Eval(tc)
                dq = xt.Evald(tc)
                while t < xt.duration:
                        Qori.append(xt.Eval(t))
                        Qdori.append(xt.Evald(t))
                        Q0.append(q)
                        QD0.append(dq)
                        #dt = topp.durationVector[ictr]
                        ictr+=1

                        qprev = copy.copy(q)
                        dqprev = copy.copy(dq)
                        F = traj.get_forces_at_waypoints(q, env)

                        ### get component of force orthogonal to path
                        #dqn = xt.Evald(t)/np.linalg.norm(xt.Evald(t))
                        #FN = F - np.dot(F,dqn)*dqn

                        ### the control we would apply if we execute the
                        ### geometrical feasible path without forces
                        ddq = xt.Evaldd(t)

                        ### Get a control which follows ddq as close as possible
                        ### while rejecting the forces as much as possible.
                        [qnew, ddq] = params.GetBestControlPathInvariant( q, dq, ddq, xt.Eval(t+dt), xt.Evald(t+dt), F, dt)

                        ### simulate forward with new control
                        dt2 = 0.5*dt*dt
                        q = q + dt*dq + dt2*ddq
                        dq = dq + dt*ddq

                        ### check that forward simulation is really reaching the
                        ### next point
                        if np.linalg.norm(q-qnew)>1e-10:
                                print "q mismatch dist",np.linalg.norm(q-qnew)
                                sys.exit(0)

                        t += dt
                        D+=np.linalg.norm(q-qprev)

                print "overall dist",D
                ### do not move waypoint derivatives
                Q0 = np.array(Q0).T
                Qori = np.array(Qori).T
                self.projectedWaypoints = Q0

                #if np.linalg.norm(q-xt.Eval(xt.duration))<0.1:
                #        print "success"
                #        return True
                #else:
                #        return False

                ### check if topp finds solution
                t2 = Trajectory(Q0)
                Nc = t2.getCriticalPoint(env)

                t3 = Trajectory(Qori)
                Nc3 = t3.getCriticalPoint(env)


                traj2 = t2.topp.traj0
                traj3 = t3.topp.traj0

                Npts = 1e5
                tvect2 = np.linspace(0,traj2.duration, Npts)
                tvect3 = np.linspace(0,traj3.duration, Npts)
                qvect2 = np.array([traj2.Eval(t) for t in tvect2]).T
                qvect3 = np.array([traj3.Eval(t) for t in tvect3]).T

                plt.plot(qvect2[0,:],qvect2[1,:],'-r',linewidth=3)
                plt.plot(qvect3[0,:],qvect3[1,:],'-g',linewidth=3)
                plt.plot(Q0[0,Nc],Q0[1,Nc],'or',markersize=5)
                plt.plot(Qori[0,Nc3],Qori[1,Nc3],'og',markersize=5)

                print "BEFORE SIMPLE PROJECTABILITY:",Nc3,"/",Qori.shape[1]
                print "AFTER  SIMPLE PROJECTABILITY:",Nc,"/",Q0.shape[1]

                tn = Trajectory(Q0[:,0:Nc])
                tn.getCriticalPoint(env)

                for k in range(1,10):
                        Nk = Nc - k
                        print Nk,tn.topp.getSpeedIntervalAtPoint(Nk)
                plt.show()

                if Nc >= Q0.shape[1]-1:
                        return True
                else:
                        return False

