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

                dt = 0.01
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
                QD0 = np.array(QD0).T
                Qori = np.array(Qori).T
                Qdori = np.array(Qdori).T

                self.projectedWaypoints = Q0

                #plt.plot(Qori[0,:],Qori[1,:],'-ok',linewidth=3,markersize=5)
                #plt.plot(Q0[0,:],Q0[1,:],'-or',linewidth=3)
                #plt.show()

                #if np.linalg.norm(q-xt.Eval(xt.duration))<0.1:
                #        print "success"
                #        return True
                #else:
                #        return False

                ### check if topp finds solution
                print "waypoints:",Q0.shape[1]

                t2 = Trajectory(Q0)
                t2.info()
                Nc = t2.getCriticalPoint(env)
                if Nc >= Q0.shape[1]-1:
                        return True
                else:
                        return False

        def IsProjectableReachableSetDraw(self, DeformInfo):
                env = DeformInfo['env']
                traj = DeformInfo['traj']
                Wori = DeformInfo['Wori']
                dWori = DeformInfo['dWori']
                ddWori = DeformInfo['ddWori']
                Ndim = DeformInfo['Ndim']
                Nwaypoints = DeformInfo['Nwaypoints']
                critical_pt = DeformInfo['critical_pt']
                #FT = traj.get_tangent_forces_at_waypoints(Wori, dWori, env)
                #topp = traj.getTOPPTrajectoryWithTangentForceField(env, Wori, dWori)
                topp = traj.getTOPPTrajectoryWithoutForceField(env, Wori, dWori)

                #topp.PlotTrajectory(env)
                ### execute blindly!?
                xt = topp.traj1

                print xt.duration
                #print topp.durationVector
                #print topp.durationVector
                #print topp.durationVector.shape
                #print np.sum(topp.durationVector)
                #print topp.durationVector[0]

                dt = 0.01
                for eta in np.linspace(0.9,1.5,5):
                        Q0 = []
                        QD0 = []
                        Qori = []
                        Qdori = []
                        D = 0
                        t = 0.0

                        q = xt.Eval(0)
                        dq = xt.Evald(0)
                        ictr=0
                        while t < xt.duration:
                                Qori.append(xt.Eval(t))
                                Qdori.append(xt.Evald(t))
                                Q0.append(q)
                                QD0.append(dq)
                                #dt = topp.durationVector[ictr]
                                ictr+=1

                                qprev = copy.copy(q)
                                dqprev = copy.copy(dq)

                                F = eta*traj.get_forces_at_waypoints(q, env)

                                ### get component of force orthogonal to path
                                #dqn = xt.Evald(t)/np.linalg.norm(xt.Evald(t))
                                #FN = F - np.dot(F,dqn)*dqn

                                ### the control we would apply if we execute the
                                ### geometrical feasible path without forces
                                ddq = xt.Evaldd(t)

                                ### Get a control which follows ddq as close as possible
                                ### while rejecting the forces as much as possible.
                                [qnew, ddq] = params.GetBestControlPathInvariant( q, dq, ddq, xt.Eval(xt.duration), xt.Evald(t+dt), F, dt)

                                ### simulate forward with new control
                                dt2 = 0.5*dt*dt
                                q = q + dt*dq + dt2*ddq
                                dq = dq + dt*ddq

                                if np.linalg.norm(q-qnew)>1e-10:
                                        print "q mismatch dist",np.linalg.norm(q-qnew)
                                        sys.exit(0)
                                t += dt
                                D+=np.linalg.norm(q-qprev)

                        print "overall dist",D
                        ### do not move waypoint derivatives
                        Q0 = np.array(Q0).T
                        QD0 = np.array(QD0).T
                        Qori = np.array(Qori).T
                        Qdori = np.array(Qdori).T

                        plt.plot(Qori[0,:],Qori[1,:],'-ok',linewidth=3,markersize=5)
                        #plt.plot(Wori[0,:],Wori[1,:],'og',markersize=12)
                        plt.plot(Q0[0,:],Q0[1,:],'-or',linewidth=3)
                plt.show()


                print "waypoints:",Q0.shape[1]
                N = Q0.shape[1]
                t2 = Trajectory(Q0)
                Nc = t2.getCriticalPoint(env)
                if Nc >= Nwaypoints-1:
                        return True
                else:
                        return False

        def IsProjectableTOPP(self, DeformInfo):

                env = DeformInfo['env']
                traj = DeformInfo['traj']
                Wori = DeformInfo['Wori']
                dWori = DeformInfo['dWori']
                ddWori = DeformInfo['ddWori']
                Ndim = DeformInfo['Ndim']
                Nwaypoints = DeformInfo['Nwaypoints']
                critical_pt = DeformInfo['critical_pt']

                #FT = traj.get_tangent_forces_at_waypoints(Wori, dWori, env)
                #topp = traj.getTOPPTrajectoryWithTangentForceField(env, Wori, dWori)
                topp = traj.getTOPPTrajectoryWithoutForceField(env, Wori, dWori)
                #topp.PlotTrajectory(env)
                ### execute blindly!?
                xt = topp.traj1

                print xt.duration
                #print topp.durationVector
                #print topp.durationVector
                #print topp.durationVector.shape
                #print np.sum(topp.durationVector)
                #print topp.durationVector[0]

                Q0 = []
                QD0 = []
                Qori = []
                Qdori = []
                D = 0
                t = 0.0

                q = xt.Eval(0)
                dq = xt.Evald(0)
                ictr=0
                while t < xt.duration:
                        Qori.append(xt.Eval(t))
                        Qdori.append(xt.Evald(t))
                        Q0.append(q)
                        QD0.append(dq)
                        dt = topp.durationVector[ictr]
                        ictr+=1

                        if ictr>2:
                                W0 = np.array(Q0).T
                                dW0 = np.array(QD0).T
                                traj0 = Trajectory(W0)
                                [smin2,smax2] = \
                                traj0.GetSpeedIntervalAtCriticalPoint(env, W0, dW0, ictr-1, dt)
                                [smin,smax] = \
                                traj0.GetSpeedIntervalAtCriticalPoint(env, W0, dW0, ictr, dt)

                                print "t",t,"/",xt.duration,"speed",np.linalg.norm(dq),"tspeed",smin,smax,"q",q,"F",F

                                if smax < 0.001:
                                        print ictr
                                        print "W",W0[:,-1],W0[:,-2]
                                        plot([q0[0],q[0]],[q0[1],q[1]],'-ok',markersize=10)
                                        plot([xt.Eval(t)[0],q0[0]],[xt.Eval(t)[1],q0[1]],'-om',markersize=10)
                                        #plot([traj0.Eval(t)[0],q0[0]],[traj0.Eval(t)[1],q0[1]],'-om',markersize=10)
                                        #[smin,smax] = traj0.GetSpeedIntervalAtCriticalPoint(env, W0, dW0, ictr-1, dt)
                                        print smax

                                        plot([q0[0],q0[0]+dt * dq0[0]],[q0[1],q0[1]+dt * dq0[1]],'-g')
                                        plot([q0[0],q0[0]+smax2*dt * dq0[0]],[q0[1],q0[1]+smax2*dt * dq0[1]],'-og',linewidth=4)
                                        plot([q0[0],q0[0]+dt*dt*0.5*F[0]],[q0[1],q0[1]+dt*dt*0.5*F[1]],'-or',linewidth=4)

                                        #k=0
                                        #while t > k*topp.durationQ:
                                        #        k+=1
                                        #print topp.a[k,:]
                                        #print topp.b[k,:]
                                        #print topp.c[k,:]

                                        t=0
                                        M = ictr-2

                                        #xt = traj0.topp
                                        #dv = xt.durationVector
                                        tstep = dt/15.0

                                        QN = []
                                        QT = []
                                        while t<=traj0.topp.durationVector[M]:
                                                print t,"/",traj0.topp.durationVector[M]
                                                [qn,dqn] = traj0.topp.EvalPoly(traj0.topp.poly,M,t)
                                                QN.append(qn)
                                                t+=tstep
                                                Msamples = 100
                                                k=0
                                                while k<Msamples:
                                                        u = params.GetRandomControl()
                                                        G = params.GetControlMatrixAtWaypoint(q0)
                                                        ddq = np.dot(G,u) + F
                                                        t2 = 0.5*t*t
                                                        qn = q0 + t*dq0 + t2*ddq
                                                        QT.append(qn)
                                                        k+=1

                                        Msamples = 4
                                        k=-Msamples
                                        QK = []
                                        QR = []

                                        Wstart = copy.copy(Q0[-1])

                                        dt = traj0.topp.durationVector[M]
                                        while k <= Msamples:
                                                m=-Msamples
                                                while m <= Msamples:
                                                        ddqn = copy.copy(ddq)
                                                        #ddqn[0] = np.random.normal(F[0], 0.0, 1)
                                                        #ddqn[1] = np.random.normal(F[1], 0.0, 1)
                                                        #ddqn[3] = np.random.normal(F[3], 0.0, 1)
                                                        ddqn = k*2*F + m*dq0 
                                                        dt2 = 0.5*dt*dt
                                                        q = q0 + smax2*dt*dq0 + dt2*ddqn
                                                        dq = smax2*dq0 + dt*ddqn
                                                        #Wnext = copy.copy(Wstart)
                                                        #Wnext[0] = np.random.normal(Wstart[0], sigma, 1)
                                                        #Wnext[1] = np.random.normal(Wstart[1], sigma, 1)
                                                        #Wnext[3] = np.random.normal(Wstart[3], sigma, 1)
                                                        del Q0[-1]
                                                        del QD0[-1]

                                                        Q0.append(q)
                                                        QD0.append(dq)
                                                        W0 = np.array(Q0).T
                                                        dW0 = np.array(QD0).T
                                                        traj0 = Trajectory(W0)
                                                        [smin,smax] = \
                                                        traj0.GetSpeedIntervalAtCriticalPoint(env, W0, dW0, ictr, dt)
                                                        print "SAMPLE",k,"/",Msamples,smax
                                                        if smax > 0.001:
                                                                QK.append(q)
                                                        else:
                                                                QR.append(q)
                                                        m+=1
                                                k+=1

                                        #plot([q0[0],x[0]],[q0[1],x[1]],'-om',markersize=10)

                                        if QK:
                                                QK = np.array(QK).T
                                                plt.plot(QK[0,:],QK[1,:],'og',markersize=10)
                                        if QR:
                                                QR = np.array(QR).T
                                                plt.plot(QR[0,:],QR[1,:],'or',markersize=10)

                                        QT = np.array(QT).T
                                        plt.plot(QT[0,:],QT[1,:],'or')

                                        QN = np.array(QN).T
                                        plt.plot(QN[0,:],QN[1,:],'-om',markersize=10)
                                        plt.show()
                                        sys.exit(0)


                        q0 = copy.copy(q)
                        dq0 = copy.copy(dq)

                        F = traj.get_forces_at_waypoints(q0, env)
                        ### get component of force orthogonal to path
                        #F = self.GetForcesAtWaypoints(q0)
                        #dqn = xt.Evald(t)/np.linalg.norm(xt.Evald(t))
                        #FN = F - np.dot(F,dqn)*dqn

                        ddq = xt.Evaldd(t)
                        #k = 0
                        #while t > k*topp.durationQ:
                        #        k+=1

                        #print topp.a[0:k,:]
                        #print topp.b[0:k,:]
                        #print topp.c[0:k,:]
                        [qnew, ddq] = params.GetBestControlPathInvariant( q, dq, ddq, xt.Eval(t+dt), F, dt)

                        #ds = np.dot(W[:,j+1]-q0,qd0n)
                        dt2 = 0.5*dt*dt
                        q = q + dt*dq + dt2*ddq

                        #q = q + dt*dq + dt2*ddq_adjust
                        ### compute here the best control to reject normal
                        ### component of force
                        dq = dq + dt*ddq
                        #dq = dq + dt*ddq_adjust

                        if np.linalg.norm(q-qnew)>1e-10:
                                print "q mismatch dist",np.linalg.norm(q-qnew)
                                sys.exit(0)
                        t += dt
                        D+=np.linalg.norm(q-q0)

                print "overall dist",D
                ### do not move waypoint derivatives
                Q0 = np.array(Q0).T
                QD0 = np.array(QD0).T
                Qori = np.array(Qori).T
                Qdori = np.array(Qdori).T

                plt.plot(Qori[0,:],Qori[1,:],'-ok',linewidth=3,markersize=5)
                #plt.plot(Wori[0,:],Wori[1,:],'og',markersize=12)
                plt.plot(Q0[0,:],Q0[1,:],'-or',linewidth=3)


                print "waypoints:",Q0.shape[1]
                N = Q0.shape[1]
                t2 = Trajectory(Q0)
                Nc = t2.getCriticalPoint(env)
                if Nc >= Nwaypoints-1:
                        return True
                else:
                        return False
                #plt.show()

                #t2.PlotParametrization(env)
                #topp = t2.getTOPPTrajectoryWithoutForceField(env, Q0, QD0)
                #topp.PlotTrajectory(env)

                #t2.PlotParametrization2(env,Q0,QD0)

                #sys.exit(0)
