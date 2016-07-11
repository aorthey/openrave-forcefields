import sys
import time
import numpy as np
from trajectory import *
import copy
import parameters_dynamical_system as params

class ProjectorSimple():

        def __init__(self):
                pass

        def IsProjectable(self, DeformInfo):

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
                #sys.exit(0)

                Q0 = []
                QD0 = []
                Qori = []
                Qdori = []
                D = 0
                t = 0.0
                dt = 0.01
                #dt = topp.durationQ
                #print dt
                #sys.exit(0)

                q = xt.Eval(0)
                dq = xt.Evald(0)
                N = int(xt.duration/dt)
                ictr=0
                F = 0
                while t < xt.duration:
                        Qori.append(xt.Eval(t))
                        Qdori.append(xt.Evald(t))
                        Q0.append(q)
                        QD0.append(dq)
                        ictr+=1

                        if ictr>2:
                                W0 = np.array(Q0).T
                                dW0 = np.array(QD0).T
                                traj0 = Trajectory(W0)
                                [smin,smax] = \
                                traj0.GetSpeedIntervalAtCriticalPoint(env, W0, dW0, ictr, dt)
                                print "t",t,"/",xt.duration,"speed",np.linalg.norm(dq),"tspeed",smin,smax,"F",F,traj0.topp.c[-1,:]
                                if smax < 0.001:
                                        dt2 = 0.5*dt*dt
                                        plot([q0[0],q[0]],[q0[1],q[1]],'-ok',markersize=10)
                                        plot([xt.Eval(t)[0],q0[0]],[xt.Eval(t)[1],q0[1]],'-om',markersize=10)

                                        print "PS",q


                                        t=0
                                        M = ictr-2

                                        xt = traj0.topp
                                        dv = xt.durationVector
                                        tstep = dv[M]/15.0
                                        print dv,dv[M]
                                        QN = []
                                        print W0.shape,traj0.topp.c.shape
                                        print traj0.topp.c
                                        QT = []
                                        while t<=traj0.topp.durationVector[M]:
                                                [qn,dqn] = traj0.topp.EvalPoly(traj0.topp.poly,M,t)
                                                QN.append(qn)
                                                t+=tstep
                                                Msamples = 1000
                                                k=0
                                                while k<Msamples:
                                                        u = params.GetRandomControl()
                                                        G = params.GetControlMatrixAtWaypoint(q0)
                                                        ddq = np.dot(G,u) + F
                                                        t2 = 0.5*t*t
                                                        qn = q0 + t*dq0 + t2*ddq
                                                        QT.append(qn)
                                                        k+=1


                                        #P =np.array( [[-2.6479999424817375, -0.9996635848610663, 0.0, 0.0], [-0.0006290807563600611, 0.02590624508881563, 1.7122297269190248e-17, 0.0], [0.1, 0.0, 0.0, 0.0], [-3.141042628523347, 0.001258399030477833, -1.0701435793243905e-18, 0.0]] )
                                        #[x,dx] = EvalPoly(P,dt)

                                        #plot([q0[0],x[0]],[q0[1],x[1]],'-om',markersize=10)

                                        QT = np.array(QT).T
                                        plt.plot(QT[0,:],QT[1,:],'or')
                                        QN = np.array(QN).T
                                        plt.plot(QN[0,:],QN[1,:],'-om',markersize=10)
                                        plt.show()
                                        sys.exit(0)


                        q0 = copy.copy(q)
                        dq0 = copy.copy(dq)

                        #if ictr>2 and traj0.topp.c[-1,1]>0:
                                #F = np.array([0,1.5,0,0])
                        #else:

                        F = traj.get_forces_at_waypoints(q0, env)
                        #print "PS",q0,F
                        ### get component of force orthogonal to path
                        #F = self.GetForcesAtWaypoints(q0)
                        #dqn = xt.Evald(t)/np.linalg.norm(xt.Evald(t))
                        #FN = F - np.dot(F,dqn)*dqn

                        ddq = xt.Evaldd(t)
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
