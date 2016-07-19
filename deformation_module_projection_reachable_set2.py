import abc
import sys
import copy
import parameters_dynamical_system as params
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
from deformation_module import *
from cvxopt import matrix, solvers
from util import *

class DeformationModuleProjectionReachableSet2(DeformationModule):

        def get_gradient(self, lambda_coeff):

                traj = self.DeformInfo['traj']
                env = self.DeformInfo['env']
                eta = self.DeformInfo['eta']
                Wori = self.DeformInfo['Wori']
                dWori = self.DeformInfo['dWori']
                Ndim = self.DeformInfo['Ndim']
                Nwaypoints = self.DeformInfo['Nwaypoints']
                critical_pt = self.DeformInfo['critical_pt']
                F = self.DeformInfo['F']
                dpmin = self.DeformInfo['dpmin']
                dpmax = self.DeformInfo['dpmax']

                Wmove = np.zeros((Ndim,Nwaypoints))
                dUtmp = np.zeros((Ndim,Nwaypoints))
                if critical_pt > Nwaypoints-1:
                        return dUtmp


                [cp_semin, cp_semax] = traj.GetSpeedIntervalAtCriticalPoint(env, Wori, dWori, critical_pt)
                print "CP",critical_pt,"vel:[",cp_semin,",",cp_semax,"]"

                ### get max velocity at waypoint
                #### compute the intersection of a forward simulation and the
                #### ball with radius ds

                i=critical_pt

                p = Wori[:,i]
                dp = dWori[:,i]
                smax = cp_semax
                ds = traj.DISCRETIZATION_TIME_STEP

                Wnext = copy.copy(Wori)
                dWnext = copy.copy(dWori)

                offset = np.zeros((Ndim))
                Mforward = 100
                print "CP+1:",Wori[:,i+1]

                Wk = copy.copy(Wori[0:4,critical_pt:critical_pt+Mforward+1])
                Qk = []
                Qk.append(Wori[0:4,critical_pt])
                while i < critical_pt+Mforward and i<Nwaypoints-1:
                        #p = Wori[:,i]
                        #dp = dWori[:,i]
                        #dpnext = dWori[:,i+1]

                        pnext = Wori[:,i+1]
                        Fcur = F[:,i]
                        ds = np.linalg.norm(p-pnext)
                        [qnext,dqnext,qdd,tt] = params.ForwardSimulate(p, dp, smax, ds, Fcur, pnext=Wori[:,i+1])
                        Wmove[:,i+1] = qnext - pnext

                        #DEBUG=1
                        #if DEBUG:
                        #        Wori[:,i+1] = copy.copy(qnext)
                        #        dWori[:,i+1] = qnext-p
                        #        dp = qnext-p
                        #        ddd = np.linalg.norm(dp)
                        #        dpn = dp/ddd

                        #        M=3
                        #        plt.figure().gca(projection='3d')
                        #        X  = copy.copy(Wori[0,i-M/2:i+M])
                        #        Y  = copy.copy(Wori[1,i-M/2:i+M])
                        #        TX = copy.copy(Wori[3,i-M/2:i+M])
                        #        plt.plot(X,Y,'-ok',linewidth=5,markersize=5)

                        #        for theta in np.linspace(-pi/4,pi/4,10):


                        #                dp[0:3] = np.dot(Rz(theta),dpn[0:3])

                        #                #print dp[0:3]
                        #                #Wori[:,i+1] = p + ds*dp
                        #                #Wori[:,i+2] = p + 2*ds*dp
                        #                #dWori[:,i+2] = Wori[:,i+2]-Wori[:,i+1]
                        #                #for k in range(i+1,i+M):
                        #                k = i+1
                        #                pnext = copy.copy(Wori[:,k])
                        #                Wori[:,k] = p + ds*dp
                        #                dWori[:,k] = dp

                        #                #for eta in np.linspace(-pi/2,pi/2,10):
                        #                #        k = i+2

                        #                #        dpeta = copy.copy(dp)
                        #                #        dpeta[0:3] = np.dot(Rz(eta),dp[0:3])

                        #                #        pnext = copy.copy(Wori[:,k])
                        #                #        Wori[:,k] = Wori[:,k-1] + ds*dpeta
                        #                #        dWori[:,k] = dpeta

                        #                X  = copy.copy(Wori[0,i:i+M])
                        #                Y  = copy.copy(Wori[1,i:i+M])
                        #                TX = copy.copy(Wori[3,i:i+M])

                        #                [cp_semin, cp_semax] = traj.GetSpeedIntervalAtCriticalPoint(env, Wori, dWori, critical_pt+M+1)
                        #                print "ctr",critical_pt+M,"vel:[",cp_semin,",",cp_semax,"]"
                        #                CP2 = traj.getCriticalPointFromWaypoints(env, Wori, dWori, None)
                        #                print "[CP2]:",CP2


                        #                if cp_semax > 0.01:
                        #                        plt.plot(X,Y,'-og',linewidth=2)
                        #                else:
                        #                        plt.plot(X,Y,'-or',linewidth=2)
                                ##dWori[0:3,i+1] = np.dot(Rz(theta),dp[0:3])#copy.copy(qnext-p)
                                ##dWori[:,i+1] = copy.copy(dqnext)
                                #[cp_semin, cp_semax] = traj.GetSpeedIntervalAtCriticalPoint(env, Wori, dWori, i)
                                ##CP2 = traj.getCriticalPointFromWaypoints(env, Wori, dWori, None)
                                #print "[CP2]:",CP2
                                #plt.show()
                                #sys.exit(0)
                                #smax = cp_semax
                                #print "ctr",i,"vel:[",cp_semin,",",cp_semax,"]","MVC"
                                #print "ds",np.linalg.norm(qnext-p)
                                #print "####################################################"
                        ##raw_input('Press <ENTER> to continue.')
                        #params.VisualizeReachableSet3D(p, dp, dp, smax, ds, Fcur)
                        #print "Projectability",i-critical_pt,"/",Mforward,Wmove[:,i+1]
                        #if i==critical_pt+int(Mforward/2):
                        #        PrintNumpy('p', p)
                        #        PrintNumpy('dp', dp)
                        #        PrintNumpy('force', Fcur)
                        #        print "speed=",smax
                        #        params.VisualizeReachableSet3D(p, dp, Wmove[:,i+1], smax, ds, Fcur)

                        #Wnext[:,i+1]= qnext
                        #dWnext[:,i+1]= dqnext
                        #[cp_semin, cp_semax] = traj.GetSpeedIntervalAtCriticalPoint(env, Wnext, dWnext, i)
                        #Wori[:,i+1] = copy.copy(qnext)
                        #dWori[:,i+1] = copy.copy(dqnext)
                        #[cp_semin, cp_semax] = traj.GetSpeedIntervalAtCriticalPoint(env, Wori, dWori, i)
                        #CP2 = traj.getCriticalPointFromWaypoints(env, Wori, dWori, None)
                        #CP2 = traj.getCriticalPointFromWaypoints(env, Wori, dWori, None)
                        #print "SPEED",cp_semin,cp_semax
                        #smax = cp_semax
                        #print "[CP2]:",CP2

                        #break
                        Qk.append(qnext)
                        p = copy.copy(qnext)
                        dp = copy.copy(dqnext)
                        i=i+1

                #Qk = np.array(Qk).T
                #print Wk.shape #2,10
                #print Qk.shape #2,11
                #plt.figure().gca(projection='3d')
                #plt.plot(Qk[0,:],Qk[1,:],Qk[3,:],'-or',linewidth=2)
                #plt.plot(Wk[0,:],Wk[1,:],Wk[3,:],'-ok',linewidth=2)
                #plt.show()
                #sys.exit(0)

                print Wmove[:,critical_pt:critical_pt+Mforward]
                #print dUtmp[:,critical_pt:critical_pt+Mforward]
                sf = 5.0
                for i in range(0,Nwaypoints):
                        if i>critical_pt:
                               #A = self.SmoothVector(traj,i,Wori,smoothing_factor=sf)
                                A = self.SmoothVectorProjection(traj,critical_pt,i,Wori,smoothing_factor=sf)
                                dUtmp[:,i] = np.dot(A, (Wmove.T))
                        #dUtmp[:,i] = Wmove[:,i]#np.dot(A, (lambda_coeff * Wmove.T))

                #print dUtmp
                #sys.exit(0)
                return dUtmp

        def get_name(self):
                return "reachable-set projection"

        def SmoothVectorProjection(self, traj, Ncritical, Ncur, W, smoothing_factor=None):
                if smoothing_factor is None:
                        smoothing_factor = self.DeformInfo['smoothing_factor']

                [Ndim, Nwaypoints] = traj.getWaypointDim(W)
                A = np.zeros(Nwaypoints)

                assert(Ncritical<Nwaypoints)
                i = Nwaypoints-1

                while i > 0:
                        A[i] = self.avalue(Ncur, i, smoothing_factor)
                        #if i > Ncur:
                        #else:
                        #        A[i] = 0.0
                        i -= 1
                return A

