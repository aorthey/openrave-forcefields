from deformation_factory import *
from util_mvc import *

import sys

DEBUG = 0
def computeCostFunction(q, dq, ddq, F, amax, l1, l2, dt):
        dt2 = dt*dt*0.5
        trajNormal = np.array((-q[1],q[0],0.0,0.0))
        trajNormal = trajNormal/np.linalg.norm(trajNormal)
        Fproj = np.dot(F.flatten().T, trajNormal)
        Qlambda1 = -l1*Fproj*trajNormal
        Qlambda2 = sqrt(l2*amax[0]*2)*dt*dq
        qt = q + dt2*F.flatten() + Qlambda1 + Qlambda2

        L = dq/np.linalg.norm(dq)
        Q = (qt-q)/np.linalg.norm(qt-q)
        dtravel = np.dot((qt-q),L)
        dangle = acos(np.dot(L,Q))

        d = 1*dangle-0.1*dtravel
        d=dtravel
        return d


def compute_lambda_updates(q,dq,ddq,F,amax,dt):

        ### compute minimal direction in l1/l2 direction

        dold=10000+1
        d = 10000

        theta = 0.0
        thetaBest = 0.0
        tstep = 0.01
        epsilon = 0.05
        while d<dold:
                dold = d
                l1 = epsilon*cos(theta+tstep)
                l2 = epsilon*sin(theta+tstep)
                dplus = computeCostFunction( q, dq, ddq, F, amax, l1, l2, dt)
                l1 = epsilon*cos(theta-tstep)
                l2 = epsilon*sin(theta-tstep)
                dminus = computeCostFunction( q, dq, ddq, F, amax, l1, l2, dt)

                if dplus > dminus:
                        d = dminus
                        thetaBest = theta
                        theta = theta - tstep
                else:
                        d = dplus
                        thetaBest = theta
                        theta = theta + tstep

        return [epsilon*cos(thetaBest),epsilon*sin(thetaBest)]


def avalue(Ncritical, i):
        c = 5.0
        return np.exp(-((Ncritical-i)*(Ncritical-i))/(2*c*c))

def A1matrix(traj, Ncritical, W):
        [Ndim, Nwaypoints] = traj.getWaypointDim(W)
        A1 = np.zeros(Nwaypoints)

        assert(Ncritical<Nwaypoints)

        M = 20
        i = Nwaypoints
        while i > 0:
                if abs(i-Ncritical)<M and i<Nwaypoints:
                        A1[i] = avalue(Ncritical, i)
                i -= 1
        return A1

def A2matrix(traj, Ncritical, W):
        [Ndim, Nwaypoints] = traj.getWaypointDim(W)
        A2 = np.zeros(Nwaypoints)

        assert(Ncritical<Nwaypoints)

        M = 20
        i = Ncritical-1
        while i > 0:
                if abs(i-Ncritical-1)<M:
                        A2[i] = avalue(Ncritical, i)
                i -= 1
        return A2

class DeformationStretchPull(Deformation):

        ## change only traj_deformed here
        def deform_onestep(self):
                #dt = 0.01
                traj = self.traj_deformed
                L = traj.get_length()
                #Nwaypoints = int(L/dt)
                dt = 0.05
                Nwaypoints = int(L/dt)
                print "WAYPOINTS:",Nwaypoints,"LENGTH:",L
                [Wori,dWori,ddWori] = traj.get_waypoints_second_order(N=Nwaypoints)
                [Ndim, Nwaypoints] = traj.getWaypointDim(Wori)
                F = traj.get_forces_at_waypoints(Wori, self.env)
                [R,amin,amax] = traj.getControlMatrix(Wori)

                ### FORWARD PASS UNTIL CRITICAL POINT IS HIT
                Nc = traj.getCriticalPointFromWaypoints(self.env, Wori, dWori, ddWori)

                if Nc >= Nwaypoints:
                        ## no critical points => trajectory is dynamically
                        print "No deformation necessary=>Trajectory dynamically feasible"
                        ## feasible
                        return

                ### compute new waypoints
                #W = Wori[:,0:Nc+1]
                #dW = dWori[:,0:Nc+1]
                #ddW = ddWori[:,0:Nc+1]
                #F = traj.get_forces_at_waypoints(W, self.env)
                #[R,amin,amax] = traj.getControlMatrix(W)

                ###############################################################
                ## for all points where dF = 0, we can say that lambda_1 = lambda_2 = 0
                ###############################################################

                dt = 0.1
                [lambda_1, lambda_2] = compute_lambda_updates(Wori[:,Nc], dWori[:,Nc], ddWori[:,Nc], F[:,Nc], amax, dt)
                lambda_2=1

                ###############################################################
                ## update trajectory into lambda directions
                ###############################################################
                q = Wori[:,Nc]
                trajNormal = np.array((-q[1],q[0],0.0,0.0))
                trajNormal = trajNormal/np.linalg.norm(trajNormal)
                Fproj = np.dot(F[:,Nc].flatten().T, trajNormal)
                A1 = A1matrix(traj,Nc,Wori)
                A2 = A2matrix(traj,Nc,Wori)

                dU = np.zeros((Ndim,Nwaypoints))
                #map(lambda a: a*(-lambda_1*Fproj*trajNormal)
                #print np.around(A1,decimals=2)
                #print np.around(A2,decimals=2)
                #print A1[Nc],A2[Nc]
                for i in range(0,Nwaypoints):
                        dU[:,i] = A1[i]*(-lambda_1 * Fproj * trajNormal)
                        dU[:,i] += A2[i]*(- lambda_2 * dWori[:,Nc])

                #print dU
                #print lambda_1,lambda_2
                eta = 1.0
                print Wori[0:2,0]
                Wnext = Wori + eta*dU
                print Wori[0:2,0]
                self.traj_deformed.new_from_waypoints(Wnext)

