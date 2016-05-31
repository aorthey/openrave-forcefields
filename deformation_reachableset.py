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

        try:
                L = dq/np.linalg.norm(dq)
                Q = (qt-q)/np.linalg.norm(qt-q)
                dtravel = np.dot((qt-q),L)
                dangle = acos(np.dot(L,Q))
                        
        except ValueError as e:
                print "#######################################"
                print "ValueError:",e
                print "#######################################"
                print "Q    = ",Q
                print "qt   = ",qt
                print "q    = ",q
                print "L    = ", L
                print "qt-q = ",qt-q
                sys.exit(1)

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


k = 0.5
def avalue(Ncritical, i, c=k*20.0):
        return np.exp(-((Ncritical-i)*(Ncritical-i))/(2*c*c))

def B1matrix(traj, Ncritical, W):
        [Ndim, Nwaypoints] = traj.getWaypointDim(W)
        A1 = np.zeros(Nwaypoints)

        assert(Ncritical<Nwaypoints)

        i = Nwaypoints
        while i > 0:
                #if abs(i-Ncritical)<M and i<Nwaypoints:
                if i<Nwaypoints-1:
                        A1[i] = avalue(Ncritical, i,k*30.0)
                i -= 1
        return A1

def A1matrix(traj, Ncritical, W):
        [Ndim, Nwaypoints] = traj.getWaypointDim(W)
        A1 = np.zeros(Nwaypoints)

        assert(Ncritical<Nwaypoints)

        i = Nwaypoints
        while i > 0:
                #if abs(i-Ncritical)<M and i<Nwaypoints:
                if i<Nwaypoints-1:
                        A1[i] = avalue(Ncritical, i, k*10.0)
                i -= 1
        return A1

def A2matrix(traj, Ncritical, W):
        [Ndim, Nwaypoints] = traj.getWaypointDim(W)
        A2 = np.zeros(Nwaypoints)
        assert(Ncritical<Nwaypoints)

        i = Ncritical-1
        while i > 0:
                A2[i] = avalue(Ncritical, i)
                i -= 1
        return A2

def A3matrix(traj, Ncritical, W):
        [Ndim, Nwaypoints] = traj.getWaypointDim(W)
        A3 = np.zeros(Nwaypoints)

        assert(Ncritical<Nwaypoints)

        #M = 100
        i = Nwaypoints-1
        while i > 0:
                #if abs(i-Ncritical)<M and i<Nwaypoints:
                if i>Ncritical:
                        A3[i] = avalue(Ncritical, i, k*10.0)
                else:
                        A3[i] = 1.0
                        #A3[i] = 1.0+0.1*(Ncritical-i)
                i -= 1
        return A3
def A4matrix(traj, W):
        [Ndim, Nwaypoints] = traj.getWaypointDim(W)
        A4 = np.zeros(Nwaypoints)

        #M = 100
        i = Nwaypoints-1
        while i > 0:
                A4[i] = avalue(Nwaypoints-1, i, k*20.0)
                i -= 1
        return A4

def A5matrix(traj, W):
        [Ndim, Nwaypoints] = traj.getWaypointDim(W)
        A5 = np.zeros(Nwaypoints)

        #M = 100
        i = Nwaypoints-1
        while i > 0:
                A5[i] = avalue(0, i, k*30.0)
                i -= 1
        return A5

DEBUG=0

class DeformationReachableSet(Deformation):

        def getForceNormalComponent(self, F, dWori):
                Nwaypoints = dWori.shape[1]
                FN = np.zeros((F.shape))
                for i in range(0,Nwaypoints):
                        trajNormal = np.array((-dWori[1,i],dWori[0,i],0,0))
                        trajNormal = trajNormal/np.linalg.norm(trajNormal)
                        FN[:,i] = np.dot(F[:,i].flatten().T, trajNormal)*trajNormal
                return FN

        def getTorqueComponent(self, F, Wori, dWori):
                Nwaypoints = Wori.shape[1]
                FN = np.zeros((F.shape))
                for i in range(0,Nwaypoints):
                        trajNormal = np.array((-dWori[1,i],dWori[0,i],0,0))
                        trajNormal = trajNormal/np.linalg.norm(trajNormal)
                        FN[:,i] = np.dot(F[:,i].flatten().T, trajNormal)*trajNormal
                        df = np.dot(F[:,i].flatten().T, trajNormal)
                        #if df > 0:

                return FN




        ## change only traj_deformed here
        Nc_handle = []

        ###############################################################
        ### LAMBDA1: (smooth) move against force
        ### LAMBDA2: (smooth) project onto reachable set
        ###############################################################
        lambda_1 = 0.001
        lambda_2 = 0.00

        def deform_onestep(self, computeNewCriticalPoint = True):
                eta = 1.0

                traj = self.traj_deformed
                L = traj.get_length()
                [Wori,dWori,ddWori] = traj.get_waypoints_second_order()
                [Ndim, Nwaypoints] = traj.getWaypointDim(Wori)
                F = traj.get_forces_at_waypoints(Wori, self.env)
                [R,amin,amax] = traj.getControlMatrix(Wori)

                ### FORWARD PASS UNTIL CRITICAL POINT IS HIT
                if computeNewCriticalPoint:
                        Nc = traj.getCriticalPointFromWaypoints(self.env, Wori, dWori, ddWori, self.critical_pt)
                else:
                        Nc = self.critical_pt

                print "###########################################"
                print "CRITICAL WAYPOINT: ",Nc,"/",Nwaypoints," oldNc=",self.critical_pt
                self.critical_pt = Nc
                #self.Nc_handle=self.env.env.plot3(Wori[0:3,Nc]+np.array((0,0,0.1)),pointsize=0.07,colors=np.array((1.0,1.0,0)),drawstyle=1)

                if Nc >= Nwaypoints:
                        print "No deformation necessary => Trajectory dynamically feasible"
                        print traj.getCriticalPointFromWaypoints(self.env, Wori, dWori, ddWori, self.critical_pt)
                        print "###########################################"
                        traj.PlotParametrization(self.env)
                        return DEFORM_NONE

                traj = self.traj_deformed
                L = traj.get_length()
                [Wori,dWori,ddWori] = traj.get_waypoints_second_order()
                [Ndim, Nwaypoints] = traj.getWaypointDim(Wori)
                F = traj.get_forces_at_waypoints(Wori, self.env)
                [R,amin,amax] = traj.getControlMatrix(Wori)

                #print "POS   :",Wori[0:2,Nc]
                #print "FORCE :",F[0:2,Nc]
                print "###########################################"
                ###############################################################
                ## for all points where dF = 0, we can say that lambda_1 = lambda_2 = 0
                ###############################################################

                ###############################################################
                ## update trajectory into lambda directions
                ###############################################################
                q = Wori[:,Nc]
                self.FN_critical = F[:,Nc].flatten()

                FNxy = self.getForceNormalComponent(F, dWori)
                #Tz = self.getTorqueComponent(F, Wori, dWori)
                lambda_1 = self.lambda_1
                lambda_2 = self.lambda_2

                A1 = A1matrix(traj,Nc,Wori)
                A2 = A2matrix(traj,Nc,Wori)
                dU = np.zeros((Ndim,Nwaypoints))

                print "## LAMBDAS: ",lambda_1,lambda_2

                #################################################################
                ## lambda 1 update
                ## move against the normal component of the force field
                #################################################################
                dUtmp = np.zeros((Ndim,Nwaypoints))
                for i in range(0,Nwaypoints):
                        #dUtmp[:,i] += A1[i]*(-lambda_1 * F[:,Nc])
                        B1 = B1matrix(traj,i,Wori)
                        #dUtmp[:,i] += np.dot(B1,( -lambda_1 * FNxy.T))
                        dUtmp[:,i] += np.dot(B1,( -lambda_1 * FNxy.T))



                Wnext = Wori + eta*dUtmp
                if not self.traj_deformed.IsInCollision(self.env, Wnext):
                        dU += dUtmp
                else:
                        print "## $> lambda1 collision (contra-force movement)"

                #################################################################
                ## lambda 2 update
                ## project onto reachable set
                #################################################################
                dt = 0.001 ##same as TOPP
                dt2 = dt*dt/2
                Wmove = np.zeros((Ndim,Nwaypoints))

                for i in range(1,Nwaypoints):
                        Wori[:,i]
                        vn = dWori[:,i-1]/np.linalg.norm(dWori[:,i-1])
                        p = dt2*F[:,i-1]
                        ##new point: p+oproj
                        Wproj = p + np.dot(Wori[:,i]-p,vn)*vn
                        ## move toward oproj
                        Wmove[:,i] = Wproj - Wori[:,i]

                dUtmp = np.zeros((Ndim,Nwaypoints))
                for i in range(1,Nwaypoints):
                        A3 = A3matrix(traj,i,Wori)
                        #dUtmp[:,i] += np.dot(A3,( lambda_2 * Wmove))

                Wnext = Wori + eta*dUtmp

                if not self.traj_deformed.IsInCollision(self.env, Wnext):
                        dU += dUtmp
                else:
                        print "## $> lambda2 collision (reachable-set projection)"

                #################################################################
                ## projection onto fixed endpoint subspace
                #################################################################
                dEnd = dU[:,-1]
                dStart = dU[:,0]
                A4 = A4matrix(traj, Wori)
                A5 = A5matrix(traj, Wori)
                for i in range(0,Nwaypoints):
                        dU[:,i] += -A4[i]*dEnd
                        dU[:,i] += -A5[i]*dStart

                Wnext = Wori + eta*dU

                if np.linalg.norm(Wori-Wnext)<1e-10:
                        print "no deformation achieved with current critical point"
                        return DEFORM_NOPROGRESS

                if self.traj_deformed.IsInCollision(self.env, Wnext):
                        print "no deformation possible -> collision"
                        return DEFORM_COLLISION
                else:
                        self.traj_deformed.new_from_waypoints(Wnext)
                        return DEFORM_OK

