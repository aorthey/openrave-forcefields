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
                        A1[i] = avalue(Ncritical, i, k*20.0)
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

def AReachMatrix(traj, Ncritical, W):
        [Ndim, Nwaypoints] = traj.getWaypointDim(W)
        A = np.zeros(Nwaypoints)

        assert(Ncritical<Nwaypoints)
        i = Nwaypoints-1
        while i > 0:
                if i>Ncritical:
                        A[i] = avalue(Ncritical, i, k*10.0)
                else:
                        A[i] = 2.0*(Ncritical-i)
                        #A[i] = 2.0**(Ncritical-i)
                i -= 1
        return A
def A3matrix(traj, Ncritical, W):
        [Ndim, Nwaypoints] = traj.getWaypointDim(W)
        A3 = np.zeros(Nwaypoints)

        assert(Ncritical<Nwaypoints)

        #M = 100
        i = Nwaypoints-1
        while i > 0:
                #if abs(i-Ncritical)<M and i<Nwaypoints:
                if i>Ncritical:
                        #A3[i] = avalue(Ncritical, i, k*20.0)
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
                A4[i] = avalue(Nwaypoints-1, i, k*10.0)
                i -= 1
        return A4

def A5matrix(traj, W):
        [Ndim, Nwaypoints] = traj.getWaypointDim(W)
        A5 = np.zeros(Nwaypoints)

        #M = 100
        i = Nwaypoints-1
        while i > 0:
                A5[i] = avalue(0, i, k*40.0)
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
        #lambda_1 = 0.001
        #lambda_1 = 0.0005
        lambda_1 = 0.000
        lambda_2 = 1e-3
        lambda_3 = 0.2

        def deform_onestep(self, computeNewCriticalPoint = True):
                COLLISION_ENABLED = True
                eta = 1.0

                traj = self.traj_deformed
                L = traj.get_length()

                ### traj.DISCRETIZATION_TIME_STEP controls ds
                [Wori,dWori,ddWori] = traj.get_waypoints_second_order()
                [Ndim, Nwaypoints] = traj.getWaypointDim(Wori)
                F = traj.get_forces_at_waypoints(Wori, self.env)
                [R,amin,amax] = traj.getControlMatrix(Wori)

                ###############################################################
                ## check if path dynamically feasible => return if yes
                ###############################################################

                self.critical_pt = traj.getCriticalPointFromWaypoints(self.env, Wori, dWori, ddWori, self.critical_pt)
                print "###########################################"
                print "CRITICAL WAYPOINT: ",self.critical_pt,"/",Nwaypoints


                if self.critical_pt >= Nwaypoints:
                        print "No deformation necessary => Trajectory dynamically feasible"
                        print traj.getCriticalPointFromWaypoints(self.env, Wori, dWori, ddWori, self.critical_pt)
                        print "###########################################"
                        #traj.PlotParametrization(self.env)
                        return DEFORM_NONE
                else:
                        ##plot critical pt
                        Nc_handle = self.env.env.plot3(points=Wori[0:3,self.critical_pt],
                                        pointsize=0.02,
                                        colors=np.array(((1.0,0.2,0.2))),
                                        drawstyle=1)
                print "###########################################"

                ###############################################################
                ## preliminary computations
                ###############################################################

                ### get forces in normal direction to trajectory
                FNxy = self.getForceNormalComponent(F, dWori)


                self.traj_velprofile = traj.getVelocityIntervalWithoutForceField(self.env, Wori, dWori, ddWori)
                Tend = self.traj_velprofile.duration
                Tstart = 0.0
                Tstep = Tend/1e3

                dpmin = np.zeros((1,Nwaypoints))
                dpmax = np.zeros((1,Nwaypoints))
                for i in range(0,Nwaypoints):
                        Tcur =Tstart
                        p = Wori[:,i]
                        q = self.traj_velprofile.Eval(Tcur)
                        while np.linalg.norm(p-q)>2*Tstep:
                                Tcur += Tstep
                                q = self.traj_velprofile.Eval(Tcur)
                        dq = self.traj_velprofile.Evald(Tcur)
                        dpmax[:,i] = np.linalg.norm(dq)
                        Tstart = Tcur

                #Tz = self.getTorqueComponent(F, Wori, dWori)
                lambda_1 = self.lambda_1
                lambda_2 = self.lambda_2
                lambda_3 = self.lambda_3
                print "## LAMBDAS: ",lambda_1,lambda_2

                dU = np.zeros((Ndim,Nwaypoints))

                #################################################################
                ## lambda 1 update
                ## move against the normal component of the force field
                #################################################################
                dUtmp = np.zeros((Ndim,Nwaypoints))
                for i in range(0,Nwaypoints):
                        B1 = B1matrix(traj,i,Wori)
                        dUtmp[:,i] += np.dot(B1,( -lambda_1 * FNxy.T))

                Wnext = Wori + eta*dUtmp
                if COLLISION_ENABLED:
                        if not self.traj_deformed.IsInCollision(self.env, Wnext):
                                dU += dUtmp
                        else:
                                print "## $> lambda1 collision (contra-force movement)"
                else:
                        dU += dUtmp

                print "## LAMBDAS: ",lambda_1,lambda_2
                #################################################################
                ## lambda 2 update
                ## project onto reachable set
                #################################################################
                Wmove = np.zeros((Ndim,Nwaypoints))
                dUtmp = np.zeros((Ndim,Nwaypoints))

                ds = traj.DISCRETIZATION_TIME_STEP
                #### compute the intersection of a forward simulation and the
                #### ball with radius ds

                for i in range(0,Nwaypoints-1):
                        p = Wori[:,i]
                        dp = dWori[:,i]
                        pnext = Wori[:,i+1]
                        #smax = dpmax[:,i]/2
                        smax = 0

                        if np.linalg.norm(F[:,i])>1e-3:
                                qnext = np.zeros((Ndim))
                                qnext = p
                                tstep = 1e-4
                                dt = 0.0
                                while abs(np.linalg.norm(p-qnext) - ds) > 4*tstep:
                                        dt += tstep
                                        dt2 = dt*dt/2
                                        qnext = p + dt*smax*dp + dt2*F[:,i]

                                #print "found qnext at :",qnext-p,"time",dt
                                #print "ds:",ds, \
                                        #"d(p,pnext):",np.linalg.norm(p-pnext), \
                                        #"d(qnext,pnext):",np.linalg.norm(pnext-qnext)
                                #sys.exit(0)
                                dpq = np.dot(np.linalg.norm(qnext),np.linalg.norm(pnext))
                                if dpq < 0.1:
                                        print "WaRNING: dpq:",dpq
                                        sys.exit(0)

                                ## project onto orthogonal direction to path
                                dw = qnext - pnext
                                #dw = dw - np.dot(dw,pnext-p)*(pnext-p)
                                Wmove[:,i+1] = dw


                for i in range(0,Nwaypoints):
                        A = AReachMatrix(traj,i,Wori)
                        #dUtmp[:,i+1] = np.dot(A3, (lambda_2 * Wmove.T))
                        #if i==int(Nwaypoints/2):
                                #print i,np.around(A,decimals=2)
                                #sys.exit(0)
                        dUtmp[:,i] = np.dot(A, (lambda_2 * Wmove.T))

                Wnext = Wori + eta*dUtmp

                if COLLISION_ENABLED:
                        if not self.traj_deformed.IsInCollision(self.env, Wnext):
                                dU += dUtmp
                        else:
                                print "## $> lambda2 collision (reachable-set projection)"
                else:
                        dU += dUtmp
                print "## LAMBDAS: ",lambda_1,lambda_2
                #################################################################
                ## lambda 3 update
                ## orient into path direction
                #################################################################
                Tdir = np.zeros((1,Nwaypoints))
                dUtmp = np.zeros((Ndim,Nwaypoints))
                for i in range(0,Nwaypoints):
                        theta = Wori[3,i]

                        dqxy = dWori[0:2,i]/np.linalg.norm(dWori[0:2,i])
                        thetapath = acos(np.dot(dqxy,ex[0:2]))

                        d1 = abs(theta-thetapath)
                        d2 = 2*pi-abs(theta-thetapath)

                        thetanew = min(d1,d2)
                        Tdir[:,i] = thetanew

                for i in range(0,Nwaypoints):
                        A1 = A1matrix(traj,i,Wori)
                        dUtmp[3,i] = np.dot(A1, (lambda_3 * Tdir.T))

                Wnext = Wori + eta*dUtmp

                if COLLISION_ENABLED:
                        if not self.traj_deformed.IsInCollision(self.env, Wnext):
                                dU += dUtmp
                        else:
                                print "## $> lambda3 collision (orientation)"
                else:
                        dU += dUtmp

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

                if COLLISION_ENABLED:
                        if self.traj_deformed.IsInCollision(self.env, Wnext):
                                print "no deformation possible -> collision"
                                return DEFORM_COLLISION

                self.traj_deformed.new_from_waypoints(Wnext)
                return DEFORM_OK
