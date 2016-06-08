from deformation_factory import *
from util_mvc import *

import sys

DEBUG = 0

def avalue(Ncritical, i, c=10.0):
        return np.exp(-((Ncritical-i)*(Ncritical-i))/(2*c*c))

def A1matrix(traj, Ncritical, W):
        [Ndim, Nwaypoints] = traj.getWaypointDim(W)
        A1 = np.zeros(Nwaypoints)

        assert(Ncritical<Nwaypoints)

        i = Nwaypoints
        while i > 0:
                #if abs(i-Ncritical)<M and i<Nwaypoints:
                if i<Nwaypoints-1:
                        A1[i] = avalue(Ncritical, i,15.0)
                i -= 1
        return A1
def AOrientationMatrix(traj, Ncritical, W):
        [Ndim, Nwaypoints] = traj.getWaypointDim(W)
        A1 = np.zeros(Nwaypoints)

        assert(Ncritical<Nwaypoints)

        i = Nwaypoints
        while i > 0:
                #if abs(i-Ncritical)<M and i<Nwaypoints:
                if i<Nwaypoints-1:
                        A1[i] = avalue(Ncritical, i,15.0)
                i -= 1
        return A1

def ABackwardMatrix(traj, Ncritical, W):
        [Ndim, Nwaypoints] = traj.getWaypointDim(W)
        A = np.zeros(Nwaypoints)

        assert(Ncritical<Nwaypoints)
        i = Nwaypoints-1
        while i > 0:
                if i>Ncritical:
                        #A[i] = avalue(Ncritical, i, k*20.0)
                        A[i] = (i-Ncritical+1)*traj.DISCRETIZATION_TIME_STEP
                else:
                        A[i] = 0.0
                        #A[i] = 2.0*(Ncritical-i)
                        #A[i] = 2.0**(Ncritical-i)
                i -= 1
        return A
def AForwardMatrix(traj, Ncritical, W):
        [Ndim, Nwaypoints] = traj.getWaypointDim(W)
        A = np.zeros(Nwaypoints)

        assert(Ncritical<Nwaypoints)
        i = Nwaypoints-1
        while i > 0:
                if i>Ncritical:
                        #A[i] = avalue(Ncritical, i, k*20.0)
                        A[i] = 0.0
                else:
                        #A[i] = 2.0*(Ncritical-i)
                        A[i] = (Ncritical-i+1)*traj.DISCRETIZATION_TIME_STEP
                        #A[i] = 2.0**(Ncritical-i)
                i -= 1
        return A
def AEndMatrix(traj, W):
        [Ndim, Nwaypoints] = traj.getWaypointDim(W)
        A4 = np.zeros(Nwaypoints)

        #M = 100
        i = Nwaypoints-1
        while i >= 0:
                #A4[i] = avalue(Nwaypoints-1, i, 10.0)
                A4[i] = avalue(Nwaypoints-1, i, 50.0)
                i -= 1
        return A4

def AStartMatrix(traj, W):
        [Ndim, Nwaypoints] = traj.getWaypointDim(W)
        A5 = np.zeros(Nwaypoints)

        #M = 100
        i = Nwaypoints-1
        while i >= 0:
                A5[i] = avalue(0, i, 50.0)
                i -= 1
        return A5


class DeformationReachableSet(Deformation):

        def getForceNormalComponent(self, F, dWori):
                Nwaypoints = dWori.shape[1]
                FN = np.zeros((F.shape))
                for i in range(0,Nwaypoints):
                        trajNormal = np.array((-dWori[1,i],dWori[0,i],0,0))
                        trajNormal = trajNormal/np.linalg.norm(trajNormal)
                        FN[:,i] = np.dot(F[:,i].flatten().T, trajNormal)*trajNormal
                return FN

        def getTorqueNormalComponent(self, F, Wori, dWori):
                Nwaypoints = Wori.shape[1]
                FNtorque = np.zeros((F.shape))
                for i in range(0,Nwaypoints):
                        #trajNormal = np.array((-dWori[1,i],dWori[0,i],0,0))
                        #trajNormal = trajNormal/np.linalg.norm(trajNormal)
                        FNtorque[3,i] = F[3,i]
                        #FNtorque[:,i] = np.dot(F[:,i].flatten().T, trajNormal)*trajNormal

                return FNtorque

        def ComputeNextWaypoints(self, Wori, eta, dUtmp):
                Wnext = Wori + eta*dUtmp
                return Wnext

        ## change only traj_deformed here
        Nc_handle = []

        ###############################################################
        ### LAMBDA1: (smooth) move against force
        ### LAMBDA2: (smooth) project onto reachable set
        ###############################################################
        #lambda_1 = 0.001
        #lambda_1 = 0.0005
        lambda_1 = 0.0005
        lambda_2 = 1
        lambda_3 = 0.5*1e-2

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
                        return DEFORM_SUCCESS
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
                FNtorque = self.getTorqueNormalComponent(F, Wori, dWori)

                #### compute min/max velocity profile from path without forces
                #### (if available). otherwise use [0,0]

                self.traj_velprofile = traj.getVelocityIntervalWithoutForceField(self.env, Wori, dWori, ddWori)
                dpmin = np.zeros((1,Nwaypoints))
                dpmax = np.zeros((1,Nwaypoints))
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


                #Tz = self.getTorqueComponent(F, Wori, dWori)
                print "## LAMBDAS: ",self.lambda_1,self.lambda_2,self.lambda_3

                dU = np.zeros((Ndim,Nwaypoints))

                DeformInfo = {}
                DeformInfo['Ndim'] = Ndim
                DeformInfo['Nwaypoints'] = Nwaypoints
                DeformInfo['lambda_1'] = self.lambda_1
                DeformInfo['lambda_2'] = self.lambda_2
                DeformInfo['lambda_3'] = self.lambda_3
                DeformInfo['traj'] = self.traj_deformed
                DeformInfo['Wori'] = Wori
                DeformInfo['dWori'] = dWori
                DeformInfo['F'] = F
                DeformInfo['R'] = R
                DeformInfo['amin'] = amin
                DeformInfo['amax'] = amax
                DeformInfo['eta'] = eta
                DeformInfo['env'] = self.env

                #################################################################
                ## lambda 1 update
                ## move against the normal component of the force field
                #################################################################

                from deformation_module_counterwrench import *
                d1 = DeformationModuleCounterWrench( DeformInfo )
                dU += d1.get_update( self.lambda_1 )

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
                        smax = dpmax[:,i]/2
                        #smax = 0

                        if np.linalg.norm(F[:,i])>1e-3:
                                qnext = np.zeros((Ndim))
                                qnext = p
                                tstep = 1e-4
                                dt = 0.0

                                dold = 1e5
                                dnew = abs(np.linalg.norm(p-qnext) - ds)
                                while dnew < dold:
                                        dold = dnew
                                        dt += tstep
                                        dt2 = dt*dt/2
                                        qnext = p + dt*smax*dp + dt2*F[:,i]
                                        dnew = abs(np.linalg.norm(p-qnext) - ds)

                                dpq = np.dot(np.linalg.norm(qnext),np.linalg.norm(pnext))
                                if dpq < 0.1:
                                        print "WaRNING: dpq:",dpq
                                        sys.exit(0)

                                ## project onto orthogonal direction to path
                                dp = qnext - pnext
                                dp = dp - np.dot(dp,pnext-p)*(pnext-p)
                                Wmove[:,i] = dp
                                #Wmove[3,i] *= 5


                for i in range(0,Nwaypoints):
                        A = AForwardMatrix(traj,i,Wori)
                        dUtmp[:,i] = np.dot(A, (self.lambda_2 * Wmove.T))

                Wnext = self.ComputeNextWaypoints(Wori, eta, dUtmp)

                if COLLISION_ENABLED:
                        if not self.traj_deformed.IsInCollision(self.env, Wnext):
                                dU += dUtmp
                        else:
                                print "## $> lambda2 collision (reachable-set projection)"
                                ### make lambda smaller to see if step size is
                                ### an issue
                                self.lambda_2 /= 2.0
                else:
                        dU += dUtmp

                #################################################################
                ## lambda 3 update
                ## orient into path direction
                #################################################################
                Tdir = np.zeros((1,Nwaypoints))
                dUtmp = np.zeros((Ndim,Nwaypoints))

                for i in range(0,Nwaypoints-1):
                        p = Wori[:,i]
                        theta = p[3]
                        dp = dWori[:,i]
                        pnext = Wori[:,i+1]
                        smax = dpmax[:,i]/2

                        if np.linalg.norm(F[:,i])>1e-3:
                                qnext = np.zeros((Ndim))
                                qnext = p
                                tstep = 1e-4
                                dt = 0.0

                                dold = 1e5
                                dnew = abs(np.linalg.norm(p-qnext) - ds)
                                while dnew < dold:
                                        dold = dnew
                                        dt += tstep
                                        dt2 = dt*dt/2
                                        qnext = p + dt*smax*dp + dt2*F[:,i]
                                        dnew = abs(np.linalg.norm(p-qnext) - ds)


                                pori = np.zeros((Ndim))
                                dline = np.zeros((Ndim))
                                pori[0:3] = p[0:3]

                                dline[0:3] = np.dot(Rz(pi/2),(dWori[:,i])[0:3])
                                theta_step = 1e-2

                                pori[3] = p[3]+theta_step
                                [R,atmp,a2tmp] = self.traj_deformed.getControlMatrix(pori)
                                R = R[:,:,0]
                                #qnext1 = pori + dt*smax*dp + dt2*F[:,i]
                                qproj_f1 = np.dot(np.dot(Rz(pori[3]),ex),FNxy[0:3,i])
                                vol1 = self.GetReachableSetProjectedVolume( dt, ds, qnext, dline, pori, smax, dp, F[:,i], R,amin, amax)

                                pori[3] = p[3]-theta_step
                                [R,atmp,a2tmp] = self.traj_deformed.getControlMatrix(pori)
                                R = R[:,:,0]
                                #qnext2 = pori + dt*smax*dp + dt2*F[:,i]
                                qproj_f2 = np.dot(np.dot(Rz(pori[3]),ex),FNxy[0:3,i])
                                vol2 = self.GetReachableSetProjectedVolume( dt, ds, qnext, dline, pori, smax, dp, F[:,i], R,amin, amax)

                                pori[3] = p[3]
                                [R,atmp,a2tmp] = self.traj_deformed.getControlMatrix(pori)
                                R = R[:,:,0]
                                vol3 = self.GetReachableSetProjectedVolume( dt, ds, qnext, dline, pori, smax, dp, F[:,i], R,amin, amax)

                                if vol1 <= vol3 and vol2 <= vol3:
                                        Tdir[:,i]=0
                                else:
                                        if vol1 > vol3 and vol2 > vol3:
                                                ## both directions are good, so
                                                ## go against force field to
                                                ## increase reachability
                                                if qproj_f1 > qproj_f2:
                                                        ## move to v2
                                                        Tdir[:,i] = -1
                                                else:
                                                        Tdir[:,i] = 1

                                        elif vol1 > vol2:
                                                ##move into + direction
                                                Tdir[:,i] = 1
                                        else:
                                                Tdir[:,i] = -1

                for i in range(0,Nwaypoints):
                        A1 = AOrientationMatrix(traj,i,Wori)
                        dUtmp[3,i] = np.dot(A1, (self.lambda_3*Tdir.T))

                Wnext = self.ComputeNextWaypoints(Wori, eta, dUtmp)

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
                A4 = AEndMatrix(traj, Wori)
                A5 = AStartMatrix(traj, Wori)
                for i in range(0,Nwaypoints):
                        dU[:,i] += -A4[i]*dEnd
                        dU[:,i] += -A5[i]*dStart

                Wnext = Wori + eta*dU
                #Wnext = self.ComputeNextWaypoints(Wori, eta, dUtmp)

                if not (dU[:,0]==0).all():
                        print "ERROR: start point deformation -> not allowed"
                        sys.exit(1)
                if not (dU[:,-1]==0).all():
                        print "ERROR: end point deformation -> not allowed"
                        sys.exit(1)

                if np.linalg.norm(Wori-Wnext)<1e-10:
                        print "no deformation achieved with current critical point"
                        return DEFORM_NOPROGRESS

                if COLLISION_ENABLED:
                        if self.traj_deformed.IsInCollision(self.env, Wnext):
                                print "no deformation possible -> collision"
                                return DEFORM_COLLISION

                self.traj_deformed.new_from_waypoints(Wnext)
                return DEFORM_OK

        def GetReachableSetVerticesSE2(self, dt, p, s, dp, F, R, amin, amax):
                Ndim = p.shape[0]
                dt2 = dt*dt*0.5
                a = np.zeros((amin.shape))

                rs_vertices = np.zeros((Ndim,4))

                ## create ordered set (clockwise)
                a[0] = amin[0]
                a[1] = amin[1]
                rs_vertices[:,0] = p + dt*s*dp + dt2*F + dt2*np.dot(R,a)
                a[0] = amin[0]
                a[1] = amax[1]
                rs_vertices[:,1] = p + dt*s*dp + dt2*F + dt2*np.dot(R,a)
                a[0] = amax[0]
                a[1] = amax[1]
                rs_vertices[:,2] = p + dt*s*dp + dt2*F + dt2*np.dot(R,a)
                a[0] = amax[0]
                a[1] = amin[1]
                rs_vertices[:,3] = p + dt*s*dp + dt2*F + dt2*np.dot(R,a)

                return rs_vertices

        def GetReachableSetProjectedVolume( self, dt, ds, pline, dline, p, smax, dp, F, R, amin, amax):

                rs_vertices = self.GetReachableSetVerticesSE2( dt, p, smax, dp, F, R, amin, amax)
                Nvertices = rs_vertices.shape[1]

                ##project onto ds-ball around p

                rs_vertices_proj = np.zeros(rs_vertices.shape)

                dline = dline/np.linalg.norm(dline)

                dp = np.zeros(Nvertices)

                for i in range(0,Nvertices):
                        pvertex = rs_vertices[:,i]-pline
                        dp[i] = np.dot(pvertex[0:2], dline[0:2])

                mindp = np.min(dp)
                maxdp = np.max(dp)
                return maxdp - mindp
