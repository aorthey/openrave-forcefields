import numpy as np
import sys
import copy
from numpy import sqrt,sin,cos,pi
from qcqp import *
import scipy
import cvxopt
import matplotlib.pyplot as plt
from cvxpy import *
from cvxopt import matrix, solvers
from util import PrintNumpy,inf

FILENAME = 'virus2'
DYNAMICAL_SYTEM_NAME = 'non-holonomic differential drive (deform)'
system_name = DYNAMICAL_SYTEM_NAME

AM = 0.1
DEBUG = True
### car/sailboat
#amin = np.array((-AM,-AM,-AM,0))
#amax = np.array((AM,AM,AM,0))
amin = np.array((-AM,-AM,-0.5*AM))
amax = np.array((AM,AM,0.5*AM))

def ControlPerWaypoint(W, Ndim, Nwaypoints):
        assert(Ndim==4)
        Kdim = 3
        R = np.zeros((Ndim,Kdim,Nwaypoints))
        for i in range(0,Nwaypoints):
                if Nwaypoints>1:
                        t = W[3,i]
                else:
                        t = W[3]

                R[0,:,i] = np.array((cos(t),-sin(t),0.0))
                R[1,:,i] = np.array((sin(t),cos(t),0.0))
                R[2,:,i] = np.array((0.0,0.0,0.0))
                R[3,:,i] = np.array((0.0,0.0,1.0))

        return R

### input: p -- waypoint on manifold
###     F -- force on manifold at p
###
### output: [A,b] such that A*qdd + b <= 0 

def GetControlConstraintMatrices(p, F):
        #Ndim = p.shape[0]
        #R = ControlPerWaypoint(p, Ndim, 1)[:,:,0]
        #return GetControlConstraintMatricesFromControl(R,F)
        return GetControlConstraintMatricesAdjust(p, F, epsilon=0)

def GetControlConstraintMatricesAdjust(p, F, epsilon=0.1):
        Ndim = p.shape[0]
        R = ControlPerWaypoint(p, Ndim, 1)[:,:,0]
        return GetControlConstraintMatricesFromControlAdjust(R,F,epsilon=epsilon)

def GetControlConstraintMatricesFromControl(R, F):
        return GetControlConstraintMatricesFromControlAdjust(R, F, epsilon=0)

def GetControlConstraintMatricesFromControlAdjust(R, F, epsilon=0):
        Adim = amin.shape[0]
        aminE = np.zeros((Adim))
        amaxE = np.zeros((Adim))
        for i in range(0,Adim):
                aminE[i]=amin[i]+epsilon
                amaxE[i]=amax[i]-epsilon
                if aminE[i]>amaxE[i]:
                        amid = amin[i]+0.5*(amax[i]-amin[i])
                        aminE[i] = amid
                        amaxE[i] = amid

        #print "amin",amin,"->",aminE
        #print "amax",amax,"->",amaxE
        Ndim = R.shape[0]
        Rmin = np.minimum(np.dot(R,aminE),np.dot(R,amaxE))
        Rmax = np.maximum(np.dot(R,aminE),np.dot(R,amaxE))
        H1 = F - Rmax
        H2 = -F + Rmin
        for j in range(Ndim):
                if H2[j] > -H1[j]:
                        print H2[j],"<= q[",j,"]<=",-H1[j]
                        sys.exit(1)
        b = np.hstack((H1,H2)).flatten()
        I = np.identity(Ndim)
        A = np.vstack((I,-I))
        return [A,b]

def GetNearestControlPoint(p, dp, pnext, F, dt, speed, Acontrol, bcontrol):
        ds_next = np.linalg.norm(p-pnext)
        dt2 = dt*dt/2
        Ndim = p.shape[0]

        if dt < 1e-10:
                qnext = p + dt*speed*dp + dt2*F
                dqnext = speed*dp + dt*F
                ds_next = np.linalg.norm(qnext-pnext)
                return [qnext,dqnext,ds_next]


        A = np.zeros((Acontrol.shape))
        b = np.zeros((bcontrol.shape))

        A = (-2*Acontrol)/(dt*dt)
        b = bcontrol + np.dot(A,(-p - dt*speed*dp))

###############################################################################
##### QCQP
###############################################################################
        Id = np.eye(Ndim)
        x = Variable(Ndim)
        objective = Minimize( norm(x - pnext))
        constraints = [ A*x <= -b,
                        quad_form(x-p,Id) <= ds_next]
        prob = Problem(objective, constraints)
        qcontrol = None
        try:
                #dnew = np.abs(prob.solve(solver=ECOS, verbose=False, max_iters=500, feastol_inacc=1e-20))
                dnew = np.abs(prob.solve(solver=ECOS, max_iters=500, feastol_inacc=1e-3))
                #dnew = np.abs(prob.solve(solver=SCS, max_iters=5000, eps=1e-200))
                #dnew = np.abs(prob.solve(solver=CVXOPT))
                #dnew = np.abs(prob.solve(solver=ECOS_BB, mi_rel_eps=1e-100))
                if dnew < inf:
                        qcontrol =np.array(x.value).flatten()
                        qdd = 2*(qcontrol-p-dt*speed*dp)/(dt*dt)
                        qdnext = speed*dp + dt*qdd
                        qnext = p + dt*speed*dp + dt2*qdd
                else:
                        qnext = None
                        qdnext = None
        except Exception as e:
                dnew = inf
                pass

###############################################################################
##### LP
###############################################################################
        #A = matrix(A)
        #b = matrix(b)
        #c = matrix( (pnext-p) )
        ##sol=solvers.lp(c,A,-b)
        #try:
        #        solvers.options['show_progress'] = False
        #        solvers.options['maxiters'] = 100 ##default: 100
        #        solvers.options['abstol'] = 1e-03 ## below 1e-11 some weird behavior
        #        solvers.options['reltol'] = 1e-03 ## below 1e-11 some weird behavior
        #        solvers.options['feastol'] = 1e-03 ## below 1e-11 some weird behavior
        #        sol=solvers.lp(c,A,-b)
        #        qcontrol = np.array(sol['x']).flatten()
        #if qcontrol is not None:

        #        qdd = 2*(qcontrol-p-dt*speed*dp)/(dt*dt)
        #        qdnext = speed*dp + dt*qdd
        #        qnext = p + dt*speed*dp + dt2*qdd

        #dnew = np.linalg.norm(qnext-pnext)

        #except Exception as e:
        #        print e
        #        sol=solvers.lp(c,A,-b)
        #        
        #        PrintNumpy('A', np.array(A))
        #        PrintNumpy('b', np.array(b))
        #        PrintNumpy('c', np.array(c))
        #        PrintNumpy('p', p)
        #        PrintNumpy('dp', dp)
        #        PrintNumpy('force', F)
        #        print "ds=",ds_next
        #        print "speed=",speed
        #        sys.exit(0)


###############################################################################
###############################################################################
        #if dnew < inf:
                #sys.exit(0)
        #print dnew,dt,speed,p,dp
        #for dt in self.expspace(tstart,tend,tsamples):
        #        dt2 = dt*dt*0.5
        #        M_speed = 1
        #        speed = np.linspace(0,s,M_speed)
        #        speed=[s]
        #        q = np.zeros((Ndim,M_speed*M_time))
        #        for k in range(0,M_speed):
        #                for i in range(0,M_time):
        #                        control = np.dot(R,A[i])
        #                        q[:,i+k*M_time] = p + dt*speed[k]*dp + dt2 * force + dt2 * control

        #if dnew < inf:
        #        qcontrol =np.array(x.value).flatten()
        #        qdd = 2*(qcontrol-p-dt*speed*dp)/(dt*dt)
        #        qdnext = speed*dp + dt*qdd
        #        qnext = p + dt*speed*dp + dt2*qdd
        #else:
        #        qnext = None
        #        qdnext = None

        return [qnext, qdnext, dnew]

def ForwardSimulate(p, dp, speed, ds, F, pnext=None):
        ### should best follow path!

        if pnext is None:
                pnext = p + ds*dp/np.linalg.norm(dp)

        ### estimate tstep
        dds = speed*np.linalg.norm(dp)
        tv = ds / (dds)
        tf = np.sqrt(2*ds/np.linalg.norm(F))
        tc = np.sqrt(2*ds/np.linalg.norm(AM))

        tall_predict = np.minimum(np.minimum(tv,tf),tc)
        tstep = tall_predict/10.0
        dt = 0.0
        if DEBUG:
                print tv,tf,tall_predict,tstep
        #sys.exit(0)


        boundary_distance = 0

        ### slide along dynamical path until ds-ball is hit with some
        ### tolerance

        [Acontrol,bcontrol] = GetControlConstraintMatricesAdjust(p,F,epsilon=boundary_distance)
        Ndim = p.shape[0]
        dbest = 1e5
        dtbest = 0

        ####################

        qcontrol = None
        #print "speed=",speed
        #dt = ds/(np.linalg.norm(dp)*speed)
        #while qcontrol is None:
        #        [qcontrol, qdcontrol, dcontrol_to_pnext] = GetNearestControlPoint(p, dp, pnext, F, dt, speed, Acontrol, bcontrol)
        #        print dt,dcontrol_to_pnext
        #        dt *= 0.5
        #dtbest = 2*dt
        #[qcontrol, qdcontrol, dcontrol_to_pnext] = GetNearestControlPoint(p, dp, pnext, F, dtbest, speed, Acontrol, bcontrol)
        #sys.exit(0)
        ####################

        #PrintNumpy('pnext2', pnext)
        #PrintNumpy('p', p)
        #PrintNumpy('dp', dp)
        #PrintNumpy('force', F)
        #print "ds=",ds
        #print "speed=",speed

        ictr = 0
        
        ICTR_STOP = 20
        while True:
                #### Solve QCQP -> get new distance
                [qcontrol, qdcontrol, dcontrol_to_pnext] = GetNearestControlPoint(p, dp, pnext, F, dt, speed, Acontrol, bcontrol)

                if qcontrol is not None:
                        ddd = np.linalg.norm(qcontrol-p)

                        if ddd > ds-ds/10:
                                dt -= tstep
                                if DEBUG:
                                        print "dtend:",dt,dbest
                                break
                        else:
                                dbest = dcontrol_to_pnext
                                dtbest = dt
                                if DEBUG:
                                        print "dnext",ddd,"dt:",dt,"d",dcontrol_to_pnext,"dbest",dbest,"dtbest:",dtbest
                                dt += tstep
                else:
                        dt -= tstep
                        tstep /=2
                        
                        
                #dt2 = 0.5*dt*dt
                #p + speed*dp*dt + dt2 * F + dt2 * control


                #if qcontrol is None:
                #        dt -= tstep
                #        tstep /= 2
                #else:
                #        if dcontrol_to_pnext <= dbest:
                #                dbest = dcontrol_to_pnext
                #                dtbest = dt
                #        else:
                #                dt -= tstep
                #                print "dtend:",dt,dbest
                #                break

                #dt += tstep


                #if ictr > ICTR_STOP:
                #        if qcontrol is None:
                #                print "ICTR:",ICTR_STOP
                #                ICTR_STOP += ICTR_STOP
                #        else:
                #                dt -= tstep
                #                break
                #ictr+=1

        [qcontrol, qdcontrol, dcontrol_to_pnext] = GetNearestControlPoint(p, dp, pnext, F, dtbest, speed, Acontrol, bcontrol)
        #ictr = 0

        if dtbest < 1e-100:
                print "dtbest",dtbest
                PrintNumpy('p',p)
                PrintNumpy('dp',dp)
                PrintNumpy('F',F)
                print "ds=",ds
                print "s=",speed
                #sys.exit(0)

        if dtbest > 5:
                print "dtbest",dtbest
                PrintNumpy('p',p)
                PrintNumpy('dp',dp)
                PrintNumpy('F',F)
                print "ds=",ds
                print "s=",speed
                #VisualizeReachableSet3D(p,dp,dp,speed,ds,F)
                sys.exit(0)

        #dold = dcontrol_to_pnext 
        #while qcontrol is not None and np.linalg.norm(qcontrol-p) <ds:
        #        [qcontrol, qdcontrol, dcontrol_to_pnext] = GetNearestControlPoint(p, dp, pnext, F, dtbest, speed, Acontrol, bcontrol)
        #        dtbest+=tstep
        #        ictr+=1
        #        if np.linalg.norm(dold-dcontrol_to_pnext)<1e-40:
        #                break

        #if ictr>0:
        #        dtbest-=tstep
        #        [qcontrol, qdcontrol, dcontrol_to_pnext] = GetNearestControlPoint(p, dp, pnext, F, dtbest, speed, Acontrol, bcontrol)

        #print dtbest,qcontrol,p,dp
        #plt.plot(qcontrol[0],qcontrol[1],'or',markersize=10)
        #plt.show()
        return [qcontrol,qdcontrol,dtbest]

def VisualizeReachableSet3D(p, dp, dwori, speed, ds, F):
        from reachable_set3d import ReachableSet3D
        Ndim = p.shape[0]
        R = ControlPerWaypoint(p, Ndim, 1)[:,:,0]
        reach = ReachableSet3D( ds, p, speed, dp, F, R, amin, amax)
        reach.Plot()
        reach.PlotShow()

