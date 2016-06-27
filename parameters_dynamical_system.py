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

AM = 1
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

#def GetNearestControlPoint(p, dp, speed, ds, F):
#
#        [qnext,dqnext,dt] = ForwardSimulate(p, dp, speed, ds, F)
#        return qnext
#
#        dt2 = dt*dt*0.5
#        Ndim = qnext.shape[0]
#
#        ## solve LP
#        ### min x^T*c
#        ### s.t. A*x + b <= 0
#
#        [A,bcontrol] = GetControlConstraintMatricesAdjust(p,F,epsilon=0.1)
#
#        A = -A
#        A = (2*A)/(dt*dt)
#        b = bcontrol + np.dot(A,(-p - dt*speed*dp))
#        A = matrix(A)
#        b = matrix(b)
#
#        nF = F/np.linalg.norm(F)
#        ndP = dp/np.linalg.norm(dp)
#        gamma = 0.1
#        c = matrix(nF - gamma*ndP)
#
#        try:
#                solvers.options['show_progress'] = False
#                solvers.options['maxiters'] = 100 ##default: 100
#                solvers.options['abstol'] = 1e-03 ## below 1e-11 some weird behavior
#                solvers.options['reltol'] = 1e-03 ## below 1e-11 some weird behavior
#                solvers.options['feastol'] = 1e-03 ## below 1e-11 some weird behavior
#                sol=solvers.lp(c,A,-b)
#                x = np.array(sol['x']).flatten()
#
#        except ValueError:
#                return qnext
#
#        except Exception as e:
#                print e
#                print "A:",A
#                print "b:",b
#                print "c:",c
#                print "qnext:",qnext
#                PrintNumpy('p', p)
#                PrintNumpy('dp', dp)
#                PrintNumpy('force', F)
#                print "ds=",ds
#                print "speed=",speed
#                #VisualizeReachableSet3D(p, dp, dp, speed, ds, F)
#                sys.exit(0)
#
#        return x

#def ForwardSimulate_deprecated(p, dp, smax, ds, F):
#        ### should best follow path!
#        if np.linalg.norm(F)>1e-3:
#                qnext = copy.copy(p)
#                tstep = 1e-3
#                dt = 0.0
#                dnew = 0.0
#                #ictr=0
#
#                tolerance = 1e-6
#
#                ### slide along dynamical path until ds-ball is hit with some
#                ### tolerance
#
#
#                while np.linalg.norm(dnew-ds) > tolerance:
#                        if dnew < ds:
#                                dt += tstep
#                        else:
#                                ## overshoot, step back
#                                tstep = tstep/2.0
#                                dt -= tstep
#
#                        dt2 = dt*dt/2
#                        qnext = p + dt*smax*dp + dt2*F
#                        dnew = np.linalg.norm(p-qnext)
#                dqnext = smax*dp + dt*F
#                #print qnext
#                return [qnext,dqnext,dt]
#        return None
def GetNearestControlPoint(p, dp, pnext, F, dt, speed, Acontrol, bcontrol):
        ds_next = np.linalg.norm(p-pnext)
        Ndim = p.shape[0]
        if dt < 1e-100:
                #print "dt<1e-100",dt
                return [p,ds_next]

        dt2 = dt*dt/2

        A = np.zeros((Acontrol.shape))
        b = np.zeros((bcontrol.shape))

        A = (-2*Acontrol)/(dt*dt)
        b = bcontrol + np.dot(A,(-p - dt*speed*dp))

        Id = np.eye(Ndim)
        x = Variable(Ndim)

        objective = Minimize( norm(x - pnext))
        constraints = [ A*x <= -b,
                        quad_form(x-p,Id) <= ds_next]

        prob = Problem(objective, constraints)
        dnew = np.abs(prob.solve(solver=ECOS, max_iters=100, feastol_inacc=1e-5))

        if dnew < inf:
                qcontrol =np.array(x.value).flatten()
        else:
                qcontrol = None

        return [qcontrol, dnew]

def ForwardSimulate(p, dp, speed, ds, F):
        ### should best follow path!
        pnext = p + ds*dp/np.linalg.norm(dp)

        ### estimate tstep
        dds = np.linalg.norm(speed*dp)
        tstep = ds/(10*dds)
        #tstep = 1e-3

        dt = 0.0

        ## distance from RS boundary
        boundary_distance = 0
        ### slide along dynamical path until ds-ball is hit with some
        ### tolerance

        [Acontrol,bcontrol] = GetControlConstraintMatricesAdjust(p,F,epsilon=boundary_distance)
        Ndim = p.shape[0]
        dbest = 1e5
        dtbest = 0

        #PrintNumpy('pnext2', pnext)
        #PrintNumpy('p', p)
        #PrintNumpy('dp', dp)
        #PrintNumpy('force', F)
        #print "ds=",ds
        #print "speed=",speed

        ictr = 0
        
        #tstep = 1e-3
        #plt.figure(100)
        #plt.plot(p[0],p[1],'ok',markersize=20)
        #plt.plot(pnext[0],pnext[1],'og',markersize=20)
        while True:
                #### Solve QCQP -> get new distance
                [qcontrol, dcontrol_to_pnext] = GetNearestControlPoint(p, dp, pnext, F, dt, speed, Acontrol, bcontrol)
                #print ictr,"dt:",dt,"d",dcontrol_to_pnext

                if qcontrol is None:
                        dt -= tstep
                        tstep /= 2.0
                else:
                        #plt.plot(qcontrol[0],qcontrol[1],'or')
                        if dcontrol_to_pnext < dbest:
                                dbest = dcontrol_to_pnext
                                dtbest = dt
                        else:
                                break

                dt += tstep

                if ictr > 20:
                        break
                ictr+=1
                #dt += tstep

        #print "DT:",dtbest,dbest
        [qcontrol, dcontrol_to_pnext] = GetNearestControlPoint(p, dp, pnext, F, dtbest, speed, Acontrol, bcontrol)
        #plt.plot(qcontrol[0],qcontrol[1],'or',markersize=10)
        #plt.show()
        return [qcontrol,None,dtbest]

def VisualizeReachableSet3D(p, dp, dwori, speed, ds, F):
        from reachable_set3d import ReachableSet3D
        Ndim = p.shape[0]
        R = ControlPerWaypoint(p, Ndim, 1)[:,:,0]
        reach = ReachableSet3D( ds, p, speed, dp, F, R, amin, amax)
        reach.Plot()
        reach.PlotShow()

