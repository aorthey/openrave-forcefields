import numpy as np
import sys
import copy
from numpy import sqrt,sin,cos,pi
from qcqp import *
import scipy
import cvxopt
from cvxpy import *
from cvxopt import matrix, solvers
from util import PrintNumpy,inf

FILENAME = 'virus2'
DYNAMICAL_SYTEM_NAME = 'non-holonomic differential drive (deform)'
system_name = DYNAMICAL_SYTEM_NAME

AM = 2
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
        #Ndim = R.shape[0]
        #Rmin = np.minimum(np.dot(R,amin),np.dot(R,amax))
        #Rmax = np.maximum(np.dot(R,amin),np.dot(R,amax))
        #H1 = F - Rmax
        #H2 = -F + Rmin
        #for j in range(Ndim):
        #        if H2[j] > -H1[j]:
        #                print H2[j],"<= q[",j,"]<=",-H1[j]
        #                sys.exit(1)
        #b = np.hstack((H1,H2)).flatten()
        #I = np.identity(Ndim)
        #A = np.vstack((I,-I))
        #return [A,b]
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

def GetNearestControlPoint(p, dp, speed, ds, F):

        [qnext,dqnext,dt] = ForwardSimulate(p, dp, speed, ds, F)
        return qnext

        dt2 = dt*dt*0.5
        Ndim = qnext.shape[0]

        ## solve LP
        ### min x^T*c
        ### s.t. A*x + b <= 0

        [A,bcontrol] = GetControlConstraintMatricesAdjust(p,F,epsilon=0.1)

        A = -A
        A = (2*A)/(dt*dt)
        b = bcontrol + np.dot(A,(-p - dt*speed*dp))
        A = matrix(A)
        b = matrix(b)

        nF = F/np.linalg.norm(F)
        ndP = dp/np.linalg.norm(dp)
        gamma = 0.1
        c = matrix(nF - gamma*ndP)

        try:
                solvers.options['show_progress'] = False
                solvers.options['maxiters'] = 100 ##default: 100
                solvers.options['abstol'] = 1e-03 ## below 1e-11 some weird behavior
                solvers.options['reltol'] = 1e-03 ## below 1e-11 some weird behavior
                solvers.options['feastol'] = 1e-03 ## below 1e-11 some weird behavior
                sol=solvers.lp(c,A,-b)
                x = np.array(sol['x']).flatten()

        except ValueError:
                return qnext

        except Exception as e:
                print e
                print "A:",A
                print "b:",b
                print "c:",c
                print "qnext:",qnext
                PrintNumpy('p', p)
                PrintNumpy('dp', dp)
                PrintNumpy('force', F)
                print "ds=",ds
                print "speed=",speed
                #VisualizeReachableSet3D(p, dp, dp, speed, ds, F)
                sys.exit(0)

        return x

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

import pylab as plt
def ForwardSimulate(p, dp, speed, ds, F):
        ### should best follow path!
        qnext = copy.copy(p)
        pnext = p + ds*dp/np.linalg.norm(dp)
        tstep = 1e-3
        dt = 0.0 + tstep

        tolerance = 1e-5

        ## distance from RS boundary
        boundary_distance = 0.5

        ### slide along dynamical path until ds-ball is hit with some
        ### tolerance

        [Acontrol,bcontrol] = GetControlConstraintMatricesAdjust(p,F,epsilon=boundary_distance)
        Ndim = qnext.shape[0]
        #for i in range(0,Ndim):
                #print bcontrol[i+Ndim],"<=","x[",i,"] <=",-bcontrol[i] 

        #while np.linalg.norm(dnew-ds) > tolerance:
        dold = 1e5

        #print "START:",p
        while True:
                #SolveQCQP(dt, Acontrol, bcontrol, p, dp, speed, ds, F)

                #### Solve QCQP -> get new distance
                dt2 = dt*dt/2
                qnext = p + dt*speed*dp + dt2*F
        
                A = np.zeros((Acontrol.shape))
                b = np.zeros((bcontrol.shape))

                A = (-2*Acontrol)/(dt*dt)
                b = bcontrol + np.dot(A,(-p - dt*speed*dp))

                #PrintNumpy('pnext',pnext)
                #PrintNumpy('A',A)
                #PrintNumpy('b',b)
                #print "ds=",ds

                Id = np.eye(Ndim)
                x = Variable(Ndim)

                objective = Minimize( norm(x - pnext)) #-c.T*x )
                constraints = [ 0.5*quad_form(x,Id) <= ds*ds,
                                A*x + b <= 0 ]

                prob = Problem(objective, constraints)
                dnew = prob.solve(solver=SCS)

                if dnew < inf:
                        #print ds,dt,dnew,dold
                        qcontrol =np.array(x.value).flatten()
                        dt += tstep
                        if dnew >= dold:
                                break
                        else:
                                #print "dold:",dold,"->",dnew
                                dold = copy.copy(dnew)
                else:
                        ## overshoot, step back
                        tstep = tstep/2.0
                        dt -= tstep


                        #print "No Solution QCQP"
                        #print "qnext:",ds
                        #print "pnext:",pnext
                        #print "A:",A
                        #print "b:",b
                        #qcontrol = qnext 
                        #sys.exit(0)

                #dnew = np.linalg.norm(p-qcontrol)
                #try:
                #        #res = scipy.optimize.linprog(c,A,-b)
                #        #if res.success:
                #        #        qcontrol = np.array(res.x).flatten()
                #        #else:
                #        #        raise ValueError
                #        #solvers.options['show_progress'] = True
                #        #solvers.options['maxiters'] = 100 ##default: 100
                #        solvers.options['abstol'] = 1e-20 ## below 1e-11 some weird behavior
                #        solvers.options['reltol'] = 1e-20 ## below 1e-11 some weird behavior
                #        #solvers.options['feastol'] = 1e-03 ## below 1e-11 some weird behavior
                #        A=cvxopt.matrix(A)
                #        b=cvxopt.matrix(b)
                #        #c=cvxopt.matrix(c)
                #        #sol=solvers.lp(c,A,-b)
                #        #qcontrol = np.array(sol['x']).flatten()

                #        #Id = cvxopt.matrix(np.eye(Ndim))
                #        #ZZ = cvxopt.matrix(np.zeros((Ndim,Ndim)))


                #        CQ = np.eye(Ndim)
                #        CQ = cvxopt.matrix(CQ)
                #        Cp = -2*pnext
                #        Cp = cvxopt.matrix(Cp)
                #        solvers.options['feastol'] = 1e-100 ## below 1e-11 some weird behavior
                #        solvers.options['show_progress'] = False
                #        sol=solvers.qp(CQ,Cp, A,-b,None,None)
                #        qcontrol = np.array(sol['x']).flatten()

                #        print qcontrol,pnext


                #except ValueError:
                #        print "ValueError"
                #        print "qnext:",qnext
                #        print "pnext:",pnext
                #        print "A:",A
                #        print "b:",b
                #        print "c:",c
                #        qcontrol = qnext 
                #        sys.exit(0)



        #dqnext = speed*dp + dt*F
        return [qcontrol,None,dt]

def VisualizeReachableSet3D(p, dp, dwori, speed, ds, F):
        from reachable_set3d import ReachableSet3D
        Ndim = p.shape[0]
        R = ControlPerWaypoint(p, Ndim, 1)[:,:,0]
        reach = ReachableSet3D( ds, p, speed, dp, F, R, amin, amax)
        reach.Plot()
        reach.PlotShow()

