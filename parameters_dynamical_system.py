import numpy as np
import copy
from numpy import sqrt,sin,cos,pi
from cvxopt import matrix, solvers

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
        Ndim = p.shape[0]
        R = ControlPerWaypoint(p, Ndim, 1)[:,:,0]
        return GetControlConstraintMatricesFromControl(R,F)

def GetControlConstraintMatricesFromControl(R, F):
        Ndim = R.shape[0]
        Rmin = np.minimum(np.dot(R,amin),np.dot(R,amax))
        Rmax = np.maximum(np.dot(R,amin),np.dot(R,amax))
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
        dt2 = dt*dt*0.5
        Ndim = qnext.shape[0]

        ## solve LP
        ### min x^T*c
        ### s.t. A*x + b <= 0

        [A,bcontrol] = GetControlConstraintMatrices(p,F)

        #for i in range(0,Ndim):
                #print bcontrol[i+Ndim],"<=","x[",i,"] <=",-bcontrol[i]

        #M = 1000
        #qq = np.zeros((Ndim,M))
        #qq = []
        #A = matrix(A)
        #b = matrix(bcontrol)
        #for i in range(0,M):
        #        qx = np.random.uniform(-5,5)
        #        qy = np.random.uniform(-5,5)
        #        qz = np.random.uniform(-5,5)
        #        qt = np.random.uniform(-6,5)
        #        q = np.array((qx,qy,qz,qt))
        #        #qc = np.dot(A,q) + bcontrol 
        #        q = matrix(q)
        #        sol=solvers.lp(q,-A,-b)
        #        qc = np.array(sol['x']).flatten()
        #        qqc = p + dt*speed*dp + dt2*qc
        #        qq.append(qqc)

        ##qq = p + dt*speed*dp + dt2*F
        ##qq = np.dot(A,qnext) + np.dot(A, (-p-dt*speed*dp))+bcontrol
        #qq = np.array(qq)
        #ax.scatter(qq[:,0],qq[:,1],qq[:,3], 'or',
        #                s=50)

        A = -A
        A = (2*A)/(dt*dt)
        b = bcontrol + np.dot(A,(-p - dt*speed*dp))
        A = matrix(A)
        b = matrix(b)

        #for i in range(0,Ndim):
                #print b[i+Ndim],"<=","x[",i,"] <=",-b[i]

        nF = F/np.linalg.norm(F)
        ndP = dp/np.linalg.norm(dp)
        gamma = 0.1
        c = matrix(nF - gamma*ndP)

        solvers.options['show_progress'] = False
        try:
                sol=solvers.lp(c,A,-b)
                x = np.array(sol['x']).flatten()
        except Exception as e:
                print e
                print "A:",A
                print "b:",b
                print "c:",c
                print "p:",p
                print "dp:",dp
                print "ds:",ds
                print "speed:",speed
                print "F:",F
                #VisualizeReachableSet3D(p, dp, dp, speed, ds, F)
                sys.exit(0)

        return x

def ForwardSimulate(p, dp, smax, ds, F):
        if np.linalg.norm(F)>1e-3:
                qnext = copy.copy(p)
                tstep = 1e-3
                dt = 0.0

                dnew = 0.0
                ictr=0
                while dnew < ds:
                        dt += tstep
                        dt2 = dt*dt/2
                        qnext = p + dt*smax*dp + dt2*F
                        dnew = np.linalg.norm(p-qnext)
                        ictr+=1

                dqnext = smax*dp + dt*F
                return [qnext,dqnext,dt]
        return None

def VisualizeReachableSet3D(p, dp, dwori, speed, ds, F):
        from reachable_set3d import ReachableSet3D
        Ndim = p.shape[0]
        R = ControlPerWaypoint(p, Ndim, 1)[:,:,0]
        reach = ReachableSet3D( ds, p, speed, dp, F, R, amin, amax)
        reach.Plot()
        reach.PlotShow()

