import numpy as np
from numpy import sqrt,sin,cos,pi

FILENAME = 'virus2'
DYNAMICAL_SYTEM_NAME = 'non-holonomic differential drive (deform)'
system_name = DYNAMICAL_SYTEM_NAME

AM = 2
### car/sailboat
#amin = np.array((-AM,-AM,-AM,0))
#amax = np.array((AM,AM,AM,0))
amin = np.array((-AM,-AM,-AM))
amax = np.array((AM,AM,AM))

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
        R = ControlPerWaypoint(p, Ndim, 1)
        return GetControlConstraintMatricesFromControl(R,F)

def GetControlConstraintMatricesFromControl(R, F):
        Ndim = R.shape[0]
        Rmax = np.maximum(np.dot(R,amin),np.dot(R,amax))
        Rmin = np.minimum(np.dot(R,amin),np.dot(R,amax))
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
