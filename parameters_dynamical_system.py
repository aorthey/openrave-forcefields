import numpy as np
from numpy import sqrt,sin,cos,pi

FILENAME = 'virus2'
DYNAMICAL_SYTEM_NAME = 'non-holonomic differential drive (deform)'
system_name = DYNAMICAL_SYTEM_NAME

AM = 1
### car/sailboat
amin = np.array((-AM,-AM,-0.5*AM,0))
amax = np.array((AM,AM,0.5*AM,0))

def ControlPerWaypoint(W, Ndim, Nwaypoints):
        assert(Ndim==4)
        Kdim = 4
        R = np.zeros((Ndim,Kdim,Nwaypoints))
        for i in range(0,Nwaypoints):
                if Nwaypoints>1:
                        t = W[3,i]
                else:
                        t = W[3]

                R[0,:,i] = np.array((cos(t),-sin(t),0.0,0.0))
                R[1,:,i] = np.array((sin(t),cos(t),0.0,0.0))
                R[2,:,i] = np.array((0.0,0.0,0.0,1.0))
                R[3,:,i] = np.array((0.0,0.0,1.0,0.0))

        return R

