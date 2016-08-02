#!/usr/bin/env python
import time
import numpy as np
from TOPP import Trajectory as TOPPTrajectory
from TOPP import Utilities
import TOPP
from pylab import *
from scipy import interpolate

np.set_printoptions(precision=2)

def TOPPInterface(traj0, W):
        Ndim = 2
        Adim = 2
        discrtimestep = 1e-3
        ndiscrsteps = int((traj0.duration + 1e-10) / discrtimestep) + 1
        q = np.zeros((Ndim,ndiscrsteps))
        qs = np.zeros((Ndim,ndiscrsteps))
        qss = np.zeros((Ndim,ndiscrsteps))
        a = np.zeros((ndiscrsteps, 2*Adim))
        b = np.zeros((ndiscrsteps, 2*Adim))
        c = np.zeros((ndiscrsteps, 2*Adim))
        F = np.zeros((ndiscrsteps, Ndim))

        tvect = arange(0, traj0.duration + discrtimestep, discrtimestep)
        umax = np.array((1,1))
        umin = np.array((-1,-1))
        Fborder = 0.2
        for i in range(ndiscrsteps):
                t = i*discrtimestep

                q[:,i] = traj0.Eval(t)
                qs[:,i] = traj0.Evald(t)
                qss[:,i] = traj0.Evaldd(t)

                Fi = np.array((0,0))
                if q[0,i] < Fborder:
                        Fi = np.array((1.5,0))

                I = np.eye(Ndim)
                G = np.vstack((I,-I))
                h = np.hstack((-Fi-umax,Fi+umin))
                a[i,:] = np.dot(G,qs[:,i]).flatten()
                b[i,:] = np.dot(G,qss[:,i]).flatten()
                c[i,:] = h

        vmax = 1e5*np.ones(Ndim)
        topp_inst = TOPP.QuadraticConstraints(traj0, discrtimestep, vmax, list(a), list(b), list(c))
        x = topp_inst.solver
        #ret = x.RunComputeProfiles(0,0)
        ret = x.RunVIP(0,0)
        if ret==1:
                print "TOPP VIP solution:",x.sdendmin,x.sdendmax
        else:
                print "TOPP VIP failure"

def InterpolateViapointsCustom(path):
    nviapoints = len(path[0,:])

    tv = np.zeros(nviapoints)
    for i in range(nviapoints-1):
            tv[i+1] = tv[i]+np.linalg.norm(path[:,i]-path[:,i+1])
    tcklist = []
    for idof in range(0,path.shape[0]):
        tcklist.append(interpolate.splrep(tv,path[idof,:],s=0,k=3))
    t = tcklist[0][0]
    chunkslist = []
    for i in range(len(t)-1):
        polylist = []
        if abs(t[i+1]-t[i])>1e-5:
            for tck in tcklist:
                a = 1/6. * interpolate.splev(t[i],tck,der=3)
                b = 0.5 * interpolate.splev(t[i],tck,der=2)
                c = interpolate.splev(t[i],tck,der=1)
                d = interpolate.splev(t[i],tck,der=0)
                polylist.append(TOPPTrajectory.Polynomial([d,c,b,a]))
            chunkslist.append(TOPPTrajectory.Chunk(t[i+1]-t[i],polylist))
    return TOPPTrajectory.PiecewisePolynomialTrajectory(chunkslist)

def PlotTraj(W):
        traj0 = Utilities.InterpolateViapoints(W)

        print "Interpolate Viapoints Default",
        M = 5000
        q = np.array([traj0.Eval(t) for t in np.linspace(0,traj0.duration,M)]).T
        plt.plot(q[0,:],q[1,:],'-r',linewidth=2)
        plt.plot(W[0,:],W[1,:],'ok')
        TOPPInterface(traj0,W)

        print "Interpolate Viapoints Adjust",
        traj0 = InterpolateViapointsCustom(W)
        q = np.array([traj0.Eval(t) for t in np.linspace(0,traj0.duration,M)]).T
        plt.plot(q[0,:],q[1,:],'-g',linewidth=2)
        plt.plot(W[0,:],W[1,:],'ok')
        TOPPInterface(traj0,W)

        print "quiver it"
        X, Y = np.mgrid[-0.1:0.5:15j, -0.005:0.005:15j]
        V = X+Y
        U = X+Y
        U[X>0.19]=0
        U[X<0.19]=1.5
        V = 0.0*V
        speed = np.sqrt(U**2 + V**2)
        UN = U/speed
        VN = V/speed

        plt.quiver(X, Y, UN, VN, color='Teal', headlength=7)
        plt.show()

if __name__ == "__main__":

        W=np.array(((0,0),(0.05,0.001),(0.35,-0.001),(0.4,0))).T
        PlotTraj(W)
