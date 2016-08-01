#!/usr/bin/env python
import time
import scipy
import sys
import numpy as np
import openravepy
from openravepy import *
from math import *
from environment_force_the_stream import *
from environment_force_the_counterstream import *
from environment_force_the_ray import *
from environment_force_blood_stream import *
from environment_periodic_force_the_hideout import *
from environment_periodic_force_triple_stream import *
from environment_periodic_force_crossroad_stream import *
import numpy as np
#from trajectory import Trajectory
from TOPP import Trajectory as TOPPTrajectory
from TOPP import Utilities
import TOPP
from pylab import *
from util import *
from scipy import interpolate
import parameters_dynamical_system as params

#import statsmodels.api as sm
np.set_printoptions(precision=2)


def InterpolateViapointsLinear(path):
    nviapoints = len(path[0,:])
    ### assumption: distance between waypoints is always the same and can be scaled to 1.0/nviapoints

    tv = np.zeros(nviapoints)
    for i in range(nviapoints-1):
            tv[i+1] = tv[i]+np.linalg.norm(path[:,i]-path[:,i+1])

    tv = linspace(0,1,nviapoints)
    tcklist = []
    for idof in range(0,path.shape[0]):
        tcklist.append(interpolate.splrep(tv,path[idof,:],s=0,k=1))
    t = tcklist[0][0]
    chunkslist = []
    for i in range(len(t)-1):
        polylist = []
        if abs(t[i+1]-t[i])>1e-5:
            for tck in tcklist:
                a = 0
                b = 0
                c = interpolate.splev(t[i],tck,der=1)
                d = interpolate.splev(t[i],tck,der=0)
                polylist.append(TOPPTrajectory.Polynomial([d,c,b,a]))
            chunkslist.append(TOPPTrajectory.Chunk(t[i+1]-t[i],polylist))
    return TOPPTrajectory.PiecewisePolynomialTrajectory(chunkslist)

def InterpolateViapointsLinearFixedDerivative(path):
    nviapoints = len(path[0,:])
    tv = linspace(0,1,nviapoints)
    tcklist = []
    for idof in range(0,path.shape[0]):
        tcklist.append(interpolate.splrep(tv,path[idof,:],s=0,k=1))

    t = tcklist[0][0]
    chunkslist = []
    for i in range(len(t)-1):
        polylist = []
        if abs(t[i+1]-t[i])>1e-5:
            for tck in tcklist:
                a = 0
                b = 0
                c = interpolate.splev(t[i],tck,der=1)
                d = interpolate.splev(t[i],tck,der=0)
                polylist.append(TOPPTrajectory.Polynomial([d,c,b,a]))
            chunkslist.append(TOPPTrajectory.Chunk(t[i+1]-t[i],polylist))
    return TOPPTrajectory.PiecewisePolynomialTrajectory(chunkslist)

def VisualizeTrajectory(traj0, W, Fx, Fy, offset, PostMakeup=True):

        Ndim = W.shape[0]
        #Ndim = 2
        Adim = 3
        #######################################################################
        #######################################################################
        #discrtimestep = 1.0/(W.shape[1]-1)#0.05
        discrtimestep = 0.5*1e-3
        #discrtimestep = 1e-2
        ndiscrsteps = int((traj0.duration + 1e-10) / discrtimestep) + 1
        q = np.zeros((Ndim,ndiscrsteps))
        qs = np.zeros((Ndim,ndiscrsteps))
        qss = np.zeros((Ndim,ndiscrsteps))
        a = np.zeros((ndiscrsteps, 2*Adim))
        b = np.zeros((ndiscrsteps, 2*Adim))
        c = np.zeros((ndiscrsteps, 2*Adim))
        F = np.zeros((ndiscrsteps, Ndim))

        tvect = arange(0, traj0.duration + discrtimestep, discrtimestep)
        for i in range(ndiscrsteps):
                t = i*discrtimestep

                q[:,i] = traj0.Eval(t)
                qs[:,i] = traj0.Evald(t)
                qss[:,i] = traj0.Evaldd(t)

                if q[0,i] >= 2:
                        F[i,:] = np.array((Fx,Fy,0,0))

                Ri = params.GetControlMatrixAtWaypoint(q[:,i])

                [G,h] = params.GetControlConstraintMatricesFromControl(Ri, F[i,:])
                a[i,:] = np.dot(G,qs[:,i]).flatten()
                b[i,:] = np.dot(G,qss[:,i]).flatten()
                c[i,:] = h

        #print "discrtimestep=",discrtimestep
        #PrintNumpy("a",a[:,[0,1,3,4]])
        #PrintNumpy("b",b[:,[0,1,3,4]])
        #PrintNumpy("c",c[:,[0,1,3,4]])
        #print a.shape
        #print b.shape
        #print c.shape
        ######################################################
        t0 = time.time()
        trajectorystring = str(traj0)
        vmax = 1e5*np.ones(Ndim)
        topp_inst = TOPP.QuadraticConstraints(traj0, discrtimestep, vmax, list(a), list(b), list(c))

        x = topp_inst.solver
        t1 = time.time()
        #ret = x.RunComputeProfiles(0,0)
        ret = x.RunVIP(0,0)
        print "TOPP:",ret
        #print t1-t0
        if ret==1:
                x.ReparameterizeTrajectory()
                t2 = time.time()
                x.WriteResultTrajectory()
                msstyle = 'g'
                try:
                        #print x.restrajectorystring
                        #traj1 = TOPPTrajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
                        semin = x.sdendmin
                        semax = x.sdendmax
                        print "speed profile [",semin,",",semax,"]"
                        #plot(qddvect[:,0], qddvect[:,1]-offset, 'r')
                        #if PostMakeup:
                                #tvect = arange(0, traj1.duration + discrtimestep, discrtimestep)
                                #qddvect = np.array([traj1.Evaldd(t) for t in tvect])
                                #for i in range(ndiscrsteps):
                                        #plot([q[0,i],q[0,i]+qddvect[i,0]],[q[1,i]-offset,q[1,i]+qddvect[i,1]-offset], '-om', linewidth=2)

                        if semin > 0:
                                msstyle = 'm'

                except Exception as e:
                        print "Exception:",e
                        msstyle = 'k'
                        sys.exit(0)
                        pass
        else:
                msstyle = 'r'

        if PostMakeup:
                for i in range(ndiscrsteps):
                        plot([q[0,i],q[0,i]+F[i,0]],[q[1,i]-offset,q[1,i]+F[i,1]-offset], '-ob', linewidth=2)

                plot(q[0,:], q[1,:]-offset, '|'+msstyle, markersize=40,markeredgewidth=3)
        Npts = 1000
        tvect = np.linspace(0,traj0.duration, Npts)
        qvect = np.array([traj0.Eval(t) for t in tvect])

        plot(qvect[:,0], qvect[:,1]-offset, '-'+msstyle)
        plot(W[0,:], W[1,:]-offset, 'o'+msstyle, markersize=15)

def PlotTrajXY(x,y, PostMakeup = False):
        Fx = 0.0
        #Fy = 2.144
        Fy = 1.1
        #W = np.array(((0,0,0,0),(1,0,0,0),(2,0,0,0),(x,y,0,0))).T
        W = np.array(((0,0,0,0),(0.2,0,0,0),(0.4,0,0,0),(0.6,0,0,0),(0.8,0,0,0),(1,0,0,0),(1.2,0,0,0),(1.4,0,0,0),(1.6,0,0,0),(1.8,0,0,0),(2,0,0,0),(x,y,0,0))).T
        #W = np.array(((1.6,0,0,0),(1.8,0,0,0),(2,0,0,0),(x,y,0,0))).T
        traj0 = Utilities.InterpolateViapoints(W)
        #traj0 = InterpolateViapointsLinear(W)
        #print "#######################"
        #print "trajectorystring=",traj0
        #print "#######################"
        VisualizeTrajectory(traj0, W, Fx, Fy, PostMakeup = PostMakeup, offset=1.0)
        return traj0

if __name__ == "__main__":

        env = EnvironmentBloodStream()


        ### Visualize singularities imposed by interpolation method
        ###W = np.array(((0,0,0,0),(1,0,0,0),(2,0,0,0),(3.5,0.5,0,0))).T
        ###traj0 = InterpolateViapoints(W)
        ####VisualizeTrajectory(traj0, W, Fx, Fy, PostMakeup = False, offset=0.0)
        ###W = np.array(((1.05,-0.1,0,0),(1,0,0,0),(1.02,0,0,0),(0.95,0.1,0,0))).T
        ###traj0 = InterpolateViapoints(W)
        ###VisualizeTrajectory(traj0, W, Fx, Fy, PostMakeup = False, offset=2.0)
        ###plt.show()
        ### sys.exit(0)

        #### visualize all 
        r=0.2
        xc = 2.0
        yc = 0.0
        epsilon = pi/32
        for r in np.linspace(0.2,2.0,10):
                for theta in -np.linspace(-pi/2+epsilon,pi/2-epsilon,10):
                        x = r*cos(theta)+xc
                        y = r*sin(theta)+yc
                        PlotTrajXY(x,y,PostMakeup=False)

        #x = 2.13160997913
        #y = -0.150594865098
        #traj0 = PlotTrajXY(x,y, PostMakeup=False)
        #PlotTrajXY(2.0,0.0, PostMakeup=True)
        #PlotTrajXY(4.0,-4.0)
        ###PlotTrajXY(3.5,-2.8)


        #### visualize interpolation breakdown
        #x = -3.0
        #y = 4.0
        #W = np.array(((0,0,0,0),(1,0,0,0),(2,0,0,pi/16),(x,y,0,pi/4))).T
        #traj0 = Utilities.InterpolateViapoints(W)
        #VisualizeTrajectory(traj0, W, Fx, Fy, PostMakeup = False, offset=1.0)

        #plt.axis('equal')
        plt.show()
