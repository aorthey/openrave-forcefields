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

from deformation_naive import *
from deformation_potentials import *
from deformation_stretchpull import *
from deformation_reachableset import *
import numpy as np
from motion_planner_geometrical import MotionPlannerGeometrical
from motion_planner_kinodynamic import MotionPlannerKinodynamic
from topp_interface import *
from topp_interface2 import *
from trajectory import Trajectory
from TOPP import Trajectory as TOPPTrajectory
#import statsmodels.api as sm
np.set_printoptions(precision=2)

if __name__ == "__main__":

        env = EnvironmentBloodStream()

        W = np.array(((0,0,0,0),(1,0,0,0),(2,0.1,0,0),(3,0.3,0,0))).T
        #W = np.array(((1,0),(1,0),(1,0.5),(1,1.0))).T

        #dW = np.array(((1,0),(1,0.5),(1,0.5))).T
        #F = np.array(((0,0),(0,2),(0,1))).T

        from TOPP import Utilities
        traj0 = Utilities.InterpolateViapoints(W)

        Ndim = W.shape[0]
        Adim = 3

        #traj0.Plotd(0.01)
        #traj0.Plotdd(0.01)
        ### DEBUG
        ####dt = 0.01
        ####tvect = arange(0, traj0.duration + dt, dt)
        ####qvect = array([traj0.Eval(t) for t in tvect]).T
        ####plot(qvect[0,:], qvect[1,:], '-or', linewidth=2)
        ####plot(tvect, qvect.T, linewidth=2)
        ####plt.show()

        #######################################################################
        #######################################################################
        discrtimestep = 0.01
        ndiscrsteps = int((traj0.duration + 1e-10) / discrtimestep) + 1
        q = np.zeros((Ndim,ndiscrsteps))
        qs = np.zeros((Ndim,ndiscrsteps))
        qss = np.zeros((Ndim,ndiscrsteps))
        a = np.zeros((ndiscrsteps, 2*Adim))
        b = np.zeros((ndiscrsteps, 2*Adim))
        c = np.zeros((ndiscrsteps, 2*Adim))
        F = np.zeros((ndiscrsteps, Ndim))


        for i in range(ndiscrsteps):
                t = i*discrtimestep

                q[:,i] = traj0.Eval(t)
                qs[:,i] = traj0.Evald(t)
                qss[:,i] = traj0.Evaldd(t)

                if q[0,i] > 2:
                        F[i,:] = np.array((0,1.1,0,0))

                Ri = params.GetControlMatrixAtWaypoint(q[:,i])

                [G,h] = params.GetControlConstraintMatricesFromControl(Ri, F[i,:])
                a[i,:] = np.dot(G,qs[:,i]).flatten()
                b[i,:] = np.dot(G,qss[:,i]).flatten()
                c[i,:] = h

        ### DEBUG PLOT
        dt = discrtimestep
        #tvect = arange(0, traj0.duration + dt, dt)
        plot(q[0,:], q[1,:], '-or', linewidth=2)

        for i in range(ndiscrsteps):
                plot([q[0,i],q[0,i]+F[i,0]],[q[1,i],q[1,i]+F[i,1]], '-ob', linewidth=2)

        plt.show()

        ######################################################
        t0 = time.time()
        trajectorystring = str(traj0)
        vmax = 1e5*np.ones(Ndim)
        topp_inst = TOPP.QuadraticConstraints(traj0, discrtimestep, vmax, list(a), list(b), list(c))

        x = topp_inst.solver
        t1 = time.time()
        ret = x.RunComputeProfiles(0,0)
        print "TOPP:",ret
        print t1-t0
        if ret==1:
                x.ReparameterizeTrajectory()
                t2 = time.time()
                x.ReparameterizeTrajectory()
                x.WriteResultTrajectory()
                traj1 = TOPPTrajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)

        #topp = TOPPInterface(env, traj)
        #print topp.durationVector
        #topp.getCriticalPoint()
