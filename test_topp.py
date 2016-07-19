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
#import statsmodels.api as sm
np.set_printoptions(precision=2)

if __name__ == "__main__":

        env = EnvironmentBloodStream()

        W = np.array(((0,0),(1,0),(2,0.5),(3,1.0))).T
        #dW = np.array(((1,0),(1,0.5),(1,0.5))).T
        #F = np.array(((0,0),(0,2),(0,1))).T

        print W.shape

        #traj = Trajectory(W)
        #topp = TOPPInterface(traj, None, None, F, R, amin, amax, W, dW)
        from TOPP import Utilities
        traj0 = Utilities.InterpolateViapoints(W)
        traj0.Plot(0.1)
        plt.show()

        #topp = TOPPInterface(env, traj)
        #print topp.durationVector
        #topp.getCriticalPoint()
