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
from environment_force_radial import *
from environment_force_blood_stream2 import *
from environment_periodic_force_the_hideout import *
from environment_periodic_force_triple_stream import *
from environment_periodic_force_crossroad_stream import *

from deformation_naive import *
from deformation_potentials import *
from deformation_stretchpull import *
from deformation_reachableset import *
from trajectory import *
import numpy as np
from motion_planner_geometrical import MotionPlannerGeometrical
from motion_planner_kinodynamic import MotionPlannerKinodynamic
from dynamical_system import *
#import statsmodels.api as sm
np.set_printoptions(precision=2)

if __name__ == "__main__":

        #######################################################################
        env = EnvironmentRadial()
        #env = EnvironmentTheCounterStream()
        #env = EnvironmentTheStream()
        env.DisplayForces()
        #######################################################################

        robot = env.GetRobot()
        env.MakeRobotInvisible()

        p = np.array((0,0,0,-pi/2))




