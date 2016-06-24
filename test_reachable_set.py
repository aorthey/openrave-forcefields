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
from trajectory_bspline import *
import numpy as np
from motion_planner_geometrical import MotionPlannerGeometrical
from motion_planner_kinodynamic import MotionPlannerKinodynamic
#import statsmodels.api as sm

if __name__ == "__main__":

        #######################################################################
        env = EnvironmentBloodStream()
        #######################################################################

        robot = env.GetRobot()

        planner = MotionPlannerGeometrical(robot, env)
        rave_path = planner.GetPath()

        #trajectory = MotionPlannerDeformation(path, robot, env)
        traj = Trajectory.from_ravetraj(rave_path)
        #traj.plot_reachableset(env)
        traj.test_domain_error(env)

        raw_input('Enter any key to quit. ')
