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
from environment_force_radial import *
from environment_force_blood_stream import *
from environment_force_blood_stream2 import *
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
        env = EnvironmentBloodStream2()
        traj = Trajectory.from_file('trajectories/bloodstream_geometric70')
        #######################################################################

        #######################################################################
        #env = EnvironmentRadial()
        #traj = Trajectory.from_file('trajectories/bloodstream_geometric80')
        #######################################################################

        #######################################################################
        #env = EnvironmentBloodStream()
        #traj = Trajectory.from_file('trajectories/bloodstream_geometric30')
        #traj = Trajectory.from_file('trajectories/the_blood_stream_deformed')
        #traj = Trajectory.from_file('trajectories/epsilon_loop_projectable')
        #######################################################################

        robot = env.GetRobot()
        env.MakeRobotInvisible()
        env.DisplayForces()
        time.sleep(0.5)

        #traj = Trajectory.from_file('trajectories/bloodstream_kinodynamic')
        #traj = Trajectory.from_file('trajectories/bloodstream_kinodynamic_deformed')
        traj.draw(env)
        traj.draw_robot_along_path(env, robot, N=15)
        #traj.PlotParametrization(env)
        #traj.execute(env, robot, tsleep=0.003, stepping=False)
