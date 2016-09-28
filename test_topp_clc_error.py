
#!/usr/bin/env python
import time
import scipy
import sys
import numpy as np
import openravepy
from openravepy import *
from math import *
from environment_force_blood_stream import *

from deformation_reachableset import *
from trajectory_bspline import *
import numpy as np
from motion_planner_geometrical import MotionPlannerGeometrical
from motion_planner_kinodynamic import MotionPlannerKinodynamic
#import statsmodels.api as sm
if __name__ == "__main__":

        #######################################################################
        env = EnvironmentBloodStream()
        traj = Trajectory.from_file('trajectories/clc_error')
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
