
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
        robot = env.GetRobot()
        env.MakeRobotInvisible()
        env.DisplayForces()
        #time.sleep(0.5)

        #traj = Trajectory.from_file('trajectories/clc_error')
        #traj.draw(env)
        #traj.draw_robot_along_path(env, robot, N=15)
        #######################################################################

        W = np.loadtxt("topp/segfault_W1")
        traj = Trajectory(W)
        traj.draw(env)
        #time.sleep(2)
                        #print "ECOS dt",dt,"dp2x",np.linalg.norm(p-qcontrol),"dx2pnext",np.linalg.norm(pnext-qcontrol),"qcontrol",qcontrol,"p",p,"pnext",pnext
        Nc = traj.getCriticalPoint(env)


        #traj = Trajectory.from_file('trajectories/bloodstream_kinodynamic')
        #traj = Trajectory.from_file('trajectories/bloodstream_kinodynamic_deformed')
        #traj.draw(env)
        #traj.draw_robot_along_path(env, robot, N=15)
        #traj.PlotParametrization(env)
        #traj.execute(env, robot, tsleep=0.003, stepping=False)
