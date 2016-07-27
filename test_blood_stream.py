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
from trajectory import *
import numpy as np
from motion_planner_geometrical import MotionPlannerGeometrical
from motion_planner_kinodynamic import MotionPlannerKinodynamic
#import statsmodels.api as sm
np.set_printoptions(precision=2)

if __name__ == "__main__":

        #######################################################################
        env = EnvironmentBloodStream()
        #env = EnvironmentTheCounterStream()
        #env = EnvironmentTheStream()
        #######################################################################

        robot = env.GetRobot()
        env.MakeRobotInvisible()
        env.DisplayForces()

        planner = MotionPlannerGeometrical(robot, env)
        ##planner = MotionPlannerKinodynamic(robot, env)
        rave_path = planner.GetPath()

        traj = Trajectory.from_ravetraj(rave_path)
        #traj = Trajectory.from_file('trajectories/bloodstream_kinodynamic1')
        #traj = Trajectory.from_file('deform1')

        traj.info()
        traj.draw(env)
        xml = env.GetName()
        #time.sleep(5.0)
        traj.draw_delete()
        traj.save('trajectories/'+xml)

        #traj.PlotParametrization(env)
        #traj.execute(env, robot, tsleep=0.005, stepping=False)
        #time.sleep(1)

        td = DeformationReachableSet(traj, env)
        deform_success = td.deform(N_iter=100)

        td.traj_deformed.save('trajectories/'+xml+'_deformed')
        if deform_success:
                #td.traj_deformed.PlotParametrization(env)
                print td.traj_deformed.waypoints.shape
                td.draw_delete()
                td.traj_deformed.draw(env,critical_pt=td.traj_deformed.waypoints.shape[1])
                td.traj_deformed.draw_robot_along_path(env, robot, N=6)
                #td.traj_deformed.execute(env, robot, tsleep=0.005, stepping=False)



