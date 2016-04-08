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
from environment_periodic_force_the_hideout import *
from environment_periodic_force_triple_stream import *
from environment_periodic_force_crossroad_stream import *

from deformation_naive import *
from deformation_potentials import *
from deformation_stretchpull import *
from trajectory_bspline import *
import numpy as np
from motion_planner_geometrical import MotionPlannerGeometrical
from motion_planner_kinodynamic import MotionPlannerKinodynamic
#import statsmodels.api as sm

if __name__ == "__main__":

        #######################################################################
        env = EnvironmentTheRay()
        #env = EnvironmentTheCounterStream()
        #env = EnvironmentTheStream()
        #######################################################################

        robot = env.GetRobot()
        env.MakeRobotInvisible()
        env.DisplayForces()
        time.sleep(0.5)

        planner = MotionPlannerGeometrical(robot, env)
        #planner = MotionPlannerKinodynamic(robot, env)

        rave_path = planner.GetPath()

        #trajectory = MotionPlannerDeformation(path, robot, env)
        traj = Trajectory.from_ravetraj(rave_path)
        traj.info()
        traj.draw(env)
        #traj.PlotParametrization(env)
        traj.draw_delete()
        td = DeformationStretchPull(traj, env)

        Nd = 5
        raw_input('Press <ENTER> to start.')


        td.deform(N_iter=100)
        #for i in range(Nd):
        #        print "DEFORMATION:",i,"/",Nd
        #        if td.deform(N_iter=10):
        #                td.draw_deformation() 
        #        else:
        #                td.draw_deformation() 
        #                break

        xt = td.traj_current.topp.traj0
        env.MakeRobotVisible()

        #raw_input('Press <ENTER> to execute trajectory.')
        #RaveSetDebugLevel(DebugLevel.Debug) # set output level to debug
        #openravepy.RaveLogInfo("Waiting for controller to finish")
        #robot.GetController().SetPath(traj)
        #robot.WaitForController(0)
        #robot.GetController().Reset()
        raw_input('Enter any key to quit. ')

