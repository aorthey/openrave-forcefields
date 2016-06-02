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
        traj.draw_delete()

        td = DeformationReachableSet(traj, env)

        Nd = 5
        #raw_input('Press <ENTER> to start.')

        td.deform(N_iter=100)
        td.traj_deformed.PlotParametrization(env)

        xt = td.traj_current.topp.traj0

        with env.env:
                robot.GetLinks()[0].SetStatic(True)
                env.env.StopSimulation() 

        t = 0.0
        tstep = 0.01

        robot.SetDOFValues(xt.Eval(t))
        env.MakeRobotVisible()

        while t < xt.duration:

                q = xt.Eval(t)
                dq = xt.Evald(t)
                ddq = xt.Evaldd(t)

                qn = q + tstep*dq + 0.5*tstep*tstep*ddq
                robot.SetDOFValues(qn)

                env.env.StepSimulation(tstep)
                time.sleep(0.01)
                t += tstep

        robot.WaitForController(0)

        #raw_input('Press <ENTER> to execute trajectory.')
        #RaveSetDebugLevel(DebugLevel.Debug) # set output level to debug
        #openravepy.RaveLogInfo("Waiting for controller to finish")
        #robot.GetController().SetPath(traj)
        #robot.WaitForController(0)
        #robot.GetController().Reset()
        raw_input('Enter any key to quit. ')
