#!/usr/bin/env python
import time
import scipy
import sys
import numpy as np
import openravepy
from openravepy import *
from math import *
from environment_force_humanoid import *

from deformation_naive import *
from deformation_potentials import *
from deformation_stretchpull import *
from trajectory_bspline import *
import numpy as np
from motion_planner_geometrical import MotionPlannerGeometrical
from motion_planner_kinodynamic import MotionPlannerKinodynamic
#import statsmodels.api as sm

if __name__ == "__main__":

        env = EnvironmentHumanoid()
        env.MakeRobotVisible()
        time.sleep(0.5)

        #planner = MotionPlannerGeometrical(robot, env)
        #planner = MotionPlannerKinodynamic(robot, env)

        #rave_path = planner.GetPath()
        tname = 'misc/trajectory.xml'

        fhandle = open(tname, 'r')
        trajdata = fhandle.read()
        fhandle.close()
        rave_path = RaveCreateTrajectory(env.env,'').deserialize(trajdata)

        robot = env.GetRobot()
        robot.GetController().SetPath(rave_path)

        N = rave_path.GetNumWaypoints()

        while not robot.GetController().IsDone():
                #robot.SetConfigurationValues(rave_path.GetWaypoint(i,robot.GetConfigurationSpecification()))
                env.env.StepSimulation(0.01)

        #trajectory = MotionPlannerDeformation(path, robot, env)
        #traj = Trajectory.from_ravetraj(rave_path)
        #traj.info()
        #traj.draw(env)
        #traj.PlotParametrization(env)
        RaveSetDebugLevel(DebugLevel.Debug) # set output level to debug

        openravepy.RaveLogInfo("Waiting for controller to finish")
        robot.GetController().SetPath(rave_path)
        robot.WaitForController(0)
        #robot.GetController().Reset()

        raw_input('Enter any key to quit. ')
