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

def waitrobot(robot):
    """busy wait for robot completion"""
    while not robot.GetController().IsDone():
        time.sleep(0.01)

if __name__ == "__main__":

        env = EnvironmentHumanoid()
        env.DisplayForces()
        #env.MakeRobotVisible()
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
        #env.env.GetPhysicsEngine().SetGravity([0,0,-9.81])
        #env.env.GetPhysicsEngine().SetGravity([0,5.0,0])
        rave.planningutils.RetimeActiveDOFTrajectory(rave_path,robot)

        raw_input('Press <ENTER> to start.')
        #robot.GetController().SetPath(rave_path)
        #waitrobot(robot)

        N = rave_path.GetNumWaypoints()
        i = 0
        while i < N:
                q = rave_path.GetWaypoint(i,robot.GetConfigurationSpecification())
                robot.SetConfigurationValues(q)
                x=robot.GetJoint('x_prismatic_joint').GetValues()
                y=robot.GetJoint('y_prismatic_joint').GetValues()
                z=robot.GetJoint('z_prismatic_joint').GetValues()
                print x,y,z
                if x > 1.0:
                        env.env.GetPhysicsEngine().SetGravity([0,0.2,0])
                #env.env.StepSimulation(0.0001)
                waitrobot(robot)
                time.sleep(0.1)
                i = i+1

        #trajectory = MotionPlannerDeformation(path, robot, env)
        #traj = Trajectory.from_ravetraj(rave_path)
        #traj.info()
        #traj.draw(env)
        #traj.PlotParametrization(env)
        #RaveSetDebugLevel(DebugLevel.Debug) # set output level to debug

        #openravepy.RaveLogInfo("Waiting for controller to finish")
        #robot.GetController().SetPath(rave_path)
        #robot.WaitForController(0)
        ##robot.GetController().Reset()

        #raw_input('Enter any key to quit. ')
