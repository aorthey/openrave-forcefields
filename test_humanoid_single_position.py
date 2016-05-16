#start!/usr/bin/env python
import time
import scipy
import sys
import numpy as np
import openravepy
from openravepy import *
from math import *
from environment_force_humanoid import *
from environment_force_humanoid_contact_world import *
from cbirrtpy import *

from deformation_naive import *
from deformation_potentials import *
from deformation_stretchpull import *
from trajectory_bspline import *
import numpy as np
from motion_planner_geometrical import MotionPlannerGeometrical
from motion_planner_kinodynamic import MotionPlannerKinodynamic
from util_humanoid import *
from gik_interface import *

#sys.path.append(os.environ["MPP_PATH"]+"mpp-robot/mpp")
#sys.path.append(os.environ["OPENRAVE_WPI_PATH"]+"ros_msgs")
#from SurfaceMeshes import *


def waitrobot(robot):
    """busy wait for robot completion"""
    while not robot.GetController().IsDone():
        time.sleep(0.1)

if __name__ == "__main__":

        env = EnvironmentHumanoidContactWorld()
        env.DrawAxes()
        robot = env.GetRobot()
        surfaces = env.GetSurfaces()
        time.sleep(0.5)

        print surfaces.shape
        for i in arange(surfaces.shape[0]):
                et = surfaces[i,4,0]
                eo = surfaces[i,5,0]
                print "Surface",i,
                if et <= 0.01:
                        print "no extension",
                if eo <= 0.01:
                        print "no extension",
                print

        from surface_module import *
        M = surfaces.shape[0]

        S = SurfaceModule(surfaces)
        robot = env.GetRobot()

        com = robot.GetCenterOfMass()

        left_leg_tf = robot.GetManipulator('l_leg').GetTransform()
        right_leg_tf = robot.GetManipulator('r_leg').GetTransform()
        left_arm_tf = robot.GetManipulator('l_arm').GetTransform()
        right_arm_tf = robot.GetManipulator('r_arm').GetTransform()

        left_leg_tf = S.GetNearestContactTransform(env, left_leg_tf, 4)
        right_leg_tf = S.GetNearestContactTransform(env, right_leg_tf, 4)
        left_arm_tf = S.GetNearestContactTransform(env, left_arm_tf, 15)
        right_arm_tf = S.GetNearestContactTransform(env, right_arm_tf, 26)

        env.DrawFootContactPatchFromTransform(left_leg_tf)
        env.DrawFootContactPatchFromTransform(right_leg_tf)
        env.DrawFootContactPatchFromTransform(left_arm_tf)
        env.DrawFootContactPatchFromTransform(right_arm_tf)

        S.SampleSurface(100,26,env)
        S.SampleSurface(100,15,env)
        S.SampleSurface(100,4,env)

        gik = GIKInterface(env)
        q = gik.fromContactTransform(robot, left_leg_tf, right_leg_tf, None, None)

        #raw_input('Press <ENTER> to execute trajectory.')
        robot.WaitForController(0)
        #env.DisplayForces()
        #time.sleep(0.5)


