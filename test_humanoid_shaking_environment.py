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
from util import *
from surface_module import *

#sys.path.append(os.environ["MPP_PATH"]+"mpp-robot/mpp")
#sys.path.append(os.environ["OPENRAVE_WPI_PATH"]+"ros_msgs")
#from SurfaceMeshes import *


def waitrobot(robot):
    """busy wait for robot completion"""
    while not robot.GetController().IsDone():
        time.sleep(0.1)

if __name__ == "__main__":

        np.random.seed(2)
        env = EnvironmentHumanoidContactWorld('environments/contactworld2.env.xml')
        env.DrawAxes()
        #env.MakeRobotTransparent(alpha=1.0)
        time.sleep(0.5)

        #######################################################################
        gik = GIKInterface(env)
        robot = env.GetRobot()
        #######################################################################

        left_leg_surface = 4
        right_leg_surface = 4
        left_hand_surface = 13
        right_hand_surface = 26

        #######################################################################
        left_leg_tf = robot.GetManipulator('l_leg').GetTransform()
        right_leg_tf = robot.GetManipulator('r_leg').GetTransform()
        left_arm_tf = robot.GetManipulator('l_arm').GetTransform()
        right_arm_tf = robot.GetManipulator('r_arm').GetTransform()

        #raw_input('Press <ENTER> to sample stances.')
        with env.env:
                env.env.StopSimulation()
                robot.GetLinks()[0].SetStatic(True)
                W = env.env.GetKinBody('world')
                env.env.GetPhysicsEngine().SetGravity(np.array((0,0,-0.5)))

        surfaces = env.GetSurfaces()
        S = SurfaceModule(surfaces)

        left_leg_tf = S.GetNearestContactTransformFoot(env, 
                        left_leg_tf, 
                        left_leg_surface)
        right_leg_tf = S.GetNearestContactTransformFoot(env, 
                        right_leg_tf, 
                        right_leg_surface)
        left_arm_tf = S.GetNearestContactTransformLeftHand(env, 
                        left_arm_tf,
                        left_hand_surface)
        right_arm_tf = S.GetNearestContactTransformRightHand(env, 
                        right_arm_tf,
                        right_hand_surface)

        left_arm_relative_surface_tf = S.GetContactRelativeToSurfaceTransform(left_arm_tf, left_hand_surface) 
        right_arm_relative_surface_tf = S.GetContactRelativeToSurfaceTransform(right_arm_tf, right_hand_surface) 
        left_leg_relative_surface_tf = S.GetContactRelativeToSurfaceTransform(left_leg_tf, left_leg_surface) 
        right_leg_relative_surface_tf = S.GetContactRelativeToSurfaceTransform(right_leg_tf, right_leg_surface) 

        ictr=0
        while True:
                with env.env:
                        #env.AddRandomDisturbance()
                        if ictr%20==0:
                                maxforce = 10
                                r = (numpy.random.rand(2)-0.5)
                                F = np.array((r[0]*maxforce,r[1]*maxforce,0))

                        env.AddForceToWorld(F)

                        #######################################################################
                        surfaces = env.GetSurfaces()
                        S = SurfaceModule(surfaces)
                        #######################################################################

                        #######################################################
                        ### draw
                        #######################################################
                        if env.static_handles:
                                env.static_handles = None
                                env.static_handles = []

                        left_arm_tf = S.AdjustToNewSurfaceTransform( 0.05, left_arm_relative_surface_tf, left_hand_surface, env)
                        right_arm_tf = S.AdjustToNewSurfaceTransform( 0.05, right_arm_relative_surface_tf, right_hand_surface, env)
                        left_leg_tf = S.AdjustToNewSurfaceTransform( 0.05, left_leg_relative_surface_tf, left_leg_surface, env)
                        right_leg_tf = S.AdjustToNewSurfaceTransform( 0.05, right_leg_relative_surface_tf, right_leg_surface, env)

                        env.DrawFootContactPatchFromTransform(left_leg_tf, cpatch=green)
                        env.DrawFootContactPatchFromTransform(right_leg_tf, cpatch=green)
                        env.DrawLeftHandContactPatchFromTransform(left_arm_tf, cpatch = green)
                        S.DrawCoordinateFrames(env)

                        q = gik.fromContactTransform(robot, left_leg_tf, right_leg_tf, left_arm_tf, None)

                        #torques = (numpy.random.rand(robot.GetDOF())-0.5)
                        #robot.SetJointTorques(torques,True)
                        #######################################################
                        ### random noise on environment
                        #######################################################

                env.env.StepSimulation(0.001)
                ictr+=1

        #gik.VisualizeFrictionCone(robot)
        #while True:
                #torques = 0*(numpy.random.rand(W.GetDOF())-0.5)
                #W.SetJointTorques(torques,False)
                #time.sleep(0.01)

        #while True:
                #torques = 0*(numpy.random.rand(robot.GetDOF())-0.5)
                #for i in range(100):
                        #robot.SetJointTorques(torques,True)
                        #time.sleep(0.01)


        robot.WaitForController(0)
