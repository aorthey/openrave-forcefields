#start!/usr/bin/env python
import os
import sys
sys.path.append(os.environ["OPENRAVE_WPI_PATH"])
import time
import scipy
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

def waitrobot(robot):
    """busy wait for robot completion"""
    while not robot.GetController().IsDone():
        time.sleep(0.1)

if __name__ == "__main__":

        env = EnvironmentHumanoidContactWorld()


        env.DrawAxes()
        robot = env.GetRobot()
        env.MakeRobotTransparent(alpha=1.0)
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

        left_leg_surface = 4
        right_leg_surface = 4
        left_hand_surface = 13
        right_hand_surface = 26
        #right_hand_surface = 19
        #S.SampleSurface(100,19,env)

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

        #env.DrawFootContactPatchFromTransform(left_leg_tf)
        #env.DrawFootContactPatchFromTransform(right_leg_tf)
        #env.DrawLeftHandContactPatchFromTransform(left_arm_tf)
        #env.DrawRightHandContactPatchFromTransform(right_arm_tf)

        #S.SampleSurface(100,right_hand_surface,env)
        #S.SampleSurface(100,left_hand_surface,env)
        #sys.exit(0)
        with env.env:
                robot.GetLinks()[0].SetStatic(True)
                W = env.env.GetKinBody('world')
                L = W.GetLinks()
                F = np.array((-1,0,0))
                for link in L:
                        P = link.GetGlobalCOM()
                        link.SetStatic(True)
                        link.SetForce(F,P,False)
                #self.env.GetPhysicsEngine().SetBodyForce(link,F,P,True)
        np.random.seed(2)
        gik = GIKInterface(env)
        
        #raw_input('Press <ENTER> to sample stances.')
        q = None
        while q is None:
                with env.env:
                        env.env.StopSimulation()

                        while q is None:
                                #left_leg_tf = S.SampleContactFoot(left_leg_surface, env)
                                #right_leg_tf = S.SampleContactFoot(right_leg_surface, env)
                                left_arm_tf = S.SampleContactLeftHand(left_hand_surface, env)
                                #right_arm_tf = S.SampleContactRightHand(right_hand_surface, env)
                                #q = gik.fromContactTransform(robot, left_leg_tf, right_leg_tf, left_arm_tf, right_arm_tf)
                                q = gik.fromContactTransform(robot, left_leg_tf, right_leg_tf, left_arm_tf, None)

                                red=np.array((1,0,0))
                                green=np.array((0,1,0))
                                if q is None:
                                        #env.DrawFootContactPatchFromTransform(left_leg_tf, cpatch=red)
                                        #env.DrawFootContactPatchFromTransform(right_leg_tf, cpatch=red)
                                        env.DrawLeftHandContactPatchFromTransform(left_arm_tf, cpatch=red)
                                        env.DrawRightHandContactPatchFromTransform(right_arm_tf, cpatch=red)
                                else:
                                        #env.DrawFootContactPatchFromTransform(left_leg_tf, cpatch=green)
                                        #env.DrawFootContactPatchFromTransform(right_leg_tf, cpatch=green)
                                        env.DrawLeftHandContactPatchFromTransform(left_arm_tf, cpatch = green)
                                        env.DrawRightHandContactPatchFromTransform(right_arm_tf, cpatch=green)
                        q = None
                        env.env.StartSimulation(0.000000001)
                time.sleep(0.0001)
                #env.env.StepSimulation(0.01)

        #env.env.StopSimulation() 
        gik.VisualizeFrictionCone(robot)
        #env.env.StartSimulation(0.000005)
        #env.env.StartSimulation(0.0000000001)

        while True:
                torques = 0*(numpy.random.rand(robot.GetDOF())-0.5)
                for i in range(100):
                        robot.SetJointTorques(torques,True)
                        time.sleep(0.01)


        robot.WaitForController(0)
