#!/usr/bin/env python

from __future__ import print_function

__author__ = 'yu-chi'

# System Imports

# 3rd-Party Imports
import numpy as np
import math
import heapq
import openravepy as rave
import copy
import random
import collections
import time

# Local Imports
import load_escher as load
from environment_handler_test import *

robot_z = 0.9

def waitrobot(robot):
    """busy wait for robot completion"""
    while not robot.GetController().IsDone():
        time.sleep(0.01)

def main():

    env_handler = environment_handler()
    env = env_handler.env
    escher = load.escher(env)

    l_arm_indices = escher.robot.GetManipulator('l_arm').GetArmIndices()
    r_arm_indices = escher.robot.GetManipulator('r_arm').GetArmIndices()
    l_leg_indices = escher.robot.GetManipulator('l_leg').GetArmIndices()
    r_leg_indices = escher.robot.GetManipulator('r_leg').GetArmIndices()

    additional_active_DOFs = ['x_prismatic_joint','y_prismatic_joint','z_prismatic_joint','roll_revolute_joint','pitch_revolute_joint'
                             ,'yaw_revolute_joint','waist_yaw']
   
    additional_active_DOF_indices = [None]*len(additional_active_DOFs)
    for index,j in enumerate(additional_active_DOFs):
        additional_active_DOF_indices[index] = escher.robot.GetJoint(j).GetDOFIndex()
    
    whole_body_indices = np.concatenate((l_arm_indices, r_arm_indices, l_leg_indices, r_leg_indices, additional_active_DOF_indices),axis=1)

    escher.robot.SetActiveDOFs(whole_body_indices)

    

    # initialize robot configuration
    escher.robot.SetTransform(np.array([[1,0,0,0],[0,1,0,0],[0,0,1,robot_z],[0,0,0,1]]))
    
    # surrender poseture
    DOFValues = escher.robot.GetDOFValues()
    DOFValues[6] = math.pi/2
    DOFValues[19] = math.pi/2
    DOFValues[24] = -math.pi/2
    DOFValues[37] = -math.pi/2
    escher.robot.SetDOFValues(DOFValues)

    OriginalDOFValues = escher.robot.GetDOFValues()

    rave_traj = rave.RaveCreateTrajectory(escher.env,'')
    traj_serial = open('trajectory.xml','r').read()
    rave_traj.deserialize(traj_serial)

    rave.planningutils.RetimeActiveDOFTrajectory(rave_traj,escher.robot)

    escher.robot.GetController().SetPath(rave_traj)
    waitrobot(escher.robot)

    import IPython; IPython.embed()

if __name__ == "__main__":
    main()
