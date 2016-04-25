#!/usr/bin/env python
import time
import scipy
import sys
import numpy as np
import openravepy
from openravepy import *
from math import *
from environment_force_humanoid import *
from cbirrtpy import *

from deformation_naive import *
from deformation_potentials import *
from deformation_stretchpull import *
from trajectory_bspline import *
import numpy as np
from motion_planner_geometrical import MotionPlannerGeometrical
from motion_planner_kinodynamic import MotionPlannerKinodynamic
from util_humanoid import *
#import statsmodels.api as sm

def waitrobot(robot):
    """busy wait for robot completion"""
    while not robot.GetController().IsDone():
        time.sleep(0.1)

if __name__ == "__main__":

        env = EnvironmentHumanoid()
        env.DisplayForces()
        time.sleep(0.5)
        robot = env.GetRobot()

        #######################################################################
        ### READ GEOMETRICAL PATH FROM XML
        #######################################################################
        tname = 'misc/trajectory.xml'

        fhandle = open(tname, 'r')
        trajdata = fhandle.read()
        fhandle.close()
        rave_path = RaveCreateTrajectory(env.env,'').deserialize(trajdata)

        rave.planningutils.RetimeActiveDOFTrajectory(rave_path,robot)
        rave.planningutils.SmoothActiveDOFTrajectory(rave_path,robot)

        cbirrt = CBiRRT(env.env, env.robot_name)
        M = rave_path.GetNumWaypoints()
        N = len(robot.GetActiveDOFValues())

        [q_original, COM_original] = COM_from_path( rave_path, robot, env)


        #######################################################################
        ### COMPUTE NEW COM
        #######################################################################

        T = arange(0,M)
        C = float(M/12)
        COM_offset = np.zeros((3,M))
        COM_offset[1,:]=-0.35*exp(-(T-M/2)**2/(2*C*C))
        COM_offset[2,:]=-0.05*exp(-(T-M/2)**2/(2*C*C))

        tmp_handle=[]
        with env.env:
                h=env.env.drawlinestrip(points=COM_original.T,linewidth=6,colors=np.array((1,0,0)))
                tmp_handle.append(h)
                h=env.env.drawlinestrip(points=COM_offset.T,linewidth=6,colors=np.array((0,1,0)))
                tmp_handle.append(h)


        #######################################################################
        ### RECOMPUTE GIK FROM NEW COM
        #######################################################################

        [q_gik, COM_gik] = GIK_from_COM( COM_offset, q_original, robot, env, recompute=False)

        tmp_handle=[]
        with env.env:
                h=env.env.drawlinestrip(points=COM_original.T,linewidth=6,colors=np.array((1,0,0)))
                tmp_handle.append(h)
                h=env.env.drawlinestrip(points=COM_gik.T,linewidth=6,colors=np.array((0,1,0)))
                tmp_handle.append(h)
                h=env.env.drawlinestrip(points=COM_offset.T,linewidth=6,colors=np.array((0.8,0,0.8,0.3)))
                tmp_handle.append(h)
                robot.SetActiveDOFValues(q_gik[:,0])

        ###############################################################################
        ##### CREATE A NEW RAVE TRAJ FROM NEW CONFIGURATIONS AND RETIME
        ###############################################################################

        traj = RaveCreateTrajectory(env.env,'')
        traj.Init(robot.GetActiveConfigurationSpecification())
        i=0
        while i < M:
                traj.Insert(i,q_gik[:,i])
                i=i+1

        #rave.planningutils.RetimeActiveDOFTrajectory(traj,robot)
        rave.planningutils.RetimeActiveDOFTrajectory(traj,robot)#plannername='ParabolicTrajectoryRetimer')
        rave.planningutils.SmoothActiveDOFTrajectory(traj,robot)

        raw_input('Press <ENTER> to start.')
        openravepy.RaveLogInfo("Waiting for controller to finish")
        #robot.GetController().SetPath(rave_path)
        robot.GetController().SetPath(traj)
        robot.WaitForController(0)
        robot.GetController().Reset()
        raw_input('Enter any key to quit. ')
