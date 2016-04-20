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

        cbirrt = CBiRRT(env.env, env.robot_name)

        raw_input('Press <ENTER> to start.')
        #robot.GetController().SetPath(rave_path)
        #waitrobot(robot)
        #traj = TrajectoryHumanoid.from_ravetraj(rave_path)
        N = len(robot.GetActiveDOFValues())
        active_dofs = robot.GetActiveConfigurationSpecification()

        M = rave_path.GetNumWaypoints()
        COM_original = np.zeros((3,M))
        COM_gik = np.zeros((3,M))
        q_gik = np.zeros((N, M))

        i = 0
        while i < N:
                q = rave_path.GetWaypoint(i,active_dofs)
                robot.SetActiveDOFValues(q)

                left_leg_tf = robot.GetManipulator('l_leg').GetTransform()
                right_leg_tf = robot.GetManipulator('r_leg').GetTransform()
                left_arm_tf = robot.GetManipulator('l_arm').GetTransform()
                right_arm_tf = robot.GetManipulator('r_arm').GetTransform()

                maniptm_list = [('l_leg',left_leg_tf),('r_leg',right_leg_tf),('l_arm',left_arm_tf),('r_arm',right_arm_tf)]

                friction_coefficient = 0.8

                #support_list = [('l_leg',friction_coefficient),('r_leg',friction_coefficient),('l_arm',friction_coefficient),('r_arm',friction_coefficient)]

                support_list = [('l_leg',friction_coefficient),('r_leg',friction_coefficient)]

                COM_original[:,i] = robot.GetCenterOfMass()

                q_gik[:,i] = cbirrt.DoGeneralIK(execute=False,
                                returnclosest=True, maniptm=maniptm_list,
                                support=support_list, movecog=COM_original[:,i],
                                printcommand=False)
                robot.SetActiveDOFValues(q_gik[:,i])
                COM_gik[:,i] = robot.GetCenterOfMass()
                print i,":",COM_original[:,i],"<->",COM_gik[:,i]

                i = i+1

        tmp_handle=[]
        with env.env:
                h=env.env.drawlinestrip(points=COM_gik.T,linewidth=6,colors=np.array((1,0,0)))
                tmp_handle.append(h)
        i=0
        while i < N:
                robot.SetActiveDOFValues(q_gik[:,i])
                F = np.array([0.0,0.0000001,0])
                for link in robot.GetLinks():
                        P = link.GetLocalCOM()
                        Pg = link.GetGlobalCOM()
                        #print P,Pg,link.GetName(),link.IsEnabled()
                        #print link.GetGlobalInertia(), 
                        print link.GetLocalInertia()
                        if Pg[0] > 1.0:
                                link.SetForce(F,P,True)
                        #sys.exit(0)

                #print x,y,z
                #sys.exit(0)
                #if x > 1.0:
                        #env.env.GetPhysicsEngine().SetGravity([0,0.2,0])
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
