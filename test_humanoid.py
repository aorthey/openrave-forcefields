#!/usr/bin/env python
import time
import scipy
import sys
import os.path
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
        time.sleep(0.1)

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
        rave.planningutils.SmoothActiveDOFTrajectory(rave_path,robot)

        cbirrt = CBiRRT(env.env, env.robot_name)
        #robot.GetController().SetPath(rave_path)
        #waitrobot(robot)
        #from trajectory_humanoid import *
        #traj = TrajectoryHumanoid.from_ravetraj(rave_path, robot)
        #traj.computeTrajectoryStringForTOPP(traj.waypoints)

        N = len(robot.GetActiveDOFValues())
        active_dofs = robot.GetActiveConfigurationSpecification()

        M = rave_path.GetNumWaypoints()
        COM_original = np.zeros((3,M))
        COM_gik = np.zeros((3,M))

        q_original = np.zeros((N, M))
        q_gik = np.zeros((N, M))

        i = 0
        [qlimL,qlimU]=robot.GetActiveDOFLimits()

        print "Waypoints:",M," - Dimension:",N
        while i < M:
                q = rave_path.GetWaypoint(i,active_dofs)
                q_original[:,i] = env.EnforceLimits(q,qlimL,qlimU)
                robot.SetActiveDOFValues(q_original[:,i])
                COM_original[:,i] = robot.GetCenterOfMass()
                i = i+1

        q_gik_fname = 'tmp/q_gik.numpy'
        COM_gik_fname = 'tmp/COM_gik.numpy'

        T = arange(0,M)
        C = float(M/12)
        COM_offset = np.zeros((3,M))
        COM_offset[1,:]=-0.30*exp(-(T-M/2)**2/(2*C*C))
        COM_offset[2,:]=-0.05*exp(-(T-M/2)**2/(2*C*C))

        def avalue(Ncritical, i, c=k*20.0):
                return np.exp(-((Ncritical-i)*(Ncritical-i))/(2*c*c))

        def Amatrix(i, N):
                A = np.zeros(N)
                assert(i<N)

                j = N
                while j > 0:
                        if j<N-1:
                                A[j] = avalue(i, j, 10.0)
                        j -= 1
                return A

        #COM_offset_tmp = np.zeros((3,M))
        #print Amatrix(M/2,M)
        #sys.exit(0)
        #for i in range(0,M):
        #        A = Amatrix(i,M)
        #        COM_offset_tmp[0,i] = np.dot(A,COM_offset[0,:])
        #        COM_offset_tmp[1,i] = np.dot(A,COM_offset[1,:])
        #        COM_offset_tmp[2,i] = np.dot(A,COM_offset[2,:])
        #COM_offset = COM_offset_tmp

        tmp_handle=[]
        with env.env:
                h=env.env.drawlinestrip(points=COM_original.T,linewidth=6,colors=np.array((1,0,0)))
                tmp_handle.append(h)
                h=env.env.drawlinestrip(points=COM_offset.T,linewidth=6,colors=np.array((0,1,0)))
                tmp_handle.append(h)
                #h=env.env.drawlinestrip(points=COM_offset_tmp.T,linewidth=6,colors=np.array((0,0,1)))
                #tmp_handle.append(h)


        #if True:
        if not os.path.isfile(q_gik_fname+'.npy'):
                i = 0
                while i < M:
                        print "------------------------------------------------------------------"
                        print "------------ WAYPOINT",i,"/",M,": recompute GIK ------------ "
                        print "------------------------------------------------------------------"
                        try:
                                with env.env:
                                        #q_original[:,i] = env.EnforceLimits(q_original[:,i],qlimL,qlimU,DEBUG=True)
                                        #robot.SetActiveDOFValues(q_original[:,i])
                                        #waitrobot(robot)
                                        robot.SetActiveDOFValues(q_original[:,i])

                                        left_leg_tf = robot.GetManipulator('l_leg').GetTransform()
                                        right_leg_tf = robot.GetManipulator('r_leg').GetTransform()
                                        left_arm_tf = robot.GetManipulator('l_arm').GetTransform()
                                        right_arm_tf = robot.GetManipulator('r_arm').GetTransform()

                                        #maniptm_list = [('l_leg',left_leg_tf),('r_leg',right_leg_tf),('l_arm',left_arm_tf),('r_arm',right_arm_tf)]
                                        maniptm_list = [('l_leg',left_leg_tf),('r_leg',right_leg_tf)]
                                        friction_coefficient = 0.8
                                        support_list = [('l_leg',friction_coefficient),
                                                        ('r_leg',friction_coefficient),
                                                        ('l_arm',friction_coefficient),
                                                        ('r_arm',friction_coefficient)]

                                        COM_original[:,i] = robot.GetCenterOfMass()

                                        cog = COM_offset[:,i]
                                        #cog = COM_original[:,i]
                                        q_res = cbirrt.DoGeneralIK(
                                                        movecog=cog,
                                                        maniptm=maniptm_list,
                                                        support=support_list,
                                                        printcommand=True)

                                        if q_res is None:
                                                print "------------ WAYPOINT",i,"/",M,": recompute GIK ------------ "
                                                print "No solution found GIK"
                                                sys.exit(0)
                                        else:
                                                q_gik[:,i] = q_res
                                        print "DIST q,q_gik:",np.linalg.norm(q_original[:,i]-q_gik[:,i])
                                        robot.SetActiveDOFValues(q_gik[:,i])
                                        COM_gik[:,i] = robot.GetCenterOfMass()

                        except Exception as e:
                                print "Exception in GIK, waypoint",i,"/",M
                                print e
                                print "q_original[:,i]:",q_original[:,i]
                                epsilon=1e-3
                                q = q_original[:,i]
                                for j in range(0,q.shape[0]):
                                        if ( q[j] <= qlimL[j] + epsilon ):
                                                qq = q[j] - epsilon
                                                q[j] = qlimL[j] + epsilon
                                                print "q[",j,"]=",qq," (<",qlimL[j],") => q[j]=",q[j]
                                        if ( q[j] >= qlimU[j] - epsilon ):
                                                qq = q[j] + epsilon
                                                q[j] = qlimU[j] - epsilon
                                                print "q[",j,"]=",qq," (>",qlimU[j],") => q[j]=",q[j]
                                sys.exit(0)

                        dcc = np.linalg.norm(COM_original[:,i]-COM_gik[:,i])
                        print i,"/",M," :",dcc,"(",COM_original[:,i],"<->",COM_gik[:,i],")"
                        print "ORIGINAL   COM :",COM_original[:,i]
                        print "INPUT GIK  COM :",cog
                        print "OUTPUT GIK COM :",COM_gik[:,i]
                        i = i+1
                        #sys.exit(0)

                np.save(q_gik_fname,q_gik)
                np.save(COM_gik_fname,COM_gik)
        else:
                q_gik = np.load(q_gik_fname+'.npy')
                COM_gik = np.load(COM_gik_fname+'.npy')

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
##### CREATE NEW COM PATH
###############################################################################

        traj = RaveCreateTrajectory(env.env,'')
        traj.Init(robot.GetActiveConfigurationSpecification())
        #traj.Insert(0,robot.GetActiveDOFValues())

        i=0
        while i < M:
                traj.Insert(i,q_gik[:,i])
                i=i+1

        #rave.planningutils.RetimeActiveDOFTrajectory(traj,robot)
        rave.planningutils.RetimeActiveDOFTrajectory(traj,robot,hastimestamps=False,maxvelmult=1)#plannername='ParabolicTrajectoryRetimer')
        rave.planningutils.SmoothActiveDOFTrajectory(traj,robot)

        raw_input('Press <ENTER> to start.')
        openravepy.RaveLogInfo("Waiting for controller to finish")
        robot.GetController().SetPath(traj)
        #robot.GetController().SetPath(rave_path)
        robot.WaitForController(0)
        robot.GetController().Reset()
   
        #while i < M:
                #print "DIST q,q_gik:",np.linalg.norm(q_original[:,i]-q_gik[:,i])
                #robot.SetActiveDOFValues(q_gik[:,i])
                #robot.WaitForController(0)
                #waitrobot(robot)
                #time.sleep(0.1)
                #robot.SetActiveDOFValues(q_original[:,i])
                #waitrobot(robot)
                #time.sleep(0.1)

                #print "OUTPUT GIK COM :",COM_gik[:,i]
                #F = np.array([0.0,0.0000001,0])
                #for link in robot.GetLinks():
                #        P = link.GetLocalCOM()
                #        Pg = link.GetGlobalCOM()
                #        #print P,Pg,link.GetName(),link.IsEnabled()
                #        #print link.GetGlobalInertia(), 
                #        #print link.GetLocalInertia()
                #        if Pg[0] > 1.0:
                #                link.SetForce(F,P,True)
                #        #sys.exit(0)

                #print x,y,z
                #sys.exit(0)
                #if x > 1.0:
                        #env.env.GetPhysicsEngine().SetGravity([0,0.2,0])
                #env.env.StepSimulation(0.1)
                #time.sleep(0.1)
                #i = i+1

        #trajectory = MotionPlannerDeformation(path, robot, env)
        #traj = Trajectory.from_ravetraj(rave_path)
        #traj.info()
        #traj.draw(env)
        #traj.PlotParametrization(env)
        #RaveSetDebugLevel(DebugLevel.Debug) # set output level to debug

        #openravepy.RaveLogInfo("Waiting for controller to finish")
        #robot.GetController().SetPath(rave_path)
        #robot.WaitForController(0)
        #robot.GetController().Reset()

        #raw_input('Enter any key to quit. ')
