#start!/usr/bin/env python
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

        M = rave_path.GetNumWaypoints()
        N = len(robot.GetActiveDOFValues())

        [q_original, COM_original] = COM_from_path( rave_path, robot, env)

        #######################################################################
        ### COMPUTE NEW COM
        #######################################################################

        COM_linear = COM_interpolate(COM_original[:,0],COM_original[:,-1],M)

        T = arange(0,M)
        C = float(M/5)

        scale=0.01

        COM_offset = np.zeros((3,M))
        COM_offset += COM_linear
        COM_offset[0,:]+=0.00*exp(-(T-M/2)**2/(2*C*C))
        COM_offset[1,:]+=-scale*exp(-(T-M/2)**2/(2*C*C))
        COM_offset[2,:]+=-scale*0.1*exp(-(T-M/2)**2/(2*C*C))

        env.DrawAxes()

        tmp_handle=[]
        with env.env:
                h=env.env.drawlinestrip(points=COM_offset.T,linewidth=6,colors=np.array((0,1,0)))
                tmp_handle.append(h)
                [COM_zig_zag, footpos, dfootpos] = COM_compute_zig_zag_motion(COM_offset, env)

        time.sleep(0.1)

        #[q_gik, COM_gik] = GIK_from_COM_and_FOOTPOS( COM_zig_zag, footpos, dfootpos, robot, env, recompute=True)
        [q_gik, COM_gik] = GIK_from_COM_and_FOOTPOS( COM_zig_zag, footpos, dfootpos, robot, env)
        #[q_gik, COM_gik] = GIK_from_COM( COM_zig_zag, q_original, robot, env, recompute=True)

        #######################################################################
        ### RECOMPUTE GIK FROM NEW COM
        #######################################################################

        #q_gik_fname = 'tmp/q_gik.numpy'
        #COM_gik_fname = 'tmp/COM_gik.numpy'
        #q_gik = np.load(q_gik_fname+'.npy')
        #COM_gik = np.load(COM_gik_fname+'.npy')

        tmp_handle=[]
        with env.env:
                #h=env.env.drawlinestrip(points=COM_original.T,linewidth=6,colors=np.array((1,0,0)))
                #tmp_handle.append(h)
                h=env.env.drawlinestrip(points=COM_gik.T,linewidth=8,colors=np.array((0,1,0)))
                tmp_handle.append(h)
                h=env.env.drawlinestrip(points=COM_offset.T,linewidth=8,colors=np.array((0.8,0,0.8,0.3)))
                tmp_handle.append(h)
                robot.SetActiveDOFValues(q_gik[:,0])

        #time.sleep(0.5)
        #visualize_configurations(q_original, robot, env)
        #sys.exit(0)
        ###############################################################################
        ##### CREATE A NEW RAVE TRAJ FROM NEW CONFIGURATIONS AND RETIME
        ###############################################################################
        M = COM_gik.shape[1]
        print "WAYPOINTS: ",M
        traj = RaveCreateTrajectory(env.env,'')
        traj.Init(robot.GetActiveConfigurationSpecification())
        i=0
        while i < M:
                traj.Insert(i,q_gik[:,i])
                i=i+1

        Nl = len(robot.GetLinks())

        with env.env:
                #robot.GetLinks()[0].SetStatic(True)
                #env.env.GetPhysicsEngine().SetGravity([0,0,-9.81])
                env.env.GetPhysicsEngine().SetGravity([0,0,-0.05])
                #env.env.GetPhysicsEngine().SetGravity([0,0,0])
                rave.planningutils.SmoothActiveDOFTrajectory(traj,robot)
                rave.planningutils.RetimeActiveDOFTrajectory(traj,robot,hastimestamps=False,maxvelmult=0.75)
                robot.GetLinks()[0].SetStatic(True)
                #print robot.GetLinks()[0].GetName()
                #controller = RaveCreateController(env.env,'odevelocity')
                #robot.SetController(controller,range(robot.GetDOF()),0)

                active_dofs = robot.GetActiveConfigurationSpecification()

                print "Active DOF:",robot.GetActiveDOF()
                print "DOF       :",robot.GetDOF()
                print "traj DOF  :",q_gik[:,0].shape[0]

        Td = traj.GetDuration()
        forcetorquemap = {x: np.zeros(6) for x in np.arange(Nl)}
        ### create forcetorquemap dictionary
        link = robot.GetLink('torso')
        iTorso = link.GetIndex()
        Ftorso = np.array((0,1.0,0,0,0,0))
        Ffoot = np.array((0,0.0,1.0,0,0,0))
        forcetorquemap[iTorso]=Ftorso

        dt = 0.1
        t = 0.0
        qold = traj.Sample(0,active_dofs)
        qext = np.zeros((robot.GetDOF()))
        zlvec = []
        zrvec = []
        ictr=0

        while t < Td:
                qt = traj.Sample(t,active_dofs)

                with env.env:
                        ##### add contact forces
                        rlink = robot.GetManipulator('r_leg')
                        llink = robot.GetManipulator('l_leg')
                        zr = rlink.GetTransform()[2,3]
                        zl = llink.GetTransform()[2,3]
                        zlvec.append(zl)
                        zrvec.append(zr)

                        lsole = robot.GetLink('l_sole')
                        zsl = lsole.GetTransform()[2,3]
                        rsole = robot.GetLink('r_sole')
                        zsr = rsole.GetTransform()[2,3]

                        iLlink = lsole.GetIndex()
                        iRlink = rsole.GetIndex()

                        llink = robot.GetLink('l_foot')
                        rlink = robot.GetLink('l_foot')
                        env.DrawArrow(robot.GetLink('torso').GetGlobalCOM(), Ftorso[0:3], deleteOld=True)

                        if zl < 0.01:
                                forcetorquemap[iLlink]=Ffoot
                                env.DrawArrow(llink.GetGlobalCOM(), -Ffoot[0:3])
                        else:
                                forcetorquemap[iLlink]=np.zeros(6)
                                env.DrawArrow(llink.GetGlobalCOM(), np.zeros(3))
                        if zr < 0.01:
                                forcetorquemap[iRlink]=Ffoot
                        else:
                                forcetorquemap[iRlink]=np.zeros(6)

                        robot.SetDOFValues(qold)
                        tm,tc,tg = robot.ComputeInverseDynamics([],forcetorquemap,returncomponents=True)

                        qdd = (tm+tc+tg)
                        qext = qdd*dt*dt*0.5 + qext

                        #robot.SetDOFValues(qext)
                        #qext_active = robot.GetActiveDOFValues()

                        q = qt + qext

                        robot.SetDOFValues(q)
                        #qt2 = traj.Sample(t+dt,active_dofs)
                        #velocities = 0.01*(qt2-qt)/dt

                        #velocities = np.zeros(robot.GetDOF())
                        #robot.GetController().SendCommand('setvelocity '+' '.join(str(f) for f in velocities))
                        qold = q
                        ictr +=1
                        #print "time:",tend-tstart,"map:",tdyn-tstart,"invdyn:",tend-tdyn
                        #if ictr > 30:
                        #        plot(zlvec,'-b')
                        #        plot(zrvec,'-r')
                        #        plt.show()
                time.sleep(dt)
                #print "time:",tsleep-t0,"sleep:",tsleep-tend,"invdyn:",tend-tdyn,"tsample:",tstart-t0,"tmutex:",tstart-tmutex
                t=t+dt

        raw_input('Press <ENTER> to execute trajectory.')

