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

        #while scale<0.8:

        #        COM_offset = np.zeros((3,M))
        #        COM_offset += COM_linear
        #        COM_offset[0,:]+=0.00*exp(-(T-M/2)**2/(2*C*C))
        #        COM_offset[1,:]+=-scale*exp(-(T-M/2)**2/(2*C*C))
        #        COM_offset[2,:]+=-scale*0.1*exp(-(T-M/2)**2/(2*C*C))

        #        env.DrawAxes()

        #        tmp_handle=[]
        #        with env.env:
        #                h=env.env.drawlinestrip(points=COM_offset.T,linewidth=6,colors=np.array((0,1,0)))
        #                tmp_handle.append(h)
        #        [COM_zig_zag, footpos, dfootpos] = COM_compute_zig_zag_motion(COM_offset, env)
        #        scale+=0.05
        #        time.sleep(0.1)

        #raw_input('Press <ENTER> to start.')

        #traj = Trajectory.from_ravetraj(rave_path)
        #traj.info()
        #traj.draw(env)
        ##traj.PlotParametrization(env)
        #traj.draw_delete()
        #td = DeformationStretchPull(traj, env)
        #td.deform(N_iter=100)

        #time.sleep(0.5)

        #Mcut = 60
        #COM_zig_zag = COM_zig_zag[:,0:Mcut]
        #footpos = footpos[0:Mcut,:]
        #dfootpos = dfootpos[0:Mcut,:]
        #q_original = q_original[:,0:Mcut]

        #C = np.around(COM_zig_zag.T,decimals=2)
        #F = np.around(footpos,decimals=2)
        #Z = np.arange(0,Mcut)
        #Z1 = np.hstack((C,F))

        #print Z1[0:40]
        #print Z1[40:Mcut]

        #sys.exit(0)

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

        ##forcetorquemap = np.zeros((Nl,6))
        #forcetorquemap = {x: np.zeros(6) for x in np.arange(Nl)}
        ### create forcetorquemap dictionary
        #link = robot.GetLink('torso')
        #iTorso = link.GetIndex()
        #Ftorso = np.array((0,0.0,0,0,0,0))
        #forcetorquemap[iTorso]=Ftorso

        with env.env:
                rave.planningutils.SmoothActiveDOFTrajectory(traj,robot)
                rave.planningutils.RetimeActiveDOFTrajectory(traj,robot,hastimestamps=False,maxvelmult=0.75)

        Td = traj.GetDuration()

        active_dofs = robot.GetActiveConfigurationSpecification()

        dt = 0.1
        t = 0.0
        while t < Td:
                qt = traj.Sample(t,active_dofs)
                t=t+dt
        sys.exit(0)
        dtsleep = 0.1

        #raw_input('Press <ENTER> to execute trajectory.')

        with env.env:
                #env.env.GetPhysicsEngine().SetGravity([0,0,-0.981])
                env.env.GetPhysicsEngine().SetGravity([0,0,0])
                robot.GetLinks()[0].SetStatic(True)
                #F = np.array((-1000,0,0))
                #env.AddForceAtTorso(robot, F)
                #robot.SetController(RaveCreateController(env.env,'odevelocity'),range(robot.GetDOF()),0)
                env.env.StopSimulation()
                env.env.StartSimulation(timestep=0.001)
                #qd = np.zeros(robot.GetDOF())
                #robot.GetController().SendCommand('setvelocity '+' '.join(str(qdi) for qdi in qd))

        robot.GetController().SetPath(traj)
        robot.WaitForController(0)
        robot.GetController().Reset()
        raw_input('Enter any key to quit. ')


        #robot.GetController().SetPath(rave_path)
        #robot.GetController().SendCommand('setvelocity '+' '.join(str(f) for f in velocities))
        dtsleep = 0.1

        def printProgress(i, d, dall):
                t=d/dall
                T = int(t*10)
                print "[",i,"]",
                for i in range(0,T):
                        ".",
                print

        i = 3
        q_old = np.zeros(robot.GetDOF())
        q_next = np.zeros(robot.GetDOF())
        q_last = np.zeros(robot.GetDOF())
        qd = np.zeros(robot.GetDOF())

        with env.env:
                robot.SetActiveDOFValues(q_gik[:,i])
                q_cur = robot.GetDOFValues()
                env.env.StopSimulation()
                env.env.StartSimulation(timestep=0.001)

        while i<M:
                ictr = 0 
                while True:
                        with env.env:
                                q_cur = robot.GetDOFValues()

                                robot.SetActiveDOFValues(q_gik[:,i+1])
                                q_next = robot.GetDOFValues()

                                robot.SetActiveDOFValues(q_gik[:,i])
                                q_last = robot.GetDOFValues()

                                robot.SetActiveDOFValues(q_cur)

                                qd = 0.001*(q_next-q_last)/dt
                                C = robot.GetController()
                                dq = np.linalg.norm(qd)

                                #qd = np.zeros(robot.GetDOF())
                                #C.SendCommand('setvelocity '+' '.join(str(qdi) for qdi in qd))

                                ### draw line from q_last to q_next, check q_cur
                                ### position on this line
                                dstart = np.linalg.norm(q_cur-q_last)
                                dtarg = np.linalg.norm(q_cur-q_next)
                                dline = np.linalg.norm(q_next-q_last)

                                if dstart > dline:
                                        ## overshoot
                                        print "overshooting"
                                        qd = np.zeros(robot.GetDOF())
                                        robot.GetController().SendCommand('setvelocity '+' '.join(str(qdi) for qdi in qd))

                                #tm,tc,tg = robot.ComputeInverseDynamics([],None,returncomponents=True)
                                #qdd = (tm+tc+tg)
                                #print qdd
                                #qd = dt*qdd + (q_new-q_old)/dt
                                #printProgress(i,d,dall)

                                print ictr,"vel:",dq,"overall dist:",dline," dist from start:",dstart,"dist to targ:",dtarg
                        time.sleep(dt)
                        ictr+=1
                        if ictr > 20:
                                sys.exit(0)
                i=i+1
                sys.exit(0)
