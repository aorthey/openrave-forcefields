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
                env.env.GetPhysicsEngine().SetGravity([0,0,-9.81])
                rave.planningutils.SmoothActiveDOFTrajectory(traj,robot)
                rave.planningutils.RetimeActiveDOFTrajectory(traj,robot,hastimestamps=False,maxvelmult=0.75)

        Td = traj.GetDuration()

        active_dofs = robot.GetActiveConfigurationSpecification()
        forcetorquemap = {x: np.zeros(6) for x in np.arange(Nl)}
        ### create forcetorquemap dictionary
        link = robot.GetLink('torso')
        iTorso = link.GetIndex()
        Ftorso = np.array((0,0.0,0,0,0,0))

        for geom in link.GetGeometries():
                c = np.array((1,0,0))
                geom.SetAmbientColor(c) 
                geom.SetDiffuseColor(c) 

        forcetorquemap[iTorso]=Ftorso

        dt = 0.1
        t = 0.0
        qold = traj.Sample(0,active_dofs)
        qext = np.zeros((robot.GetDOF()))
        zlvec = []
        zrvec = []
        ictr=0
        while t < Td:
                t0 = time.time()
                qt = traj.Sample(t,active_dofs)

                tmutex = time.time()
                with env.env:
                        #### add contact forces
                        tstart = time.time()
                        rlink = robot.GetManipulator('r_leg')
                        llink = robot.GetManipulator('l_leg')
                        zr = rlink.GetTransform()[2,3]
                        zl = llink.GetTransform()[2,3]
                        zlvec.append(zl)
                        zrvec.append(zr)
                        #iRlink = rlink.GetIndex()
                        #iLlink = llink.GetIndex()

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
                                Ffoot = np.array((0,0,-0,0,0,0))
                                forcetorquemap[iLlink]=Ffoot
                                env.DrawArrow(llink.GetGlobalCOM(), -Ffoot[0:3])
                        else:
                                forcetorquemap[iLlink]=np.zeros(6)
                                env.DrawArrow(llink.GetGlobalCOM(), np.zeros(3))
                        #        link = robot.GetLink('l_foot')
                        #        for geom in link.GetGeometries():
                        #                c = np.array((1,1,1))
                        #                geom.SetAmbientColor(c) 
                        #                geom.SetDiffuseColor(c) 
                        #if zr < 0.01:
                        #        forcetorquemap[iRlink]=-Ftorso
                        #else:
                        #        forcetorquemap[iRlink]=np.zeros(6)


                        tdyn = time.time()
                        robot.SetActiveDOFValues(qold)
                        tm,tc,tg = robot.ComputeInverseDynamics([],forcetorquemap,returncomponents=True)

                        qdd = (tm+tc+tg)
                        qext = qdd*dt*dt*0.5 + qext

                        robot.SetDOFValues(qext)
                        qext_active = robot.GetActiveDOFValues()

                        q = qt + qext_active
                        robot.SetActiveDOFValues(q)
                        qold = q
                        ictr +=1
                        tend = time.time()
                        #print "time:",tend-tstart,"map:",tdyn-tstart,"invdyn:",tend-tdyn
                        #if ictr > 30:
                        #        plot(zlvec,'-b')
                        #        plot(zrvec,'-r')
                        #        plt.show()
                tend = time.time()
                time.sleep(dt)
                tsleep = time.time()
                print "time:",tsleep-t0,"sleep:",tsleep-tend,"invdyn:",tend-tdyn,"tsample:",tstart-t0,"tmutex:",tstart-tmutex
                t=t+dt
        sys.exit(0)
        dtsleep = 0.1

        raw_input('Press <ENTER> to execute trajectory.')

