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
from environment_force_humanoid_contact_world_surfboard import *
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

def EvalSurfboardPath(t):

        initTheta = -pi/8
        goalTheta = pi/8
        initX = np.array((-1,-1,0.5))
        goalX = np.array((-1,0.5,0.5))

        Theta = t*goalTheta + (1-t)*initTheta
        X = t*goalX + (1-t)*initX

        T = np.eye(4)
        T[0:3,0:3] = np.dot(Rz(Theta),Rx(Theta),Ry(Theta))
        T[0:3,3] = X
        return T


def waitrobot(robot):
    """busy wait for robot completion"""
    while not robot.GetController().IsDone():
        time.sleep(0.1)

if __name__ == "__main__":

        np.random.seed(2)
        env = EnvironmentHumanoidContactWorldSurfboard()
        env.DrawAxes()
        #env.MakeRobotTransparent(alpha=1.0)
        time.sleep(0.5)
        env.DisplayForces()

        #######################################################################
        gik = GIKInterface(env)
        robot = env.GetRobot()
        #######################################################################

        left_leg_surface = 4
        right_leg_surface = 4


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

        left_leg_relative_surface_tf = S.GetContactRelativeToSurfaceTransform(left_leg_tf, left_leg_surface) 
        right_leg_relative_surface_tf = S.GetContactRelativeToSurfaceTransform(right_leg_tf, right_leg_surface) 

        W = env.env.GetKinBody('world')
        B = W.GetLink('surfboard')
        #T = B.GetTransform()

        cog_handler = []

        Mwaypoints = 20
        MsamplesPerStance = 2000

        tvec = linspace(0,1,Mwaypoints)

        ictr=0
        for t in tvec:
                fh = open("trajectories/surfboard_com_"+str(ictr)+".txt", 'w')
                fh2 = open("trajectories/surfboard_com_all_"+str(ictr)+".txt", 'w')
                ictr+=1
                t1 = time.time()
                for m in range(0,MsamplesPerStance):
                        with env.env:
                                T = EvalSurfboardPath(t)
                                B.SetTransform(T)

                                surfaces = env.GetSurfaces()
                                S = SurfaceModule(surfaces)
                                #######################################################################

                                #######################################################
                                ### draw
                                #######################################################

                                #right_leg_tf = S.AdjustToNewSurfaceTransform( 0.02, right_leg_relative_surface_tf, right_leg_surface, env)
                                #left_leg_tf = S.AdjustToNewSurfaceTransform( 0.02, left_leg_relative_surface_tf, left_leg_surface, env)

                                right_leg_tf = np.eye(4)
                                left_leg_tf = np.eye(4)
                                sy =  np.dot(T[0:3,0:3],ey)
                                sz =  np.dot(T[0:3,0:3],ez)
                                dl = 0.25
                                dz = 0.02

                                right_leg_tf[0:3,0:3] = T[0:3,0:3]
                                right_leg_tf[0:3,3]= T[0:3,3]-dl*sy+dz*sz
                                left_leg_tf[0:3,0:3] = T[0:3,0:3]
                                left_leg_tf[0:3,3]= T[0:3,3]+dl*sy+dz*sz


                                X = T[0:3,3]
                                xmin = X[0]-0.4
                                xmax = X[0]+0.4

                                ymin = X[1]-0.4
                                ymax = X[1]+0.4

                                zmin = X[2]-0.0
                                zmax = X[2]+1.5

                                cog_in  = np.zeros(3)
                                cog_in[0] = np.random.uniform(xmin,xmax)
                                cog_in[1] = np.random.uniform(ymin,ymax)
                                cog_in[2] = np.random.uniform(zmin,zmax)

                                res = gik.giwc(robot, cog_in, left_leg_tf, right_leg_tf, None, None)
                                if res:
                                        cog = gik.cog
                                        fh.write("%f %f %f\n" % (cog[0],cog[1],cog[2]))
                                if gik.cogIsValid:
                                        cog = cog_in
                                        fh2.write("%f %f %f\n" % (cog[0],cog[1],cog[2]))
                                #######################################################
                                ### random noise on environment
                                #######################################################
                                print "#############################################################"
                                print "WP:",ictr,"/",Mwaypoints,"sample:",m,"/",MsamplesPerStance
                                print "#############################################################"

                t2 = time.time()
                tall = t2-t1
                print "sampling",MsamplesPerStance,"samples, time:",tall
                tps = float(tall)/float(MsamplesPerStance)
                print "time per sample:",tps
                print "estimated time:",tps*Mwaypoints*MsamplesPerStance
                fh.close()
                fh2.close()


        robot.WaitForController(0)
