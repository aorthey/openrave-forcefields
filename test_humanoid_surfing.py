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

        initTheta = 0.0#-pi/8
        goalTheta = pi/8
        initX = np.array((-1,-1,0.5))
        goalX = np.array((-1,0.5,0.5))

        Theta = t*goalTheta + (1-t)*initTheta
        X = t*goalX + (1-t)*initX

        T = np.eye(4)
        T[0:3,0:3] = np.dot(Rz(Theta),Rx(Theta),Ry(Theta))
        T[0:3,3] = X
        return T

def VisualizeCOMSet(q, env):
        handler = []
        from scipy.spatial import ConvexHull
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        from matplotlib.tri import Triangulation, TriAnalyzer
        import pylab as plt
        qhull_options = 'QJ'
        hull = ConvexHull(q,qhull_options=qhull_options)    
        #q = q[hull.vertices,:]

        fig = plt.figure(facecolor='white')
        image = fig.gca(projection='3d')

        x,y,z=q.T
        print x,y

        tri = Triangulation(x, y, triangles=hull.simplices)
        print tri.triangles

        triangle_vertices = np.array([np.array([[x[T[0]], y[T[0]], z[T[0]]],
                [x[T[1]], y[T[1]], z[T[1]]],
                [x[T[2]], y[T[2]], z[T[2]]]]) for T in tri.triangles])

        #tri = Poly3DCollection(triangle_vertices)

        #ctri = np.array((0.5,0.5,1.0,0.2))
        #tri.set_color(ctri)
        #tri.set_edgecolor(np.array((0.5,0.5,0.5,0.5)))

        #tri.set_edgecolor('None')

        #image.scatter(x,y,z, 'ok', color=np.array((0,0,1.0,0.1)),s=20)
        #image.add_collection3d(tri)
        #plt.show()
        h = env.env.drawtrimesh(points=q,indices=hull.simplices,colors=array((0,1,0,0.2)))
        handler.append(h)

        h = env.env.plot3(points=q,
                pointsize=0.01,
                colors=array((1,0,1,1)),
                drawstyle=1)
        handler.append(h)
        h = env.env.plot3(points=q[hull.vertices,:],
                pointsize=0.01,
                colors=array((1,0,1,1)),
                drawstyle=1)
        handler.append(h)
        return handler


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

        ictr=0

        W = env.env.GetKinBody('world')
        B = W.GetLink('surfboard')
        #T = B.GetTransform()

        cog_handler = []

        Mwaypoints = 1
        MsamplesPerStance = 1000

        tvec = linspace(0,1,Mwaypoints)
        DEBUG=True

        qcom = np.loadtxt("trajectories/surfboard_com_test.txt",dtype='float')

        ictr=0
        with env.env:
                handler = VisualizeCOMSet(qcom, env)

        for t in tvec:
                ictr+=1
                t1 = time.time()
                for m in range(0,MsamplesPerStance):
                        with env.env:
                                #handler = VisualizeCOMSet(qcom, env)
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

                                cog = np.zeros(3)
                                cog[0] = np.random.uniform(xmin,xmax)
                                cog[1] = np.random.uniform(ymin,ymax)
                                cog[2] = np.random.uniform(zmin,zmax)

                                h = env.env.plot3(points=cog,
                                                pointsize=0.01,
                                                colors=array((1,0,0,1)),
                                                drawstyle=1)
                                cog_handler.append(h)

                                res = gik.giwc(robot, cog, left_leg_tf, right_leg_tf, None, None)

                                #res = False
                                if DEBUG:
                                        if res:
                                                print "valid com:",cog
                                                h = env.env.plot3(points=gik.cog,
                                                                pointsize=0.02,
                                                                colors=array((0,1,0,1)),
                                                                drawstyle=1)
                                                cog_handler.append(h)
                                        else:
                                                h = env.env.plot3(points=cog,
                                                                pointsize=0.01,
                                                                colors=array((1,0,0,1)),
                                                                drawstyle=1)
                                                cog_handler.append(h)

                                        if env.static_handles:
                                               env.static_handles = None
                                               env.static_handles = []
                                        env.DrawFootContactPatchFromTransform(left_leg_tf, cpatch=green)
                                        env.DrawFootContactPatchFromTransform(right_leg_tf, cpatch=green)
                                        S.DrawCoordinateFrames(env)
                                #else:
                                        #h = env.env.plot3(points=gik.cog,
                                        #                pointsize=0.01,
                                        #                colors=array((1,0,0,1)),
                                        #                drawstyle=1)
                                        #cog_handler.append(h)
                                #q = gik.fromContactTransform(robot, left_leg_tf,
                                                #right_leg_tf, None, None)

                                #gik.VisualizeFrictionCone(robot)

                                #torques = (numpy.random.rand(robot.GetDOF())-0.5)
                                #robot.SetJointTorques(torques,True)
                                #######################################################
                                ### random noise on environment
                                #######################################################

                        #if DEBUG:
                                #env.env.StepSimulation(0.002)
                                #time.sleep(0.5)
                t2 = time.time()
                tall = t2-t1
                print "sampling",MsamplesPerStance,"samples, time:",tall
                tps = float(tall)/float(MsamplesPerStance)
                print "time per sample:",tps
                print "estimated time:",tps*Mwaypoints*MsamplesPerStance


        robot.WaitForController(0)
