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
from trajectory import *
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

def VisualizeCOMSet(q, env, tricolor=np.array((1,0,1,1))):
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

        tri = Triangulation(x, y, triangles=hull.simplices)

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
        h = env.env.drawtrimesh(points=q,indices=hull.simplices,colors=tricolor)
        handler.append(h)

        #h = env.env.plot3(points=q,
        #        pointsize=0.01,
        #        colors=array((1,0,1,1)),
        #        drawstyle=1)
        #handler.append(h)
        #h = env.env.plot3(points=q[hull.vertices,:],
        #        pointsize=0.01,
        #        colors=array((1,0,1,1)),
        #        drawstyle=1)
        #handler.append(h)
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
        #env.DisplayForces()

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

        DEBUG=True

        np.random.seed(0)
        cog_in = np.loadtxt("trajectories/surfboard_com.txt",dtype='float')
        q_in = np.loadtxt("trajectories/surfboard_q.txt",dtype='float')

        h = env.env.drawlinestrip(points=cog_in,
                        linewidth=5,
                        colors=array((1,0,1,1)))
        cog_handler.append(h)

        from trajectory import Trajectory

        ctraj = Trajectory(cog_in.T)
        qtraj = Trajectory(q_in.T)

        with env.env:
                handler = []
                for i in range(0,20):
                        qcom = np.loadtxt("trajectories/surfboard_com_"+str(i)+".txt",dtype='float')
                        qcom2 = np.loadtxt("trajectories/surfboard_com_all_"+str(i)+".txt",dtype='float')
                        h = VisualizeCOMSet(qcom, env, np.array((0.0,1.0,0.0,0.1)))
                        handler.append(h)
                        h = VisualizeCOMSet(qcom2, env, np.array((1.0,0.0,1.0,0.02)))
                        handler.append(h)

        Mwaypoints=200
        tvec = linspace(0,1,Mwaypoints)

        for t in tvec:
                t1 = time.time()
                with env.env:
                        [cog,dcog] =ctraj.evaluate_at(t)
                        [q,dq] =qtraj.evaluate_at(t)
                        T = EvalSurfboardPath(t)
                        B.SetTransform(T)

                        surfaces = env.GetSurfaces()
                        S = SurfaceModule(surfaces)
                        #######################################################################

                        #######################################################
                        ### draw
                        #######################################################

                        right_leg_tf = np.eye(4)
                        left_leg_tf = np.eye(4)
                        sy =  np.dot(T[0:3,0:3],ey)
                        sz =  np.dot(T[0:3,0:3],ez)
                        dl = 0.25
                        dz = 0.03

                        right_leg_tf[0:3,0:3] = T[0:3,0:3]
                        right_leg_tf[0:3,3]= T[0:3,3]-dl*sy+dz*sz
                        left_leg_tf[0:3,0:3] = T[0:3,0:3]
                        left_leg_tf[0:3,3]= T[0:3,3]+dl*sy+dz*sz

                        #print cog
                        #res = gik.giwc(robot, cog, left_leg_tf, right_leg_tf, None, None)
                        robot.SetActiveDOFValues(q)

                        if env.static_handles:
                               env.static_handles = None
                               env.static_handles = []
                        env.DrawFootContactPatchFromTransform(left_leg_tf, cpatch=green)
                        env.DrawFootContactPatchFromTransform(right_leg_tf, cpatch=green)
                        S.DrawCoordinateFrames(env)

                env.env.StepSimulation(0.002)
                time.sleep(0.2)
                if t<=0:
                        raw_input('Press <ENTER> to execute com traj.')
                t2 = time.time()
                tall = t2-t1
                print "time:",tall


        robot.WaitForController(0)
