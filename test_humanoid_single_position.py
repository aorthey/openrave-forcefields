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
from cbirrtpy import *

from deformation_naive import *
from deformation_potentials import *
from deformation_stretchpull import *
from trajectory_bspline import *
import numpy as np
from motion_planner_geometrical import MotionPlannerGeometrical
from motion_planner_kinodynamic import MotionPlannerKinodynamic
from util_humanoid import *

#sys.path.append(os.environ["MPP_PATH"]+"mpp-robot/mpp")
#sys.path.append(os.environ["OPENRAVE_WPI_PATH"]+"ros_msgs")
#from SurfaceMeshes import *


def waitrobot(robot):
    """busy wait for robot completion"""
    while not robot.GetController().IsDone():
        time.sleep(0.1)

if __name__ == "__main__":

        env = EnvironmentHumanoidContactWorld()
        env.DrawAxes()
        surfaces = env.GetSurfaces()

        print surfaces.shape
        for i in arange(surfaces.shape[0]):
                et = surfaces[i,4,0]
                eo = surfaces[i,5,0]
                print "Surface",i,
                if et <= 0.01:
                        print "no extension",
                if eo <= 0.01:
                        print "no extension",
                print

        from surface_module import *
        M = surfaces.shape[0]

        S = SurfaceModule(surfaces)
        robot = env.GetRobot()
        com = robot.GetCenterOfMass()
        k = S.GetRelevantSurfaces(com)
        ##computeRelevantSurfaces()
        S.SampleSurface(1000,4,env)

        #env.DisplayForces()
        time.sleep(0.5)
        robot = env.GetRobot()

        raw_input('Press <ENTER> to execute trajectory.')


