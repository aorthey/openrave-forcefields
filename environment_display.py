#!/usr/bin/env python
import time
import openravepy
from math import *
from environment_the_stream import *
from environment_the_counterstream import *

if not __openravepy_build_doc__:
    from openravepy import *
    from numpy import *

if __name__ == "__main__":
    ###########################################################################
    ## environment setter
    #env = Environment()
    #env.SetViewer('qtcoin')
    #env.Reset()

    ############################################################################
    ### loader
    ##xmlenv='../../src/data/point_in_forcefield.env.xml'
    #xmlenv='environments/the_stream.env.xml'
    #xmlrobot='robots/pointrobot.robot.xml'
    #env.Add(env.ReadRobotXMLFile(xmlrobot))
    #env.Load(xmlenv)
    ############################################################################
    ### conventional variables
    #I = eye(4)
    #openravepy.misc.DrawAxes(env,I)
    ############################################################################
    ## create force field
    #physics = RaveCreatePhysicsEngine(env,'ode')
    #physics.SetGravity(array((0,0,-9.81)))
    #env.SetPhysicsEngine(physics)
    ###########################################################################

    #env = EnvironmentTheStream()
    env = EnvironmentTheCounterStream()
    cells  = env.GetCells()
    forces = env.GetForces()
    robot = env.GetRobot()
    env.DisplayForces()

    handles=[]
    Fx=0.0
    Fy=0.0
    rx=0.0
    ry=0.0

    raw_input('Enter any key to quit. ')
