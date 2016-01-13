#!/usr/bin/env python
import time
import openravepy
import scipy
from math import *
from environment_force_the_stream import *
from environment_force_the_counterstream import *
from environment_force_the_ray import *
from environment_periodic_force_the_hideout import *
from environment_periodic_force_triple_stream import *
from environment_periodic_force_crossroad_stream import *

if not __openravepy_build_doc__:
    from openravepy import *
    from numpy import *

if __name__ == "__main__":
    #env = EnvironmentTheStream()
    #env = EnvironmentTheCounterStream()
    env = EnvironmentTheRay()
    #env = EnvironmentTheHideout()
    #env = EnvironmentTripleStream()
    #env = EnvironmentCrossroadStream()
    robot = env.GetRobot()
    #env.VideoRecordStart("hideout.mpeg")

    t = 0
    while True:
            env.DisplayForces()
            time.sleep(0.01)

    #env.VideoRecordStop()

    handles=[]
    Fx=0.0
    Fy=0.0
    rx=0.0
    ry=0.0

    robot.WaitForController(0)
    raw_input('Enter any key to quit. ')
