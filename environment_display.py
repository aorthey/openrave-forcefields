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
    env = EnvironmentTheStream()
    #env = EnvironmentTheCounterStream()
    cells  = env.GetCells()
    forces = env.GetForces()
    robot = env.GetRobot()
    env.DisplayForces()

    handles=[]
    Fx=0.0
    Fy=0.0
    rx=0.0
    ry=0.0

    robot.WaitForController(0)
    raw_input('Enter any key to quit. ')
