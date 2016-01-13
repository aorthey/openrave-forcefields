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

from openravepy import *
if __name__ == "__main__":
    #env = EnvironmentTheStream()
    #env = EnvironmentTheCounterStream()
    #env = EnvironmentTheRay()
    #env = EnvironmentTheHideout()
    env = EnvironmentTripleStream()
    #env = EnvironmentCrossroadStream()
    robot = env.GetRobot()

    ### NEEDS FFMPEG -> install from source and rebuild openrave :-(

    #env.VideoRecordStart("hideout.mpeg")
    #env.DisplayForces()
    #time.sleep(1)
    #env.env.GetViewer().SendCommand('SetFiguresInCamera 1') # also shows the figures
    #I = env.env.GetViewer().GetCameraImage(640,480, env.env.GetViewer().GetCameraTransform(),[640,640,320,240])
    #scipy.misc.imsave('openrave.jpg',I)

    recorder = RaveCreateModule(env.env,'viewerrecorder')
    print recorder

    raw_input('Enter any key to quit. ')
