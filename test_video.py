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
    env = EnvironmentTheRay()
    #env = EnvironmentTheHideout()
    #env = EnvironmentTripleStream()
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

    env.env.AddModule(recorder,'')
    codecs = recorder.SendCommand('GetCodecs') # linux only
    print codecs
#-vcodec libx264
    filename = 'openrave2.mpg'
    codec = 13
    recorder.SendCommand('Start 640 480 30 codec %d timing realtime filename %s\nviewer %s'%(codec,filename,env.env.GetViewer().GetName()))
    time.sleep(5)
    recorder.SendCommand('Stop') # stop the video
    env.env.Remove(recorder) # remove the recorder

    raw_input('Enter any key to quit. ')
