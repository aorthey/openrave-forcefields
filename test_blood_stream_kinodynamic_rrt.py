#!/usr/bin/env python
import time
import scipy
import sys
import numpy as np
import openravepy
from openravepy import *
from math import *
from environment_force_the_stream import *
from environment_force_the_counterstream import *
from environment_force_the_ray import *
from environment_force_blood_stream import *
from environment_force_blood_stream2 import *
from environment_force_radial import *
from environment_periodic_force_the_hideout import *
from environment_periodic_force_triple_stream import *
from environment_periodic_force_crossroad_stream import *

from deformation_naive import *
from deformation_potentials import *
from deformation_stretchpull import *
from deformation_reachableset import *
from trajectory_bspline import *
import numpy as np
from motion_planner_geometrical import MotionPlannerGeometrical
from motion_planner_kinodynamic import MotionPlannerKinodynamic
#import statsmodels.api as sm

if __name__ == "__main__":

        #######################################################################
        #env = EnvironmentBloodStream2()
        env = EnvironmentRadial()
        #env = EnvironmentTheCounterStream()
        #env = EnvironmentTheStream()
        #######################################################################

        robot = env.GetRobot()
        env.MakeRobotInvisible()
        env.DisplayForces()
        time.sleep(0.5)

        #planner = MotionPlannerGeometrical(robot, env)
        fh = open('trajectories/time_stream_krrt_6.txt','w')

        for i in range(0,10):
                t1 = time.time()
                try:
                        planner = MotionPlannerKinodynamic(robot, env)
                        print "Starting next planning instance"
                        rave_path = planner.GetPath()
                        if rave_path is not None:
                                traj = Trajectory.from_ravetraj(rave_path)
                                traj.save('trajectories/bloodstream_kinodynamic_6'+str(i))
                except Exception as e:
                        pass
                t2 = time.time()
                tall = t2-t1
                fh.write('%d %f\n'%(i,tall))
                print "time:",tall
                time.sleep(5)
        fh.close()
