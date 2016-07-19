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
from environment_force_radial import *
from environment_force_blood_stream2 import *
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
        env = EnvironmentBloodStream2()
        #env = EnvironmentBloodStream()
        #env = EnvironmentRadial()
        #env = EnvironmentTheCounterStream()
        #env = EnvironmentTheStream()
        #######################################################################

        robot = env.GetRobot()
        env.MakeRobotInvisible()
        env.DisplayForces()
        time.sleep(0.5)

        #planner = MotionPlannerGeometrical(robot, env)
        fh = open('trajectories/time_stream_rrtdeform10.txt','w')

        for i in range(0,10):
                t1 = time.time()
                planner = MotionPlannerGeometrical(robot, env)
                rave_path = planner.GetPath()
                t3 = time.time()
                traj = Trajectory.from_ravetraj(rave_path)
                #traj.info()
                td = DeformationReachableSet(traj, env)
                deform_success = td.deform(N_iter=100)
                t2 = time.time()
                tall = t2-t1

                if deform_success:
                        #traj.draw(env)
                        fh.write('%d %f\n'%(i,tall))
                        print "time rrt:",t3-t1
                        print "time:",tall
                        td.traj_deformed.save('trajectories/bloodstream_geometric10'+str(i))
                else:
                        fh.write('%d %f\n'%(i,-1))
        fh.close()
