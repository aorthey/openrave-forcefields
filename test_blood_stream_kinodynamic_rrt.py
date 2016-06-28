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
        env = EnvironmentBloodStream()
        #env = EnvironmentTheCounterStream()
        #env = EnvironmentTheStream()
        #######################################################################

        robot = env.GetRobot()
        env.MakeRobotInvisible()
        env.DisplayForces()
        time.sleep(0.5)

        #planner = MotionPlannerGeometrical(robot, env)
        fh = open('trajectories/time_stream_krrt.txt','w')

        for i in range(0,3):
                t1 = time.time()
                planner = MotionPlannerKinodynamic(robot, env)
                t2 = time.time()
                tall = t2-t1
                fh.write('%d %f\n'%(i,tall))
                print "time:",tall
                rave_path = planner.GetPath()
                traj = Trajectory.from_ravetraj(rave_path)
                traj.info()
                #traj.draw(env)
                traj.save('trajectories/bloodstream_kinodynamic'+str(i))
        fh.close()


        raw_input('Press <ENTER> to start deforming.')
        time.sleep(1)
        traj.draw_delete()

        td = DeformationReachableSet(traj, env)
        deform_success = td.deform(N_iter=100)

        td.traj_deformed.save('trajectories/bloodstream_kinodynamic_deformed')

        if deform_success:
                td.traj_deformed.PlotParametrization(env)
                td.traj_deformed.execute(env, robot, tsleep=0.003, stepping=False)
                #td.execute(robot, tsleep=0.003, stepping=True)
