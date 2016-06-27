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
from trajectory import *
import numpy as np
from motion_planner_geometrical import MotionPlannerGeometrical
from motion_planner_kinodynamic import MotionPlannerKinodynamic
#import statsmodels.api as sm
np.set_printoptions(precision=2)

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
        ##planner = MotionPlannerKinodynamic(robot, env)

        #rave_path = planner.GetPath()

        #traj = Trajectory.from_ravetraj(rave_path)
        traj = Trajectory.from_file('deform1')
        traj.info()
        traj.draw(env)
        xml = env.GetName()
        traj.save('trajectories/'+xml)

        #raw_input('Press <ENTER> to deform.')

        time.sleep(1)
        traj.draw_delete()

        td = DeformationReachableSet(traj, env)
        deform_success = td.deform(N_iter=100)

        td.traj_deformed.save('trajectories/'+xml+'_deformed')

        if deform_success:
                td.traj_deformed.PlotParametrization(env)
                td.traj_deformed.execute(env, robot, tsleep=0.005,
                                stepping=False)
                #td.execute(robot, tsleep=0.003, stepping=True)
