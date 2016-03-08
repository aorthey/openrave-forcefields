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
from environment_periodic_force_the_hideout import *
from environment_periodic_force_triple_stream import *
from environment_periodic_force_crossroad_stream import *

from deformation_naive import *
from deformation_potentials import *
from trajectory_bspline import *
import numpy as np
import statsmodels.api as sm

if __name__ == "__main__":

        env = EnvironmentTheRay()
        #env = EnvironmentTheCounterStream()
        #env = EnvironmentTheStream()

        robot = env.GetRobot()
        env.DisplayForces()
        time.sleep(0.2)

        print robot.GetDOF()
        print "###############################################################"
        #controller = RaveCreateController(env,'MyController controller arguments here')
        #robot.SetController(controller,range(robot.GetDOF()),controltransform=1)

        #######################################################################
        ###PARAMS TEST
        #######################################################################
        params = Planner.PlannerParameters()
        params.SetRobotActiveJoints(robot)

        init=env.RobotGetInitialPosition()
        goal=env.RobotGetGoalPosition()

        params.SetInitialConfig(init)
        params.SetGoalConfig(goal)

        print env.env.GetForces()

        #######################################################################

        existing_planners=[
                'birrt',
                'OMPL_BKPIECE1',
                'OMPL_EST',
                'OMPL_KPIECE1',
                'OMPL_LazyRRT',
                'OMPL_LBKPIECE1',
                'OMPL_PDST',
                'OMPL_PRM',
                'OMPL_LazyPRM',
                'OMPL_PRMstar',
                'OMPL_pSBL',
                'OMPL_RRT',
                'OMPL_RRTConnect',
                'OMPL_RRTstar',
                'OMPL_SBL',
                'OMPL_SPARS',
                'OMPL_SPARStwo',
                'OMPL_TRRT']

        ## not working:
        #'OMPL_BITstar',
        #'OMPL_FMT',
        #'OMPL_pRRT',

        print existing_planners
        P = 'birrt'
        #P = 'kinodynamicrrt'
        #P = 'basicrrt'

        planner=RaveCreatePlanner(env.env,P)
        #print planner.SendCommand('GetParameters')
        #params.SetExtraParameters('<range>0.01</range>')

        if planner is None:
                print "###############################################################"
                print "PLANNER",P,"not implemented"
                print "###############################################################"
                sys.exit(0)

        print "###############################################################"
        print "EXECUTING PLANNER:",P
        print "###############################################################"

        #######################################################################
        planner.InitPlan(env.robot, params)

        rave_traj = RaveCreateTrajectory(env.env,'')
        result = planner.PlanPath(rave_traj)
        assert result == PlannerStatus.HasSolution


        #from util import draw_waypoints, draw_ravetraj
        #handle = draw_ravetraj(rave_traj, env)
        from trajectory_polynomial import *

        #traj = TrajectoryPolynomial.from_ravetraj(rave_traj)
        traj = TrajectoryBSpline.from_ravetraj(rave_traj)
        traj.info()
        #traj.plot_speed_profile(env)
        #traj.draw(env)

        #td = DeformationNaive(traj, env)
        td = DeformationPotentials(traj, env)
        for i in range(10):
                td.deform()
                td.draw_deformation()

        #td.traj_current.plot_speed_profile(env)

        #ravetraj = td.traj_current.to_ravetraj()
        #result = planningutils.RetimeTrajectory(ravetraj,False,0.15)

        #td.deform()
        #td.draw_deformation()

        #raw_input('Press <ENTER> to execute trajectory.')
        #RaveSetDebugLevel(DebugLevel.Debug) # set output level to debug
        #openravepy.RaveLogInfo("Waiting for controller to finish")
        #robot.GetController().SetPath(traj)
        #robot.WaitForController(0)
        #robot.GetController().Reset()
                
        raw_input('Enter any key to quit. ')

