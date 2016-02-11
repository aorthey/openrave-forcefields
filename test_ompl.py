#P!/usr/bin/env python
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

from trajectory_analyzer import *
#from trajectory_speed_profiler import *
from trajectory_deformation_naive import *

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
        P = 'kinodynamicrrt'
        P = 'birrt'
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

        traj = RaveCreateTrajectory(env.env,'')
        result = planner.PlanPath(traj)
        assert result == PlannerStatus.HasSolution

        #result = planningutils.RetimeTrajectory(traj,False,0.15)
        #assert result == PlannerStatus.HasSolution

        ########################################################################
        #print "###############################################################"
        #print "SIMPLIFYING PLAN"
        #print "###############################################################"
        ########################################################################
        
        #simplifier.InitPlan(env.robot, Planner.PlannerParameters())
        #result = simplifier.PlanPath(traj)
        #assert result == PlannerStatus.HasSolution

        #######################################################################
        #basemanip=interfaces.BaseManipulation(robot)

        N = traj.GetNumWaypoints()
        #W=[]
        #for i in range(0,N):
        #        w = array((traj.GetWaypoint(i)[0],traj.GetWaypoint(i)[1],traj.GetWaypoint(i)[2]))
        #        W.append((w))
        
        print "###############################################################"
        print N,"waypoints"
        print traj.GetWaypoint(N-1)
        print "###############################################################"
        #jwith env.env:
        #j        env.handles.append(env.env.drawlinestrip(points=array(W),
        #j                                   linewidth=10.0,
        #j                                   colors=array(((0.2,0.8,0.2)))))

        #######################################################################
        #W = np.array((1,0,0),(1.5,0,0),(2.0,0,0),(2.5,0,0))
        #W = np.array(W)[:,0:2].T
        #print W
        
        td = TrajectoryDeformationNaive.from_ravetraj(traj, env)
        td.deform()
        td.draw_trajectory_original()
        td.draw_deformation()

        #ta = TrajectoryAnalyzer(W)
        #ta.draw(env)
        #time.sleep(0.1)

        #[told, tnew] = ta.deform_onestep(env)
        #ta.draw(env)
        #ta.draw_deformation(env, told,tnew)
        #tsp = TrajectorySpeedProfiler(robot)
        #traj = tsp.retime(traj)






        #raw_input('Press <ENTER> to execute trajectory.')
        #RaveSetDebugLevel(DebugLevel.Debug) # set output level to debug
        openravepy.RaveLogInfo("Waiting for controller to finish")
        #robot.GetController().SetPath(traj)
        #robot.WaitForController(0)
        #robot.GetController().Reset()
                
        raw_input('Enter any key to quit. ')

