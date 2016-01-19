#!/usr/bin/env python
import time
import openravepy
import scipy
import numpy as np
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

###
if __name__ == "__main__":

        env = EnvironmentTheRay()
        robot = env.GetRobot()
        t = 0
        #env.DisplayForces()
        #time.sleep(0.2)

        [xi,yi]=env.RobotGetInitialPosition()

        print robot.GetDOF()
        print "###############################################################"
        print xi,yi
        print "###############################################################"
        #controller = RaveCreateController(env,'MyController controller arguments here')
        #robot.SetController(controller,range(robot.GetDOF()),controltransform=1)

        #######################################################################
        ###PARAMS TEST
        #######################################################################
        params = Planner.PlannerParameters()
        params.SetRobotActiveJoints(robot)
        goal=env.RobotGetGoalPosition()
        params.SetGoalConfig(goal)
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

        existing_planners = ['OMPL_KPIECE1']

        print existing_planners


        #######################################################################
        ### PLANNING
        #######################################################################
        for P in existing_planners:
                planner=RaveCreatePlanner(env.env,P)

                if planner is None:
                        print "###############################################################"
                        print "PLANNER",P,"not implemented"
                        print "###############################################################"
                        continue


                print "###############################################################"
                print "EXECUTING PLANNER:",P
                print "###############################################################"
                #planner.SendCommand('GetParameters')
                planner.InitPlan(env.robot, params)

                traj = RaveCreateTrajectory(env.env,'')
                planner.PlanPath(traj)

                #######################################################################
                #basemanip=interfaces.BaseManipulation(robot)
                N = traj.GetNumWaypoints()
                W=[]
                for i in range(0,N):
                        #print traj.GetWaypoint(i)
                        w = array((traj.GetWaypoint(i)[0],traj.GetWaypoint(i)[1],0.2))
                        W.append((w))
                        #traj=basemanip.MoveToHandPosition(matrices=[w],execute=True,outputtrajobj=True)

                
                #print N,"waypoints",W
                with env.env:
                        env.handles.append(env.env.drawlinestrip(points=array(W),
                                                   linewidth=5.0,
                                                   colors=array(((0.2,0.8,0.2)))))

                #######################################################################
                basemanip = interfaces.BaseManipulation(robot) # create the interface
                res = basemanip.MoveActiveJoints(goal = [goal[0],goal[1]], maxiter=3000,steplength=0.01,execute=True)

                #RaveSetDebugLevel(DebugLevel.Debug) # set output level to debug
                openravepy.RaveLogInfo("Waiting for controller to finish")
                robot.WaitForController(0)
                robot.GetController().Reset()
                time.sleep(2)
                del env.handles[-1]
                
        raw_input('Enter any key to quit. ')

