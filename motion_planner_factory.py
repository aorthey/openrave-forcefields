import abc
import sys
import time
import numpy as np
import openravepy
from openravepy import Planner, RaveCreatePlanner, RaveCreateTrajectory, PlannerStatus

class MotionPlanner():
        __metaclass__ = abc.ABCMeta
        env_ = []
        robot_ = []
        existing_planners=[]

        def __init__(self, robot, env):
                self.env_ = env
                self.robot_ = robot
                self.existing_planners=[
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


        @abc.abstractmethod 
        def GetPath(self):
                pass

        def PlanPath(self, planner_name):

                params = Planner.PlannerParameters()
                params.SetRobotActiveJoints(self.robot_)

                init=self.env_.RobotGetInitialPosition()
                goal=self.env_.RobotGetGoalPosition()

                params.SetInitialConfig(init)
                params.SetGoalConfig(goal)
                planner=RaveCreatePlanner(self.env_.env,planner_name)

                if planner is None:
                        print "###############################################################"
                        print "PLANNER",planner_name,"not implemented"
                        print "###############################################################"
                        sys.exit(0)

                #######################################################################
                planner.InitPlan(self.robot_, params)

                rave_traj = RaveCreateTrajectory(self.env_.env,'')

                t1=time.time()
                result = planner.PlanPath(rave_traj)
                t2=time.time()
                if result != PlannerStatus.HasSolution:
                        print "Could not find geometrical path"
                        print "Planner:",planner_name
                        print "Status :",result
                        sys.exit(0)
                print "Planner",planner_name," success | time:",t2-t1
                return rave_traj

