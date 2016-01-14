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

###
if __name__ == "__main__":

        env = EnvironmentTheRay()
        robot = env.GetRobot()

        t = 0
        env.DisplayForces()
        time.sleep(0.01)

        #robot.SetAffineTranslationLimits(-3,+3)
        #robot.SetActiveDOFs([],openravepy.DOFAffine.X|openravepy.DOFAffine.Y|openravepy.DOFAffine.Z,[0,0,0])
        #        goal = np.array(goal)
        #if goal.ndim == 2:
        #    #assume a 4x4 transformation matrix
        #    angle = openravepy.axisAngleFromRotationMatrix(goal)[-1]
        #    goal = np.array([goal[0,-1], goal[1,-1], angle])
        #    
        ## find the boundaries of the environment
        #envmin, envmax = utils.get_environment_limits(self.env, self.robot)        

        #with self.env:
        #    self.robot.SetAffineTranslationLimits(envmin,envmax)
        #    self.robot.SetAffineTranslationMaxVels([0.3,0.3,0.3])
        #    self.robot.SetAffineRotationAxisMaxVels(np.ones(4))
        #    self.robot.SetActiveDOFs([],
        #        openravepy.DOFAffine.X|openravepy.DOFAffine.Y|openravepy.DOFAffine.RotationAxis,
        #        [0,0,1])

        with env.env:
                #robot.SetActiveDOFs([],openravepy.DOFAffine.X|openravepy.DOFAffine.Y|openravepy.DOFAffine.Z,[0,0,1])
                robot.SetActiveDOFs([],
                                openravepy.DOFAffine.X|openravepy.DOFAffine.Y|openravepy.DOFAffine.RotationAxis,
                                [0,0,1])

        basemanip = interfaces.BaseManipulation(robot) # create the interface
        with env.env:
                goal = env.RobotGetGoalPosition()
                print goal
                res = basemanip.MoveActiveJoints(goal = [goal[0],goal[1],0.14], execute=True, outputtrajobj = True)

        N = res.GetNumWaypoints()
        W=[]
        for i in range(0,N):
                w = res.GetWaypoint(i)[0:3]
                W.append((w))
        
        with env.env:
                env.handles.append(env.env.drawlinestrip(points=array(W),
                                           linewidth=4.0,
                                           colors=array(((0.2,0.8,0.2)))))

        if res is None:
                print "biRRT found no solution -,,,-"

        RaveSetDebugLevel(DebugLevel.Debug) # set output level to debug
        openravepy.RaveLogInfo("Waiting for controller to finish")
        robot.WaitForController(0)
        robot.GetController().Reset()
        raw_input('Enter any key to quit. ')


