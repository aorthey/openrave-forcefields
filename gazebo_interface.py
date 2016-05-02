import numpy as np
import math
import sys
import copy
import time
import pickle
from openravepy import *

import rospy
import rosbag
import actionlib
import tf

from control_msgs.msg import FollowJointTrajectoryAction, FollowJointTrajectoryActionGoal
from sensor_msgs.msg import JointState
from gazebo_msgs.msg import ModelStates
from trajectory_msgs.msg import JointTrajectory
from trajectory_msgs.msg import JointTrajectoryPoint
from vigir_humanoid_control_msgs.msg import ChangeControlModeAction, ChangeControlModeActionGoal
from vigir_footstep_planning_msgs.msg import StepPlan, Step, Feet, Foot

from tf.transformations import compose_matrix
from util_humanoid import GetSignedAngle

#mode = 'whole_body'
mode = 'footstep'

class Stance():
        left_leg = []
        right_leg = []
        def __init__(self,L,dL,R,dR):

                self.left_leg = np.zeros(6)
                self.left_leg[0:3] = L
                left_yaw = GetSignedAngle(dL)
                self.left_leg[5] = left_yaw

                self.right_leg = np.zeros(6)
                self.right_leg[0:3] = R
                right_yaw = GetSignedAngle(dR)
                self.right_leg[5] = right_yaw

class GazeboInterface():

        def __init__(self):
                pass

        def executeFootStepPath(self, Lf, dLf, Rf, dRf):

                rpos = 0
                lpos = 0

                Nr = len(Rf)
                Nl = len(Lf)
                path = []

                s = Stance(Lf[lpos,:],dLf[lpos,:],Rf[rpos,:],dRf[rpos,:])
                path.append(s)

                while rpos < Nr-1 and lpos < Nl-1:
                        rpos += 1
                        s = Stance(Lf[lpos,:],dLf[lpos,:],Rf[rpos,:],dRf[rpos,:])
                        path.append(s)
                        lpos += 1
                        s = Stance(Lf[lpos,:],dLf[lpos,:],Rf[rpos,:],dRf[rpos,:])
                        path.append(s)


                for i in range(0,len(path)):
                        s = path[i]
                        print s.left_leg[0:3],s.right_leg[0:3]

        	########################################
        	# path is a list of robot stance object, 
        	# and each stance object has two member: 
        	# left_leg and right_leg, each of them 
        	# is the x y z roll pitch yaw of the 
        	# foot pose. (x,y,z in meter, roll, 
        	# pitch, yaw in degrees)
        	########################################
                rospy.init_node('whole_body_traj_client')

                print "stand up"
                self.change_mode('stand')
                time.sleep(0.5)
                print "step"
                self.change_mode('step')
                time.sleep(0.5)

                initial_state = path[0]
                goal_state = path[len(path)-1]

                initial_left_leg = Foot()
                initial_left_leg.foot_index = initial_left_leg.LEFT

                initial_right_leg = Foot()
                initial_right_leg.foot_index = initial_right_leg.RIGHT

                goal_left_leg = Foot()
                goal_left_leg.foot_index = goal_left_leg.LEFT

                goal_right_leg = Foot()
                goal_right_leg.foot_index = goal_right_leg.RIGHT

                foot_message = [initial_left_leg,initial_right_leg,goal_left_leg,goal_right_leg]
                state_object = [initial_state.left_leg,initial_state.right_leg,goal_state.left_leg,goal_state.right_leg]

                for i in range(4):
                        message = foot_message[i]
                        state = state_object[i]

                        message.pose.position.x = state[0] * 1.062 # magic number
                        message.pose.position.y = state[1] * 1.062
                        message.pose.position.z = state[2] * 1.062
                        quat = quatFromAxisAngle([state[3]*math.pi/180.0, state[4]*math.pi/180.0, state[5]*math.pi/180.0])
                        message.pose.orientation.x = quat[1]
                        message.pose.orientation.y = quat[2]
                        message.pose.orientation.z = quat[3]
                        message.pose.orientation.w = quat[0]

                initial_feet_poses = Feet()
                initial_feet_poses.left = initial_left_leg
                initial_feet_poses.right = initial_right_leg

                goal_feet_poses = Feet()
                goal_feet_poses.left = goal_left_leg
                goal_feet_poses.right = goal_right_leg

                single_step_plan = StepPlan()

                single_step_plan.start = initial_feet_poses
                single_step_plan.goal = goal_feet_poses

                the_step = Step()
                the_step.step_index = 0
                the_step.cost = 0
                the_step.risk = 0
                the_step.valid = True
                the_step.colliding = False
                the_step.locked = False
                the_step.modified = False
                the_step.sway_duration = 1.5
                the_step.step_duration = 0.5
                the_step.swing_height = 0.05

                print "create footstep path"
                for i in range(len(path)):
                        print i
                        if(i>0):
                                the_step.step_index = i-1
                                if(abs(path[i].left_leg[0] - path[i-1].left_leg[0]) > abs(path[i].right_leg[0] - path[i-1].right_leg[0]) or 
                                   abs(path[i].left_leg[1] - path[i-1].left_leg[1]) > abs(path[i].right_leg[1] - path[i-1].right_leg[1])):
                                        state = path[i].right_leg
                                        the_step.foot.foot_index = the_step.foot.RIGHT
                                        final_step = goal_left_leg
                                else:
                                        state = path[i].left_leg
                                        the_step.foot.foot_index = the_step.foot.LEFT
                                        final_step = goal_right_leg
                                    

                                the_step.foot.pose.position.x = state[0] * 1.062 # magic number
                                the_step.foot.pose.position.y = state[1] * 1.062
                                the_step.foot.pose.position.z = state[2] * 1.062
                                quat = quatFromAxisAngle([state[3]*math.pi/180.0, state[4]*math.pi/180.0, state[5]*math.pi/180.0])
                                the_step.foot.pose.orientation.x = quat[1]
                                the_step.foot.pose.orientation.y = quat[2]
                                the_step.foot.pose.orientation.z = quat[3]
                                the_step.foot.pose.orientation.w = quat[0]

                                single_step_plan.steps.append(copy.deepcopy(the_step))

                the_step.foot = final_step
                single_step_plan.steps.append(copy.deepcopy(the_step))

                single_step_plan.mode = 5

                print "trying to publish footsteps"
                pub = rospy.Publisher('/vigir/footstep_planning/add_footstep_topic', StepPlan, queue_size=10)
                time.sleep(0.5)
                pub.publish(single_step_plan)
                print "published"

        def change_mode(self,new_mode, wait_duration=3.0):
                client = actionlib.SimpleActionClient('/mode_controllers/control_mode_controller/change_control_mode',ChangeControlModeAction)

                client.wait_for_server()

                mode_msg = ChangeControlModeActionGoal()
                mode_msg.goal.mode_request = new_mode

                client.send_goal(mode_msg.goal)
                return client.wait_for_result(rospy.Duration.from_sec(wait_duration))
