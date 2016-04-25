import numpy as np
import os.path
import sys
import time
from cbirrtpy import *

def waitrobot(robot):
    """busy wait for robot completion"""
    while not robot.GetController().IsDone():
        time.sleep(0.1)

def COM_from_path(rave_path, robot, env):
        active_dofs = robot.GetActiveConfigurationSpecification()
        N = len(robot.GetActiveDOFValues())
        M = rave_path.GetNumWaypoints()
        COM_original = np.zeros((3,M))
        COM_gik = np.zeros((3,M))
        q_original = np.zeros((N, M))
        q_gik = np.zeros((N, M))

        i = 0
        [qlimL,qlimU]=robot.GetActiveDOFLimits()

        print "Waypoints:",M," - Dimension:",N
        with env.env:
                while i < M:
                        q_original[:,i] = rave_path.GetWaypoint(i,active_dofs)
                        #q_original[:,i] = env.EnforceLimits(q_original[:,i],qlimL,qlimU)
                        robot.SetActiveDOFValues(q_original[:,i])
                        COM_original[:,i] = robot.GetCenterOfMass()
                        i = i+1

        return [q_original, COM_original]

def visualize_configurations(q_original, robot, env):
        N = q_original.shape[0]
        M = q_original.shape[1]

        i = 0
        while i < M:
                with env.env:
                        print "------------------------------------------------------------------"
                        print "------------ WAYPOINT",i,"/",M,": recompute GIK ------------ "
                        print "------------------------------------------------------------------"
                        robot.SetActiveDOFValues(q_original[:,i])
                        left_leg_tf = robot.GetManipulator('l_leg').GetTransform()
                        right_leg_tf = robot.GetManipulator('r_leg').GetTransform()
                        waitrobot(robot)
                time.sleep(0.1)
                i=i+1

def GIK_from_COM(COM_path, q_original, robot, env, recompute=False, DEBUG=False):
        q_gik_fname = 'tmp/q_gik.numpy'
        COM_gik_fname = 'tmp/COM_gik.numpy'

        N = q_original.shape[0]
        M = q_original.shape[1]

        q_gik = np.zeros((N,M))
        Z_FOOT_CONTACT = 0.002
        COM_gik = np.zeros((3,M))
        if not os.path.isfile(q_gik_fname+'.npy') or recompute:
                i = 0
                with env.env:
                        cbirrt = CBiRRT(env.env, env.robot_name)
                        zfoot = np.zeros((2,M))
                        while i < M:
                                if DEBUG:
                                        print "------------------------------------------------------------------"
                                        print "------------ WAYPOINT",i,"/",M,": recompute GIK ------------ "
                                        print "------------------------------------------------------------------"
                                try:
                                        robot.SetActiveDOFValues(q_original[:,i])

                                        left_leg_tf = robot.GetManipulator('l_leg').GetTransform()
                                        right_leg_tf = robot.GetManipulator('r_leg').GetTransform()
                                        left_arm_tf = robot.GetManipulator('l_arm').GetTransform()
                                        right_arm_tf = robot.GetManipulator('r_arm').GetTransform()

                                        #maniptm_list = [('l_leg',left_leg_tf),('r_leg',right_leg_tf),('l_arm',left_arm_tf),('r_arm',right_arm_tf)]
                                        #maniptm_list = [('l_leg',left_leg_tf),('r_leg',right_leg_tf)]
                                        friction_coefficient = 0.8

                                        #support_list = [('l_leg',friction_coefficient),
                                                        #('r_leg',friction_coefficient),
                                                        #('l_arm',friction_coefficient),
                                                        #('r_arm',friction_coefficient)]

                                        support_list = []
                                        maniptm_list = []
                                        if left_leg_tf[2,3] < Z_FOOT_CONTACT:
                                                support_list.append(('l_leg',friction_coefficient))
                                                maniptm_list.append(('l_leg',left_leg_tf))
                                        if right_leg_tf[2,3] < Z_FOOT_CONTACT:
                                                support_list.append(('r_leg',friction_coefficient))
                                                maniptm_list.append(('r_leg',right_leg_tf))
                                        
                                        zfoot[0,i]=left_leg_tf[2,3]
                                        zfoot[1,i]=right_leg_tf[2,3]

                                        print "ZHEIGHT FEET:",left_leg_tf[2,3],right_leg_tf[2,3],support_list
                                        cog = COM_path[:,i]
                                        #obstacle_list = [('floor',(0,0,1))]
                                        obstacle_list = [('floor',(0,0,1))]
                                        q_res = cbirrt.DoGeneralIK(
                                                        movecog=cog,
                                                        #checkcollisionlink=['l_foot','r_foot'],
                                                        obstacles=obstacle_list,
                                                        maniptm=maniptm_list,
                                                        support=support_list,
                                                        printcommand=False)

                                        if q_res is None:
                                                print "------------ WAYPOINT",i,"/",M,": recompute GIK ------------ "
                                                print "No solution found GIK"
                                                sys.exit(0)
                                        else:
                                                q_gik[:,i] = q_res
                                        #print "DIST q,q_gik:",np.linalg.norm(q_original[:,i]-q_gik[:,i])
                                        robot.SetActiveDOFValues(q_gik[:,i])
                                        COM_gik[:,i] = robot.GetCenterOfMass()
                                        left_leg_tf_gik = robot.GetManipulator('l_leg').GetTransform()
                                        right_leg_tf_gik = robot.GetManipulator('r_leg').GetTransform()

                                        print "CHANGE IN FOOT POS (L):"
                                        print left_leg_tf[0:3,3]
                                        print left_leg_tf_gik[0:3,3]
                                        print "CHANGE IN FOOT POS (R):"
                                        print right_leg_tf[0:3,3]
                                        print right_leg_tf_gik[0:3,3]

                                except Exception as e:
                                        print "Exception in GIK, waypoint",i,"/",M
                                        print e
                                        sys.exit(0)

                                if DEBUG:
                                        dcc = np.linalg.norm(cog[:,i]-COM_gik[:,i])
                                        print "ERROR GIK      :",dcc
                                        print "INPUT GIK  COM :",cog
                                        print "OUTPUT GIK COM :",COM_gik[:,i]

                                i = i+1

                #from pylab import *
                #plot(arange(0,M),zfoot[0,:],'-g')
                #plot(arange(0,M),zfoot[1,:],'-r')
                #plt.show()
                np.save(q_gik_fname,q_gik)
                np.save(COM_gik_fname,COM_gik)

        else:
                q_gik = np.load(q_gik_fname+'.npy')
                COM_gik = np.load(COM_gik_fname+'.npy')

        return [q_gik, COM_gik]
