import numpy as np
import os.path
from util import Rz
import sys
import math
import time
from cbirrtpy import *

SURFACE_FRICTION = 0.8

class GIKInterface():

        env = []
        def __init__(env_in):
                self.env = env_in
                pass

        ### Cij \in { SE(3), True/False } True/False: is fixed or not
        def printContact( C ):
                if C is None:
                        print "free"
                else:
                        if C[1]:
                                S = C[0]
                                print "contact at",S[0:3,2]

        def GIK( Cll, Crl, Cla, Cra): 
                print "--------------------------------------------------"
                print "GIK from surface contacts"
                print "--------------------------------------------------"
                print "LEFT FOOT  : ", printContact(Cll)
                print "RIGHT FOOT : ", printContact(Crl)
                print "LEFT HAND  : ", printContact(Cla)
                print "RIGHT HAND : ", printContact(Cra)

                #N = q_original.shape[0]
                #M = q_original.shape[1]
                q_old = env.surrender_pos
                N = q_old.shape[0]

                with env.env:
                        cbirrt = CBiRRT(env.env, env.robot_name)
                        try:
                                [qlimL,qlimU]=robot.GetActiveDOFLimits()
                                env.EnforceLimits(q_old, qlimL, qlimU)
                                robot.SetActiveDOFValues(q_old)

                                if Cll is None:
                                        left_leg_tf = robot.GetManipulator('l_leg').GetTransform()
                                if Crl is None:
                                        right_leg_tf = robot.GetManipulator('r_leg').GetTransform()
                                if Cla is None:
                                        left_arm_tf = robot.GetManipulator('l_arm').GetTransform()
                                if Cra is None:
                                        right_arm_tf = robot.GetManipulator('r_arm').GetTransform()

                                support_list = []
                                maniptm_list = []

                                #[left_leg_tf, right_leg_tf] = createTransformFromPosDer( footpos[i,:], dfootpos[i,:] )

                                ### maniptm_list : all endeffectors which are fixed
                                ### (even if not on a surface)
                                ### support_list: all endeffectors in contact with a surface

                                maniptm_list = [('l_leg',left_leg_tf),('r_leg',right_leg_tf)]
                                support_list.append(('l_leg',SURFACE_FRICTION))
                                support_list.append(('r_leg',SURFACE_FRICTION))
                                #maniptm_list.append(('r_leg',right_leg_tf))
                                
                                cog = COM_path[:,i]

                                if len(support_list)<2:
                                        ## only one support foot -- make
                                        ## sure that COM is valid
                                        print "------------------------- COM adjustment made"
                                        if lfoot_contact:
                                                cog[0:2]= footpos[i,0:2]
                                        elif rfoot_contact:
                                                cog[0:2] = footpos[i,3:5]
                                        else:
                                                print "No foot contacts, not yet support"
                                                sys.exit(1)

                                print "SUPPORT       :",support_list
                                #obstacle_list = [('floor',(0,0,1))]
                                #obstacle_list = [('floor',(0,0,1))]
                                F = np.zeros((3))
                                F += np.array((0,0,-9.81))
                                #F += np.array((0,0.01,0))
                                q_res = cbirrt.DoGeneralIK(
                                                movecog=cog,
                                                gravity=F.tolist(),
                                                #returnclosest=True,
                                                #checkcollisionlink=['l_foot','r_foot'],
                                                #obstacles=obstacle_list,
                                                maniptm=maniptm_list,
                                                support=support_list,
                                                printcommand=False)

                                if q_res is None:
                                        print "------------ WAYPOINT",i,"/",M,": recompute GIK ------------ "
                                        print "No solution found GIK"
                                        print "ZHEIGHT FEET (L,R):",left_leg_tf[2,3],right_leg_tf[2,3],support_list
                                        print "COM PATH          :",COM_path[:,i]
                                        print "LEFT FOOT POS     :",left_leg_tf[0:3,3]
                                        print "RIGHT FOOT POS    :",right_leg_tf[0:3,3]
                                        sys.exit(0)
                                else:
                                        q_gik[:,i] = q_res
                                #print "DIST q,q_gik:",np.linalg.norm(q_original[:,i]-q_gik[:,i])
                                robot.SetActiveDOFValues(q_gik[:,i])
                                COM_gik[:,i] = robot.GetCenterOfMass()

                                print "CENTER OF MASS:",cog,"->",COM_gik[:,i]
                                left_leg_tf_gik = robot.GetManipulator('l_leg').GetTransform()
                                right_leg_tf_gik = robot.GetManipulator('r_leg').GetTransform()

                                print "CHANGE IN FOOT POS (L):"
                                print np.around(footpos[i,0:3],decimals=2)
                                print np.around(left_leg_tf_gik[0:3,3],decimals=2)

                                print "CHANGE IN FOOT POS (R):"
                                print np.around(footpos[i,3:6],decimals=2)
                                print np.around(right_leg_tf_gik[0:3,3],decimals=2)

                                q_old = q_gik[:,i]

                        except Exception as e:
                                print "Exception in GIK, waypoint",i,"/",M
                                print e
                                sys.exit(0)
