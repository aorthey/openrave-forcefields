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
        def __init__(self, env_in):
                self.env = env_in
                pass

        ### Cij \in { SE(3), True/False } True/False: is fixed or not
        def printContact( self, C ):
                if C is None:
                        print "free"
                else:
                        print "contact at",C[0:3,3]
                        print C

        def fromContactTransform( self, robot, Cll, Crl, Cla, Cra): 
                print "--------------------------------------------------"
                print "GIK from surface contacts"
                print "--------------------------------------------------"
                print "LEFT FOOT  : ", self.printContact(Cll)
                print "RIGHT FOOT : ", self.printContact(Crl)
                print "LEFT HAND  : ", self.printContact(Cla)
                print "RIGHT HAND : ", self.printContact(Cra)

                q_old = self.env.surrender_pos
                N = q_old.shape[0]

                with self.env.env:
                        cbirrt = CBiRRT(self.env.env, self.env.robot_name)
                        try:
                                [qlimL,qlimU]=robot.GetActiveDOFLimits()
                                self.env.EnforceLimits(q_old, qlimL, qlimU)
                                robot.SetActiveDOFValues(q_old)

                                cog = robot.GetCenterOfMass()
                                ### maniptm_list : all endeffectors which are fixed
                                ### (even if not on a surface)
                                ### support_list: all endeffectors in contact with a surface
                                support_list = []
                                maniptm_list = []

                                #maniptm_list = [('l_leg',left_leg_tf),('r_leg',right_leg_tf)]
                                #support_list.append(('l_leg',SURFACE_FRICTION))
                                #support_list.append(('r_leg',SURFACE_FRICTION))
                                #maniptm_list.append(('r_leg',right_leg_tf))

                                assert(Cll.shape == (4,4))
                                if Cll is not None:
                                        maniptm_list.append(('l_leg',Cll))
                                        support_list.append(('l_leg',SURFACE_FRICTION))
                                if Crl is not None:
                                        maniptm_list.append(('r_leg',Crl))
                                        support_list.append(('r_leg',SURFACE_FRICTION))
                                if Cla is not None:
                                        maniptm_list.append(('l_arm',Cla))
                                        support_list.append(('l_arm',SURFACE_FRICTION))
                                if Cra is not None:
                                        maniptm_list.append(('r_arm',Cra))
                                        support_list.append(('r_arm',SURFACE_FRICTION))

                                print "SUPPORT       :",support_list
                                #print "MANIP FIXED   :",maniptm_list
                                #obstacle_list = [('floor',(0,0,1))]
                                #obstacle_list = [('floor',(0,0,1))]
                                F = np.zeros((3))
                                F += np.array((0,0,-9.80)) ## boston-style gravity
                                #F += np.array((0,0.01,0))
                                q_gik = cbirrt.DoGeneralIK(
                                                #execute=False,
                                                #movecog=cog,
                                                gravity=F.tolist(),
                                                returnclosest=True,
                                                #checkcollisionlink=['l_foot','r_foot'],
                                                #obstacles=obstacle_list,
                                                maniptm=maniptm_list,
                                                support=support_list,
                                                printcommand=True)

                                if q_gik is None:
                                        print "------------------------------------------------------- "
                                        print "No solution found GIK"
                                        print "ZHEIGHT FEET (L,R):",Cll[2,3],Crl[2,3],support_list
                                        print "LEFT FOOT POS     :",Cll[0:3,3]
                                        print "RIGHT FOOT POS    :",Crl[0:3,3]
                                        print "------------------------------------------------------- "
                                        sys.exit(0)

                                robot.SetActiveDOFValues(q_gik)
                                #left_leg_tf_gik = robot.GetManipulator('l_leg').GetTransform()
                                #right_leg_tf_gik = robot.GetManipulator('r_leg').GetTransform()

                        except Exception as e:
                                print "Exception in GIK fromContactTransform"
                                print e
                                sys.exit(0)
                return q_gik
