import numpy as np
import os.path
from util import Rz, Rax
import sys
import math
from math import pi
import time
from cbirrtpy import *
from openravepy import *

SURFACE_FRICTION = 0.8

class GIKInterface():

        env = []
        handles = []
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

        def printFromToContact( self, Told, T ):
                if Told is None:
                        print "none"
                elif T is None:
                        print "none"
                else:
                        print "from: ",Told[0:3,3],
                        print "to: ",T[0:3,3]
        def checkValidTransform( self, T ):
                if T is None:
                        return
                else:
                        assert(T.shape == (4,4))
                        d = np.linalg.det(T[0:3,0:3])
                        if abs(d-1)>1e-5:
                                print "Not valid rotation matrix in transform"
                                print "rotation:"
                                print T[0:3,0:3]
                                print "det:",d
                                sys.exit(0)
                        return

        def DrawFootCone(self, p, T):
                e0 = np.array((0,0,0,1))
                ez = np.array((0,0,1,1))
                ct = np.dot(T,e0)[0:3]
                nt = np.dot(T,ez)[0:3] - ct
                A = self.env.env.drawarrow(p1=p,p2=p + 0.2*nt,linewidth=0.01,color=np.array((1,1,1)))
                self.handles.append(A)
                A = self.env.env.drawcone(pos=p,
                                direction=nt,height=0.4,aperture=0.1,color=np.array((0.6,0.3,0,0.5)))
                self.handles.append(A)

        def DrawHandCone(self, p, T, thetaNormal):
                Thand = np.eye(4)
                Thand[0:3,3] = T[0:3,3]

                e0 = np.array((0,0,0,1))
                ex = np.array((1,0,0,1))
                ey = np.array((0,1,0,1))
                ez = np.array((0,0,1,1))

                fpos = np.dot(T,e0)[0:3]
                fdir1 = np.dot(T,ex)[0:3] - fpos
                fdir2 = np.dot(T,ey)[0:3] - fpos
                fdir3 = np.dot(T,ez)[0:3] - fpos

                R = T[0:3,0:3]
                R = np.dot(Rax(pi/2, fdir3),R)
                R = np.dot(Rax(thetaNormal, fdir2),R)
                Thand[0:3,0:3]=R
                self.DrawFootCone(p, Thand)

        def DrawLeftHandCone(self, p, T):
                self.DrawHandCone(p, T, -pi/2)

        def DrawRightHandCone(self, p, T):
                self.DrawHandCone(p, T, -pi/2)

        def DrawCOM(self, com):
                ll = 0.5
                lw = 6
                xl = np.array((ll,0,0))
                yl = np.array((0,ll,0))
                zl = np.array((0,0,ll))
                ccom = np.array((1,0,0))

                P = np.array(((com+xl),(com-xl)))
                A=self.env.env.drawlinestrip(points=P,linewidth=lw,colors=ccom)
                self.handles.append(A)
                P = np.array(((com+yl),(com-yl)))
                A=self.env.env.drawlinestrip(points=P,linewidth=lw,colors=ccom)
                self.handles.append(A)
                P = np.array(((com+zl),(com-zl)))
                A=self.env.env.drawlinestrip(points=P,linewidth=lw,colors=ccom)
                self.handles.append(A)

                A = self.env.env.plot3(points=com,
                                pointsize=0.01,
                                colors=np.array((1,0,0)),
                                drawstyle=1)
                self.handles.append(A)

        def VisualizeConesCOM(self, com, Cll, Crl, Cla, Cra):
                self.DrawCOM(com)
                self.DrawFootCone(com, Cll)
                self.DrawFootCone(com, Crl)
                self.DrawLeftHandCone(com, Cla)
                self.DrawRightHandCone(com, Cra)

        def VisualizeCones(self, Cll, Crl, Cla, Cra):
                self.DrawFootCone(Cll[0:3,3], Cll)
                self.DrawFootCone(Crl[0:3,3], Crl)
                self.DrawLeftHandCone(Cla[0:3,3], Cla)
                self.DrawRightHandCone(Cra[0:3,3], Cra)


        def fromContactTransform( self, robot, Cll, Crl, Cla, Cra): 
                print "--------------------------------------------------"
                print "GIK from surface contacts"
                print "--------------------------------------------------"
                print "LEFT FOOT  : ", self.printContact(Cll)
                print "RIGHT FOOT : ", self.printContact(Crl)
                print "LEFT HAND  : ", self.printContact(Cla)
                print "RIGHT HAND : ", self.printContact(Cra)

                q_old = self.env.surrender_pos
                q_gik = None
                N = q_old.shape[0]

                self.checkValidTransform(Cll)
                self.checkValidTransform(Crl)
                self.checkValidTransform(Cla)
                self.checkValidTransform(Cra)

                with self.env.env:
                        cbirrt = CBiRRT(self.env.env, self.env.robot_name)
                        try:
                                [qlimL,qlimU]=robot.GetActiveDOFLimits()
                                self.env.EnforceLimits(q_old, qlimL, qlimU)
                                robot.SetActiveDOFValues(q_old)
                                Cll_old = robot.GetManipulator('l_leg').GetTransform()
                                Crl_old = robot.GetManipulator('r_leg').GetTransform()
                                Cla_old = robot.GetManipulator('l_arm').GetTransform()
                                Cra_old = robot.GetManipulator('r_arm').GetTransform()

                                cog = robot.GetCenterOfMass()

                                ### maniptm_list : all endeffectors which are fixed
                                ### (even if not on a surface)
                                ### support_list: all endeffectors in contact with a surface
                                support_list = []
                                maniptm_list = []

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
                                F = np.zeros((3))
                                F += np.array((0,0,-9.80)) ## boston-style gravity
                                q_gik = cbirrt.DoGeneralIK(
                                                movecog=cog,
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
                                self.env.env.GetCollisionChecker().SetCollisionOptions(CollisionOptions.Contacts)
                                report = CollisionReport()
                                ret = self.env.env.CheckCollision(robot, report=report)

                                if ret:
                                        print "------------------------------------------------------- "
                                        print "GIK found solution, but solution is in collision"
                                        print "------------------------------------------------------- "
                                        print '%d contacts'%len(report.contacts)
                                        positions = [c.pos for c in report.contacts]
                                        for c in report.contacts:
                                                p = c.pos
                                                self.handles.append(self.env.env.plot3(points=np.array(((p[0],p[1],p[2]))),
                                                                   pointsize=0.03,
                                                                   colors=np.array(((1.0,0.0,0.0,0.8))),
                                                                   drawstyle=1))
                                        print positions
                                        return None
                                        #sys.exit(0)

                                print "------------------------------------------------------- "
                                print "GIK found solution"
                                print "------------------------------------------------------- "
                                Cll = robot.GetManipulator('l_leg').GetTransform()
                                Crl = robot.GetManipulator('r_leg').GetTransform()
                                Cla = robot.GetManipulator('l_arm').GetTransform()
                                Cra = robot.GetManipulator('r_arm').GetTransform()

                                self.VisualizeCones(Cll, Crl, Cla, Cra)
                                self.VisualizeConesCOM(robot.GetCenterOfMass(), Cll, Crl, Cla, Cra)
                                #left_leg_tf_gik = robot.GetManipulator('l_leg').GetTransform()
                                #right_leg_tf_gik = robot.GetManipulator('r_leg').GetTransform()
                                print "LEFT FOOT  : ", 
                                self.printFromToContact(Cll_old,Cll)
                                print "RIGHT FOOT : ", 
                                self.printFromToContact(Crl_old, Crl)
                                print "LEFT HAND  : ",
                                self.printFromToContact(Cla_old,Cla)
                                print "RIGHT HAND : ", 
                                self.printFromToContact(Cra_old,Cra)

                        except Exception as e:
                                print "Exception in GIK fromContactTransform"
                                print e
                                sys.exit(0)
                print "return"
                return q_gik
