import numpy as np
import sys
import copy
from numpy import sqrt,sin,cos,pi
from qcqp import *
import scipy
import cvxopt
import matplotlib.pyplot as plt
from cvxpy import *
from cvxopt import matrix, solvers
from util import PrintNumpy,inf

FILENAME = 'ddrive4'
DYNAMICAL_SYTEM_NAME = 'non-holonomic differential drive (deform)'
system_name = DYNAMICAL_SYTEM_NAME

class DynamicalSystem():
        env_ptr = None

        def __init__(self, env):
                self.env_ptr = env
                self.robot = env.GetRobot()

                links = self.robot.GetLinks()
                Nl = len(links)

                F = np.array((0,-1.5,0))

                wrenchmap = {x: np.zeros(6) for x in np.arange(Nl)}
                for link in links:
                        idx = link.GetIndex()
                        print wrenchmap[idx]

                #forcetorquemap = {x: np.zeros(6) for x in np.arange(Nl)}
                #forcetorquemap[iTorso]=Ftorso
                #for joint in self.robot.GetJoints():
                        #print joint

                dofaccel = np.zeros((self.robot.GetDOF()))
                torqueconfiguration, torquecoriolis, torquegravity = self.robot.ComputeInverseDynamics(dofaccel,None,returncomponents=True)
                sys.exit(0)

        def GetControlMatrix(p):
                Ndim = p.shape[0]
                Kdim = 3
                R = np.zeros((Ndim,Kdim))

                t = p[3]
                R[0,:] = np.array((cos(t),-sin(t),0.0))
                R[1,:] = np.array((sin(t),cos(t),0.0))
                R[2,:] = np.array((0.0,0.0,0.0))
                R[3,:] = np.array((0.0,0.0,1.0))

                return R

        #def GetMassMatrix(p):
                #body.SetDOFValues(dofvalues,range(body.GetDOF()),checklimits=True)
                #M = zeros((body.GetDOF(),body.GetDOF()))
                #for index in range(body.GetDOF()):
                #        testaccel = zeros(body.GetDOF())
                #        testaccel[index] = 1.0
                #        M[:,index] = body.ComputeInverseDynamics(testaccel)

       # def GetCoriolisMatrix(p):
       #     body.SetDOFValues(dofvalues,range(body.GetDOF()),checklimits=True)
       #     M = zeros((body.GetDOF(),body.GetDOF()))
       #     for index in range(body.GetDOF()):
       #         testaccel = zeros(body.GetDOF())
       #         testaccel[index] = 1.0
       #         M[:,index] = body.ComputeInverseDynamics(testaccel)
       #         return M

        #Td = traj.GetDuration()
        #forcetorquemap = {x: np.zeros(6) for x in np.arange(Nl)}
        #### create forcetorquemap dictionary
        #link = robot.GetLink('torso')
        #iTorso = link.GetIndex()
        #Ftorso = np.array((0,1.0,0,0,0,0))
        #Ffoot = np.array((0,0.0,1.0,0,0,0))
        #forcetorquemap[iTorso]=Ftorso

        #dt = 0.1
        #t = 0.0
        #qold = traj.Sample(0,active_dofs)
        #qext = np.zeros((robot.GetDOF()))
        #zlvec = []
        #zrvec = []
        #ictr=0

        #while t < Td:
        #        qt = traj.Sample(t,active_dofs)

        #        with env.env:
        #                ##### add contact forces
        #                rlink = robot.GetManipulator('r_leg')
        #                llink = robot.GetManipulator('l_leg')
        #                zr = rlink.GetTransform()[2,3]
        #                zl = llink.GetTransform()[2,3]
        #                zlvec.append(zl)
        #                zrvec.append(zr)

        #                lsole = robot.GetLink('l_sole')
        #                zsl = lsole.GetTransform()[2,3]
        #                rsole = robot.GetLink('r_sole')
        #                zsr = rsole.GetTransform()[2,3]

        #                iLlink = lsole.GetIndex()
        #                iRlink = rsole.GetIndex()

        #                llink = robot.GetLink('l_foot')
        #                rlink = robot.GetLink('l_foot')
        #                env.DrawArrow(robot.GetLink('torso').GetGlobalCOM(), Ftorso[0:3], deleteOld=True)

        #                if zl < 0.01:
        #                        forcetorquemap[iLlink]=Ffoot
        #                        env.DrawArrow(llink.GetGlobalCOM(), -Ffoot[0:3])
        #                else:
        #                        forcetorquemap[iLlink]=np.zeros(6)
        #                        env.DrawArrow(llink.GetGlobalCOM(), np.zeros(3))
        #                if zr < 0.01:
        #                        forcetorquemap[iRlink]=Ffoot
        #                else:
        #                        forcetorquemap[iRlink]=np.zeros(6)
