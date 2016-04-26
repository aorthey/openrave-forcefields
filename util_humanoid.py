import numpy as np
import os.path
from util import Rz
import sys
import math
import time
from cbirrtpy import *

def waitrobot(robot):
    """busy wait for robot completion"""
    while not robot.GetController().IsDone():
        time.sleep(0.1)

def COM_interpolate(c1,c2,M):
        COM_linear = np.zeros((3,M))
        i=0
        while i<M:
                a = i/float(M)
                COM_linear[:,i] = c1*(1-a) + a*c2
                i=i+1
        return COM_linear

def GetStepLength(COM_project, startPos, MAX_FOOT_STEP_LENGTH):
        d=0.0
        M = COM_project.shape[1]
        i=startPos
        FOOT_STEP_LENGTH = MAX_FOOT_STEP_LENGTH
        while d<=MAX_FOOT_STEP_LENGTH:
                if i >= M:
                        FOOT_STEP_LENGTH = np.linalg.norm(COM_project[:,M-1] - COM_project[:,startPos])
                        break
                else:
                        d = np.linalg.norm(COM_project[:,i] - COM_project[:,startPos])
                        i=i+1
        return [FOOT_STEP_LENGTH,i-1]

def COM_compute_zig_zag_motion(COM_linear, env):
        ##assume that we start with both feet spaced apart at start pos

        M = COM_linear.shape[1]
        ## project COM onto surface
        COM_project = np.zeros((3,M))
        COM_project[0:2,:]=COM_linear[0:2,:]

        ## GET derivative along path
        der = np.zeros((3,M))
        nrml = np.zeros((3,M))
        Lf = np.zeros((3,2*M))
        Rf = np.zeros((3,2*M))


        ### in (cm)
        FOOT_SPACING = 0.25
        MAX_FOOT_STEP_LENGTH = 0.4
        FOOT_STEP_HEIGHT = 0.1

        #######################################################################
        ### compute der and normal to com path projhected onto floor
        #######################################################################
        for i in range(0,M-1):
                der[:,i] = COM_project[:,i+1] - COM_project[:,i]
                nrml[:,i] = np.dot(Rz(math.pi/2),der[:,i])
                nrml[:,i] /= np.linalg.norm(nrml[:,i])

        nrml[:,-1]=nrml[:,-2]
        der[:,-1]=der[:,-2]

        Lsupport = True
        Lposition = 0
        Rposition = 0
        i = 0

        leftFoot = []
        rightFoot = []

        #######################################################################
        ### compute starting foot steps
        #######################################################################
        Lf = COM_project[:,0] + 0.5*FOOT_SPACING*nrml[:,0]
        leftFoot.append(Lf)
        Rf = COM_project[:,0] - 0.5*FOOT_SPACING*nrml[:,0]
        rightFoot.append(Rf)

        #######################################################################
        ### make a half step with right foot
        #######################################################################
        d=0.0
        while d < MAX_FOOT_STEP_LENGTH/2:
                d = np.linalg.norm(COM_project[:,i] - COM_project[:,0])
                i=i+1
        Rf = COM_project[:,i-1] - 0.5*FOOT_SPACING*nrml[:,i-1]
        rightFoot.append(Rf)

        #######################################################################
        ### start alternating footsteps
        #######################################################################
        Rposition = i-1
        Lposition = 0

        Lsupport = False
        while (Rposition < M-1) and (Lposition < M-1):
                if Lsupport:
                        ## Lsupport, compute next right foot
                        [FOOT_STEP_LENGTH, Rposition] = GetStepLength(COM_project, Rposition, MAX_FOOT_STEP_LENGTH)
                        Rf = COM_project[:,Rposition] - 0.5*FOOT_SPACING*nrml[:,Rposition]
                        rightFoot.append(Rf)
                        print "RIGHT:",Rposition,Rf
                        Lsupport = False
                else:
                        ## Rsupport, compute next left foot
                        [FOOT_STEP_LENGTH, Lposition] = GetStepLength(COM_project, Lposition, MAX_FOOT_STEP_LENGTH)
                        Lf = COM_project[:,Lposition] + 0.5*FOOT_SPACING*nrml[:,Lposition]
                        print "LEFT:",Lposition,Lf
                        Lsupport = True
                        leftFoot.append(Lf)

        Lf = np.array(leftFoot)
        Rf = np.array(rightFoot)

        handles=[]
        handles.append(env.env.plot3(points=Lf,
                           pointsize=25,
                           colors=np.array(((1.0,0.0,0.0,0.8)))))
        handles.append(env.env.plot3(points=Rf,
                           pointsize=25,
                           colors=np.array(((0.0,1.0,0.0,0.8)))))

        FFn = Lf.shape[0] + Rf.shape[0]
        COM = np.zeros((FFn,3))
        ctr=0
        for i in range(0,min(Lf.shape[0],Rf.shape[0])):
                COM[ctr,:] = Rf[i,:]
                ctr+=1
                COM[ctr,:] = Lf[i,:]
                ctr+=1

        print Rf.shape[0],Lf.shape[0]
        if Rf.shape[0] > Lf.shape[0]:
                COM[ctr,:] = Rf[-1,:]
        elif Rf.shape[0]<Lf.shape[0]:
                COM[ctr,:] = Lf[-1,:]

        COM=COM.T

        #######################################################################
        ### interpolate COM
        #######################################################################

        from scipy.interpolate import interp1d,splev,splrep,splprep
        tvec = np.linspace(0,1,M)
        [trajectory,tmp] = splprep(COM,k=3,s=0.0)
        COM_traj = np.array([splev(t,trajectory) for t in tvec]).T

        print COM_traj.shape
        COM_traj[2,:] = COM_linear[2,:]
        handles.append(env.env.drawlinestrip(points=COM_linear.T,
                           linewidth=8,
                           colors=np.array((0.0,1.0,0.0,1.0))))
        handles.append(env.env.drawlinestrip(points=COM_traj.T,
                           linewidth=8,
                           colors=np.array((1.0,0.0,1.0,0.5))))

        raw_input('Press <ENTER> to continue.')


        #handles.append(env.env.drawlinestrip(points=Lf,
        #                   linewidth=6,
        #                   colors=np.array(((1.0,0.0,0.0,0.8)))))
        #handles.append(env.env.drawlinestrip(points=Rf,
        #                   linewidth=6,
        #                   colors=np.array(((0.0,1.0,0.0,0.8)))))

        #while i < M:
        #        if Lsupport:
        #                ### left foot is on support, only move right foot
        #                d= 0.0

        #                #######################################################
        #                ## CHECK THAT STEP LENGTH DOES NOT OVERSHOOT TARGET
        #                #######################################################
        #                FOOT_STEP_LENGTH = MAX_FOOT_STEP_LENGTH

        #                while d<=MAX_FOOT_STEP_LENGTH:
        #                        if i >= M:
        #                                FOOT_STEP_LENGTH = np.linalg.norm(COM_project[:,M-1] - COM_project[:,Rposition])
        #                                break
        #                        else:
        #                                d = np.linalg.norm(COM_project[:,i] - COM_project[:,Rposition])
        #                                i=i+1

        #                i = Rposition
        #                d = 0.0
        #                #######################################################
        #                ## PARABOLA ALONG STEP TO COMPUTE Z-COORDINATE
        #                #######################################################
        #                ## compute parabola
        #                ## -a*x^2 + b
        #                ## b = FOOT_STEP_HEIGHT
        #                ## a = 4*b/(FOOT_STEP_SIZE**2)

        #                while d < FOOT_STEP_LENGTH:
        #                        b = FOOT_STEP_HEIGHT
        #                        a = 4*FOOT_STEP_HEIGHT/(FOOT_STEP_LENGTH*FOOT_STEP_LENGTH)
        #                        z = max(-a*(d-FOOT_STEP_LENGTH/2)**2 + b,0.0)
        #                        Lf[:,i] = COM_project[:,Lposition] + 0.5*FOOT_SPACING*nrml[:,Lposition]
        #                        Rf[:,i] = COM_project[:,i] - 0.5*FOOT_SPACING*nrml[:,i]
        #                        Rf[2,i] = z + COM_project[2,i]

        #                        d = np.linalg.norm(COM_project[:,i] - COM_project[:,Rposition])
        #                        i=i+1
        #                ### one security step to make sure that we rest at z=0
        #                Lf[:,i-1] = COM_project[:,Lposition] + 0.5*FOOT_SPACING*nrml[:,Lposition]
        #                Rf[:,i-1] = COM_project[:,i-1] - 0.5*FOOT_SPACING*nrml[:,i-1]
        #                Lsupport = False
        #                print "RIGHT FOOT:",Rposition,"->",i
        #                Rposition = i
        #        else:
        #                ### right foot is on support, only move left foot
        #                d= 0.0
        #                i = Lposition

        #                #######################################################
        #                ## CHECK THAT STEP LENGTH DOES NOT OVERSHOOT TARGET
        #                #######################################################
        #                FOOT_STEP_LENGTH = MAX_FOOT_STEP_LENGTH

        #                while d<=MAX_FOOT_STEP_LENGTH:
        #                        if i >= M:
        #                                FOOT_STEP_LENGTH = np.linalg.norm(COM_project[:,M-1] - COM_project[:,Lposition])
        #                                break
        #                        else:
        #                                d = np.linalg.norm(COM_project[:,i] - COM_project[:,Lposition])
        #                                i=i+1

        #                i = Lposition
        #                d = 0.0
        #                #######################################################
        #                ## PARABOLA ALONG STEP TO COMPUTE Z-COORDINATE
        #                #######################################################
        #                #i = Lposition
        #                while d < FOOT_STEP_LENGTH:
        #                        b = FOOT_STEP_HEIGHT
        #                        a = 4*FOOT_STEP_HEIGHT/(FOOT_STEP_LENGTH*FOOT_STEP_LENGTH)
        #                        z = max(-a*(d-FOOT_STEP_LENGTH/2)**2 + b,0.0)
        #                        Rf[:,i] = COM_project[:,Rposition] - 0.5*FOOT_SPACING*nrml[:,Rposition]

        #                        Lf[:,i] = COM_project[:,i] + 0.5*FOOT_SPACING*nrml[:,i]
        #                        Lf[2,i] = z + COM_project[2,i]

        #                        d = np.linalg.norm(COM_project[:,i] - COM_project[:,Lposition])
        #                        i=i+1
        #                ### one security step to make sure that we rest at z=0
        #                Rf[:,i-1] = COM_project[:,Rposition] - 0.5*FOOT_SPACING*nrml[:,Rposition]
        #                Lf[:,i-1] = COM_project[:,i-1] + 0.5*FOOT_SPACING*nrml[:,i-1]
        #                Lsupport = True
        #                print "LEFT FOOT:",Lposition,"->",i
        #                Lposition = i




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
                        feet_pos = np.zeros((6,M))
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
                                        
                                        feet_pos[0:3,i]=left_leg_tf[0:3,3]
                                        feet_pos[3:6,i]=right_leg_tf[0:3,3]
                                        #zfoot[0,i]=left_leg_tf[2,3]
                                        #zfoot[1,i]=right_leg_tf[2,3]

                                        print "ZHEIGHT FEET:",left_leg_tf[2,3],right_leg_tf[2,3],support_list

                                        cog = COM_path[:,i]
                                        #obstacle_list = [('floor',(0,0,1))]
                                        obstacle_list = [('floor',(0,0,1))]
                                        F = np.zeros((3))
                                        F += np.array((0,0,-9.81))
                                        #F += np.array((0,0.01,0))
                                        q_res = cbirrt.DoGeneralIK(
                                                        movecog=cog,
                                                        gravity=F.tolist(),
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
                #plot(arange(0,M),feet_pos[0,:],'-r',linewidth=4)
                #plot(arange(0,M),feet_pos[1,:],'-g',linewidth=4)
                #plot(arange(0,M),feet_pos[2,:],'-b',linewidth=4)
                #plot(arange(0,M),feet_pos[3,:],'--r',linewidth=4)
                #plot(arange(0,M),feet_pos[4,:],'--g',linewidth=4)
                #plot(arange(0,M),feet_pos[5,:],'--b',linewidth=4)
                #plt.show()
                np.save(q_gik_fname,q_gik)
                np.save(COM_gik_fname,COM_gik)

        else:
                q_gik = np.load(q_gik_fname+'.npy')
                COM_gik = np.load(COM_gik_fname+'.npy')

        return [q_gik, COM_gik]
