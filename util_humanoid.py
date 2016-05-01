import numpy as np
import os.path
from util import Rz
import sys
import math
import time
from cbirrtpy import *

COLOR_LEFT_FOOT = np.array((1.0,0.0,0.0,0.9))
COLOR_RIGHT_FOOT = np.array((0.0,1.0,0.0,0.9))
FOOT_WIDTH = 0.07
FOOT_LENGTH = 0.12
FOOT_SPACING = 0.2
MAX_FOOT_STEP_LENGTH = 0.3
FOOT_STEP_HEIGHT = 0.05
#Z_FOOT_CONTACT = 0.002
Z_FOOT_CONTACT = 0.0001
SURFACE_FRICTION = 0.9
np.set_printoptions(precision=2)

def waitrobot(robot):
    """busy wait for robot completion"""
    while not robot.GetController().IsDone():
        time.sleep(0.1)

def findNextCOMsegment(istart, COM, Fpos):
        ## find nearest point on COM to left foot
        ictr = istart

        d = 100.0
        dold = 1000.0
        while d < dold or d > 0.05:
                dold = d
                if ictr >= COM.shape[1]:
                        return [istart, COM.shape[1]-1]
                d = np.linalg.norm(COM[:,ictr]-Fpos)
                #print d,ictr
                ictr +=1

        ictr -= 2
        d = np.linalg.norm(COM[:,ictr]-Fpos)

        return [istart,ictr]


def COM_interpolate(c1,c2,M):
        COM_linear = np.zeros((3,M))
        i=0
        while i<M:
                a = i/float(M)
                COM_linear[:,i] = c1*(1-a) + a*c2
                i=i+1
        return COM_linear

def GetStepLengthFoot(sign, COM_project, nrml, startPos):
        d=0.0
        M = COM_project.shape[1]
        i=startPos
        FOOT_STEP_LENGTH = MAX_FOOT_STEP_LENGTH
        while d<=MAX_FOOT_STEP_LENGTH:
                if i >= M:
                        pstart = COM_project[:,M-1] + sign*0.5*FOOT_SPACING*nrml[:,M-1]
                        pd = COM_project[:,startPos] + sign*0.5*FOOT_SPACING*nrml[:,startPos]
                        d = np.linalg.norm(pd - pstart)
                        return [d, M-1]
                else:
                        pstart = COM_project[:,i] + sign*0.5*FOOT_SPACING*nrml[:,i]
                        pd = COM_project[:,startPos] + sign*0.5*FOOT_SPACING*nrml[:,startPos]
                        d = np.linalg.norm(pd - pstart)
                        i=i+1

        pstart = COM_project[:,i-2] + sign*0.5*FOOT_SPACING*nrml[:,i-2]
        pd = COM_project[:,startPos] + sign*0.5*FOOT_SPACING*nrml[:,startPos]
        d = np.linalg.norm(pd - pstart)
        return [d,i-2]

def GetStepLengthRightFoot(COM_project, nrml, startPos):
        return GetStepLengthFoot(-1, COM_project, nrml, startPos)

def GetStepLengthLeftFoot(COM_project, nrml, startPos):
        return GetStepLengthFoot(+1, COM_project, nrml, startPos)


def visualizeSingleFoot(env, F, dF, colorF):

        ndF=dF/np.linalg.norm(dF)
        bdF = np.dot(Rz(math.pi/2),ndF)

        X=np.zeros((3,7))
        X[:,0] = F + FOOT_LENGTH*ndF - FOOT_WIDTH*bdF
        X[:,1] = F + FOOT_LENGTH*ndF + FOOT_WIDTH*bdF
        X[:,2] = F - FOOT_LENGTH*ndF + FOOT_WIDTH*bdF
        X[:,3] = F - FOOT_LENGTH*ndF - FOOT_WIDTH*bdF
        X[:,4] = F + FOOT_LENGTH*ndF - FOOT_WIDTH*bdF
        X[:,5] = F + FOOT_LENGTH*ndF
        X[:,6] = F

        h = env.env.drawlinestrip(points=X.T,
                        linewidth = 6,
                        colors=colorF)
        return h

def visualizeFoot(env, F, dF, colorF):
        handles = []
        for i in range(0,F.shape[0]):
                h = visualizeSingleFoot(env, F[i,:],dF[i,:],colorF)
                handles.append(h)
        return handles

def interpolateFoot(N, f1, df1, f2, df2):
        ### create N points interpolation between f1,f2
        footstep_length = np.linalg.norm(f1-f2)
        if footstep_length>MAX_FOOT_STEP_LENGTH:
                print "cannot interpolate footstep"
                print f1,f2
                print footstep_length,">",MAX_FOOT_STEP_LENGTH
                sys.exit(1)
        #######################################################
        ## PARABOLA ALONG STEP TO COMPUTE Z-COORDINATE
        #######################################################
        d=0.0
        dstep = footstep_length/(N-1)
        ftvec = []
        dftvec = []
        ictr=0
        while d <= footstep_length+1e-5:
                t = d/footstep_length
                ft = (1-t)*f1 + t*f2
                dft = (1-t)*df1 + t*df2

                b = FOOT_STEP_HEIGHT
                a = 4*FOOT_STEP_HEIGHT/(footstep_length*footstep_length)
                z = max(-a*(d-footstep_length/2)**2 + b,0.0)

                #print ictr,"/",N,t,z
                ft[2] = z
                ftvec.append(ft)
                dftvec.append(dft)
                d = d+dstep
                ictr+=1
        return [np.array(ftvec),np.array(dftvec)]

handles = []

def GetFootPositionFromProjectedCOM( COM_project ):
        M = COM_project.shape[1]
        ## GET derivative along path
        der = np.zeros((3,M))
        nrml = np.zeros((3,M))
        Lf = np.zeros((3,2*M))
        Rf = np.zeros((3,2*M))

        #######################################################################
        ### compute der and normal to com path projected onto floor
        #######################################################################
        for i in range(0,M-1):
                der[:,i] = COM_project[:,i+1] - COM_project[:,i]
                nrml[:,i] = np.dot(Rz(math.pi/2),der[:,i])
                nrml[:,i] /= np.linalg.norm(nrml[:,i])

        nrml[:,-1]=nrml[:,-2]
        der[:,-1]=der[:,-2]

        leftFoot = []
        rightFoot = []
        leftFootDer = []
        rightFootDer = []

        #######################################################################
        ### compute starting foot steps
        #######################################################################
        Lf = COM_project[:,0] + 0.5*FOOT_SPACING*nrml[:,0]
        leftFoot.append(Lf)
        leftFootDer.append(der[:,0])
        Rf = COM_project[:,0] - 0.5*FOOT_SPACING*nrml[:,0]
        rightFoot.append(Rf)
        rightFootDer.append(der[:,0])

        #######################################################################
        ### make a half step with right foot
        #######################################################################
        d=0.0
        i = 0
        while d < MAX_FOOT_STEP_LENGTH/2:
                d = np.linalg.norm(COM_project[:,i] - COM_project[:,0])
                i=i+1
        Rf = COM_project[:,i-1] - 0.5*FOOT_SPACING*nrml[:,i-1]
        rightFoot.append(Rf)
        rightFootDer.append(der[:,i-1])

        #######################################################################
        ### start alternating footsteps
        #######################################################################
        Rposition = i-1
        Lposition = 0

        Lsupport = False
        while (Rposition < M-1) and (Lposition < M-1):
                if Lsupport:
                        ## Lsupport, compute next right foot
                        [FOOT_STEP_LENGTH, Rposition] = GetStepLengthRightFoot(COM_project, nrml, Rposition)
                        Rf = COM_project[:,Rposition] - 0.5*FOOT_SPACING*nrml[:,Rposition]
                        dRf = der[:,Rposition]
                        Lsupport = False
                        rightFoot.append(Rf)
                        rightFootDer.append(dRf)
                        print "R:",Rposition,Rf,FOOT_STEP_LENGTH
                else:
                        ## Rsupport, compute next left foot
                        [FOOT_STEP_LENGTH, Lposition] = GetStepLengthLeftFoot(COM_project, nrml, Lposition)
                        Lf = COM_project[:,Lposition] + 0.5*FOOT_SPACING*nrml[:,Lposition]
                        dLf = der[:,Lposition]
                        Lsupport = True
                        leftFoot.append(Lf)
                        leftFootDer.append(dLf)
                        print "L:",Lposition,Lf,FOOT_STEP_LENGTH

        Lf = np.array(leftFoot)
        Rf = np.array(rightFoot)
        dLf = np.array(leftFootDer)
        dRf = np.array(rightFootDer)
        return [Lf, dLf, Rf, dRf]

def COM_compute_zig_zag_motion(COM_linear, env):
        del handles[:]
        ### in (cm)

        ##assume that we start with both feet spaced apart at start pos

        M = COM_linear.shape[1]
        ## project COM onto surface
        COM_project = np.zeros((3,M))
        COM_project[0:2,:]=COM_linear[0:2,:]

        [Lf, dLf, Rf, dRf] = GetFootPositionFromProjectedCOM( COM_project )

        handleL = visualizeFoot( env, Lf, dLf, COLOR_LEFT_FOOT)
        handles.append(handleL)
        handleR = visualizeFoot( env, Rf, dRf, COLOR_RIGHT_FOOT)
        handles.append(handleR)

        Rsteps = Rf.shape[0]
        Lsteps = Lf.shape[0]

        FFn = Lf.shape[0] + Rf.shape[0]
        COM = np.zeros((FFn,3))
        ctr=0
        for i in range(0,min(Lf.shape[0],Rf.shape[0])):
                COM[ctr,:] = Rf[i,:]
                ctr+=1
                COM[ctr,:] = Lf[i,:]
                ctr+=1

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
        [trajectory,tmp] = splprep(COM,k=3,s=0)
        COM_traj = np.array([splev(t,trajectory) for t in tvec]).T
        COM_traj[2,:] = COM_linear[2,:]

        ### interpolate more waypoints
        M = 1.5*M
        tvec = np.linspace(0,1,M)
        [trajectory,tmp] = splprep(COM_traj,k=3,s=0)
        COM_traj = np.array([splev(t,trajectory) for t in tvec]).T

        COM_traj_projected = np.zeros((3,M))
        COM_traj_projected[0:2,:] = COM_traj[0:2,:]

        handles.append(env.env.drawlinestrip(points=COM_linear.T,
                           linewidth=10,
                           colors=np.array((1.0,0.0,1.0,0.8))))
        handles.append(env.env.drawlinestrip(points=COM_traj.T,
                           linewidth=8,
                           colors=np.array((1.0,0.0,1.0,0.4))))
        #raw_input('Press <ENTER> to continue.')
        #######################################################################
        ### compute end effector trajectory of footsteps
        #######################################################################

        ### footpos should have same size as COM input
        footpos = np.hstack((Lf[0,:],Rf[0,:]))
        dfootpos = np.hstack((dLf[0,:],dRf[0,:]))

        ### COM_traj starts in R, then moves to L without movement
        Lposition = 0
        Rposition = 0
        COMposition = 0

        ### start (shift COM from right foot to left foot, so that L becomes SF)
        [istart, ictr] = findNextCOMsegment(COMposition, COM_traj_projected, Lf[Lposition,:])
        COMposition = ictr
        footpos = np.tile( footpos , (ictr-istart,1))
        dfootpos = np.tile( dfootpos , (ictr-istart,1))
        COM_adjusted = COM_traj[:,0:COMposition]

        Lsupport = True
        while COMposition < M:
                if Lsupport:
                        if Rposition+1 >= Rsteps:
                                ## no more right foot steps abort
                                break
                        #######################################################
                        ### left foot is on support, only move right foot
                        ### from Rf[Rposition,0] to Rf[Rposition+1,0]
                        #######################################################
                        [istart, ictr] = findNextCOMsegment(COMposition, COM_traj_projected, Rf[Rposition+1,:])
                        COMposition = ictr
                        Nwpts = ictr-istart

                        f1 = Rf[Rposition,:]
                        df1 = dRf[Rposition,:]
                        f2 = Rf[Rposition+1,:]
                        df2 = dRf[Rposition+1,:]

                        [frvec,dfrvec] = interpolateFoot(Nwpts, f1, df1, f2, df2)

                        fl = Lf[Lposition,:]
                        dfl = dLf[Lposition,:]
                        flvec = np.tile( fl, (frvec.shape[0],1) )
                        dflvec = np.tile( dfl, (dfrvec.shape[0],1) )

                        footpos = np.vstack((footpos,np.hstack((flvec, frvec))))
                        dfootpos = np.vstack((dfootpos,np.hstack((dflvec, dfrvec))))

                        #######################################################
                        ### fill up with lastfootpos
                        #######################################################
                        M_subpath1 = frvec.shape[0]
                        lastfootpos = np.zeros(3)
                        lastfootpos[0] = footpos[-1,0]
                        lastfootpos[1] = footpos[-1,1]
                        lastfootpos[2] = COM_traj[2,COMposition]
                        COM_subpath1 = np.tile(lastfootpos,(M_subpath1,1)).T
                        COM_subpath2 = COM_traj[:,istart:ictr]
                        COM_subpath = np.hstack((COM_subpath1,COM_subpath2))
                        COM_adjusted = np.hstack((COM_adjusted,COM_subpath))
                        ### during subpath2, fill up fixed footpos
                        fvec = np.tile( footpos[-1,:], (frvec.shape[0],1) )
                        footpos = np.vstack((footpos,fvec))
                        dfvec = np.tile( dfootpos[-1,:], (frvec.shape[0],1) )
                        dfootpos = np.vstack((dfootpos,dfvec))
                        #######################################################

                        Rposition = Rposition+1
                        Lsupport = False
                        #sys.exit(0)
                else:
                        if Lposition+1 >= Lsteps:
                                ## no more right foot steps abort
                                break
                        #######################################################
                        ### right foot is on support, only move left foot
                        ### from Lf[Lposition,0] to Lf[Lposition+1,0]
                        #######################################################
                        [istart, ictr] = findNextCOMsegment(COMposition, COM_traj_projected, Lf[Lposition+1,:])
                        COMposition = ictr
                        Nwpts = ictr-istart

                        print "intrplt:",istart,ictr
                        f1 = Lf[Lposition,:]
                        df1 = dLf[Lposition,:]
                        f2 = Lf[Lposition+1,:]
                        df2 = dLf[Lposition+1,:]

                        [flvec,dflvec] = interpolateFoot(Nwpts, f1, df1, f2, df2)

                        fr = Rf[Rposition,:]
                        dfr = dRf[Rposition,:]
                        frvec = np.tile( fr, (flvec.shape[0],1) )
                        dfrvec = np.tile( dfr, (dflvec.shape[0],1) )

                        footpos = np.vstack((footpos,np.hstack((flvec, frvec))))
                        dfootpos = np.vstack((dfootpos,np.hstack((dflvec, dfrvec))))

                        #######################################################
                        ### fill up with lastfootpos
                        #######################################################
                        M_subpath1 = flvec.shape[0]
                        lastfootpos = np.zeros(3)
                        lastfootpos[0] = footpos[-1,3]
                        lastfootpos[1] = footpos[-1,4]
                        lastfootpos[2] = COM_traj[2,COMposition]
                        COM_subpath1 = np.tile(lastfootpos,(M_subpath1,1)).T
                        COM_subpath2 = COM_traj[:,istart:ictr]
                        COM_subpath = np.hstack((COM_subpath1,COM_subpath2))
                        COM_adjusted = np.hstack((COM_adjusted,COM_subpath))
                        ### during subpath2, fill up fixed footpos
                        fvec = np.tile( footpos[-1,:], (flvec.shape[0],1) )
                        footpos = np.vstack((footpos,fvec))
                        dfvec = np.tile( dfootpos[-1,:], (flvec.shape[0],1) )
                        dfootpos = np.vstack((dfootpos,dfvec))
                        #######################################################
                        #sys.exit(0)

                        Lposition = Lposition+1
                        Lsupport = True


        lastfootpos = footpos[-1,:]
        dlastfootpos = dfootpos[-1,:]
        ### stack it up until finished
        while footpos.shape[0] < M:
                footpos = np.vstack((footpos,lastfootpos))
                dfootpos = np.vstack((dfootpos,dlastfootpos))

        COM_adjusted = np.hstack((COM_adjusted,COM_traj[:,COMposition:-1]))

        handles.append(env.env.drawlinestrip(points=footpos[:,0:3],
                           linewidth=8,
                           colors=COLOR_LEFT_FOOT))

        handles.append(env.env.drawlinestrip(points=footpos[:,3:6],
                           linewidth=8,
                           colors=COLOR_RIGHT_FOOT))

        #raw_input('Press <ENTER> to continue.')

        return [COM_adjusted,footpos,dfootpos]
        #return [COM_traj,footpos,dfootpos]

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

### compute signed angle between vector in xy plane and x-axis (yaw)
def GetSignedAngle(dF):
        ndF = dF/np.linalg.norm(dF)
        ex = np.array((1,0,0))
        angle = math.acos(np.dot(ndF,ex))
        sign = np.sign(np.cross(ex, ndF)[2])
        if sign >= 0:
                ## left side rotation [0,pi]
                angle = angle
        else:
                angle = -angle
        return angle

def HTfromPosDer(F, dF):
        H = np.eye(4)
        H[0:3,3] = F
        angle = GetSignedAngle(dF)
        H[0:3,0:3] = Rz(angle)
        return H

def createTransformFromPosDer( fpos, dfpos):
        HL = HTfromPosDer( fpos[0:3], dfpos[0:3])
        HR = HTfromPosDer( fpos[3:6], dfpos[3:6])
        return [HL,HR]


def GIK_from_COM(COM_path, q_original, robot, env, recompute=False, DEBUG=False):
        q_gik_fname = 'tmp/q_gik.numpy'
        COM_gik_fname = 'tmp/COM_gik.numpy'

        N = q_original.shape[0]
        M = q_original.shape[1]

        q_gik = np.zeros((N,M))
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
                                                        returnclosest=True,
                                                        #checkcollisionlink=['l_foot','r_foot'],
                                                        #obstacles=obstacle_list,
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

def GIK_from_COM_and_FOOTPOS(COM_path, footpos, dfootpos, robot, env, recompute=False, DEBUG=False):
        q_gik_fname = 'tmp/q_gik2.numpy'
        COM_gik_fname = 'tmp/COM_gik2.numpy'
        q_old = env.surrender_pos

        #N = q_original.shape[0]
        #M = q_original.shape[1]
        N = q_old.shape[0]
        M = COM_path.shape[1]

        q_gik = np.zeros((N,M))
        COM_gik = np.zeros((3,M))
        if not os.path.isfile(q_gik_fname+'.npy') or recompute:
                i = 0
                with env.env:
                        cbirrt = CBiRRT(env.env, env.robot_name)
                        while i < M:
                                print "------------------------------------------------------------------"
                                print "------------ WAYPOINT",i,"/",M,": GIK  ---------------------------"
                                print "------------------------------------------------------------------"
                                try:

                                        [qlimL,qlimU]=robot.GetActiveDOFLimits()
                                        env.EnforceLimits(q_old, qlimL, qlimU)

                                        robot.SetActiveDOFValues(q_old)
                                        #robot.SetTransform(np.array([[1,0,0,COM_path[0,i]],[0,1,0,COM_path[1,i]],[0,0,1,COM_path[2,i]],[0,0,0,1]]))

                                        left_leg_tf_ori = robot.GetManipulator('l_leg').GetTransform()
                                        right_leg_tf_ori = robot.GetManipulator('r_leg').GetTransform()
                                        #left_arm_tf = robot.GetManipulator('l_arm').GetTransform()
                                        #right_arm_tf = robot.GetManipulator('r_arm').GetTransform()

                                        support_list = []
                                        maniptm_list = []

                                        [left_leg_tf, right_leg_tf] = createTransformFromPosDer( footpos[i,:], dfootpos[i,:] )

                                        zl = left_leg_tf[2,3]
                                        zr = right_leg_tf[2,3]

                                        lfoot_contact = False
                                        rfoot_contact = False

                                        maniptm_list = [('l_leg',left_leg_tf),('r_leg',right_leg_tf)]
                                        if zl <= Z_FOOT_CONTACT:
                                                support_list.append(('l_leg',SURFACE_FRICTION))
                                                #maniptm_list.append(('l_leg',left_leg_tf))
                                                lfoot_contact = True
                                        if zr <= Z_FOOT_CONTACT:
                                                support_list.append(('r_leg',SURFACE_FRICTION))
                                                #maniptm_list.append(('r_leg',right_leg_tf))
                                                rfoot_contact = True
                                        
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


