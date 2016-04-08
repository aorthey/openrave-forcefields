import abc
import sys
import time
import numpy as np
from trajectory import *
import copy

DEFORM_NONE = 0
DEFORM_OK = 1
DEFORM_COLLISION = 2
DEFORM_NOPROGRESS = 3

class Deformation():
        DEBUG = 0
        __metaclass__ = abc.ABCMeta

        ## number of steps of a linear homotopy deformation between two given
        ## trajectories
        DEFORMATION_STEPS = 10

        env = []
        traj_ori = []
        traj_deformed = []
        traj_current = []

        traj_display = []
        handle = []
        forcehandle = []


        def __init__(self, trajectory, environment):
                self.env = environment

                ## traj_ori : original trajectory before any calls
                ## traj_current : the current deformation
                ## traj_deformed : the working copy of current for the next iteration
                ## traj_display : the trajectory which we display at the moment
                self.traj_ori = trajectory 
                self.traj_current = copy.copy(trajectory)
                self.traj_deformed = copy.copy(trajectory)
                self.traj_display = copy.copy(trajectory)

        @abc.abstractmethod 
        def deform_onestep(self):
                pass


        def deform(self, N_iter = 1):
                computeNewCriticalPoint = True
                for i in range(0,N_iter):
                        res = self.deform_onestep(computeNewCriticalPoint)
                        self.traj_current = copy.copy(self.traj_deformed)
                        self.draw_deformation() 
                        if res == DEFORM_NONE:
                                print "DEFORM_NONE"
                                return False
                        elif res== DEFORM_OK:
                                print "DEFORM_OK"
                                computeNewCriticalPoint=True
                        elif res== DEFORM_COLLISION:
                                print "DEFORM_COLLISION"
                                computeNewCriticalPoint=True
                        elif res== DEFORM_NOPROGRESS:
                                print "DEFORM_NOPROGRESS"
                                computeNewCriticalPoint=True

                return True

        def GetForcesAtWaypoints(self, W):
                Ndim = W.shape[0]
                Nwaypoints = W.shape[1]
                F = np.zeros((Ndim,Nwaypoints))
                for i in range(0,Nwaypoints):
                        pt = np.array(((W[0,i],W[1,i],-0.1,0.001)))
                        F[0:3,i] = self.env.GetForceAtX(pt)
                return F

        def draw_forces_at_waypoints(self, W, F):
                for i in range(0,W.shape[1]):
                        df = np.linalg.norm(F[:,i])
                        if df > 0.3:
                                nF = 0.3*F[:,i]/df
                        else:
                                nF = F[:,i]
                        W1 = np.array((W[0,i],W[1,i],W[2,i]))
                        W2 = np.array((W[0,i]+nF[0],W[1,i]+nF[1],W[2,i]+nF[2]))
                        P = np.array(((W1),(W2)))
                        ls = self.traj_current.linsize
                        h=self.env.env.drawlinestrip(points=P,linewidth=ls,colors=np.array(((0.8,0.2,0.2,0.9))))
                        self.forcehandle.append(h)


        def draw_deformation(self):
                ## visualize linear homotopy between two trajectories

                [W0,dW] = self.traj_display.get_waypoints() 
                Nwaypoints = W0.shape[1]
                [W1,dW] = self.traj_current.get_waypoints(Nwaypoints) 

                self.forcehandle = []
                if np.linalg.norm(W0-W1)<0.001:
                        print "No deformation"
                        self.handle = self.traj_current.draw(self.env, keep_handle=False)
                        return

                #print self.traj_display.info()
                #print self.traj_current.info()

                t1 = time.time()
                sys.stdout.write("## DRAWING LINEAR HOMOTOPY DEFORMATION ...\n")

                tdraw = 0.0
                for i in range(0,self.DEFORMATION_STEPS):
                        k = float(i)/float(self.DEFORMATION_STEPS)
                        Wk = (1-k)*W0 + k*W1
                        if self.traj_current.IsInCollision(self.env, Wk):
                                break
                        else:
                                self.traj_current.new_from_waypoints(Wk)
                                #self.handle = self.traj_current.draw_nocritical(self.env, keep_handle=False)
                                ti1 = time.time()
                                self.handle = self.traj_current.draw(self.env, keep_handle=False)
                                ti2 = time.time()
                                #if i == 0:
                                        #raw_input('Press <ENTER> to draw deformation.')
                                #time.sleep(0.001)
                        tdraw += ti2-ti1
                sys.stdout.write("DONE\n")
                t2 = time.time()
                if self.DEBUG:
                        print "DEFORMATION TIME(s): ",t2-t1," (drawing:",tdraw,")"

                self.traj_display = copy.copy(self.traj_current)


        def draw_trajectory_original(self):
                self.traj_ori.draw(self.env)

        def draw_trajectory_deformed(self):
                self.traj_current.draw(self.env)

