import abc
import time
import numpy as np
from scipy.interpolate import interp1d,splev,splrep,splprep
from scipy.misc import derivative

class Trajectory():
        __metaclass__ = abc.ABCMeta

        traj = []
        handle = []
        ptsize = 0.05
        trajectory_color = np.array((0.7,0.2,0.7,0.9))

        def __init__(self, trajectory):
                self.traj = trajectory 

        def info(self):
                print "#### TRAJECTORY CLASS ######"
                print "LENGTH: ", self.get_length()
                print "START : ", self.evaluate_at(0)
                print "GOAL  : ", self.evaluate_at(1)

        @abc.abstractmethod 
        def evaluate_at(self, t):
                pass

        @abc.abstractmethod 
        def get_length(self):
                pass

        @abc.abstractmethod 
        def get_waypoints(self, N = 100):
                pass

        @classmethod
        def from_ravetraj(cls, ravetraj):
                pass

        @classmethod
        def from_waypoints(cls, W):
                pass

        def draw(self, env, keep_handle=True):
                L = self.get_length()
                dt = 0.05
                Nwaypoints = int(L/dt)

                [W,dW] = self.get_waypoints(N=Nwaypoints)

                tmp_handle = []
                for i in range(0,Nwaypoints):
                        pt = np.array(((W[0,i],W[1,i],0.1)))
                        tmp_handle.append(env.env.plot3(points=pt,
                                           pointsize=self.ptsize,
                                           colors=self.trajectory_color,
                                           drawstyle=1))
                if keep_handle:
                        self.handle = tmp_handle
                else:
                        return tmp_handle

        def draw_delete(self, env):
                self.handle = []
