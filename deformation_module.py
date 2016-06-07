import numpy as np
import abc

class DeformationModule():
        __metaclass__ = abc.ABCMeta
        DeformInfo = {}
        COLLISION_ENABLED = True

        def __init__(self, DeformInfoIn):
                self.DeformInfo = DeformInfoIn

        def get_update(self, lambda_coeff):

                traj_deformed = self.DeformInfo['traj']
                env = self.DeformInfo['env']
                eta = self.DeformInfo['eta']
                Wori = self.DeformInfo['Wori']

                dUtmp = self.get_gradient(lambda_coeff)

                Wnext = Wori + eta*dUtmp

                dU = np.zeros((dUtmp.shape))
                if self.COLLISION_ENABLED:
                        if not traj_deformed.IsInCollision(env, Wnext):
                                dU += dUtmp
                        else:
                                print "## $> collision in module:",self.get_name()
                else:
                        dU += dUtmp

                return dU

        @abc.abstractmethod 
        def get_gradient(self, lambda_coeff):
                pass

        @abc.abstractmethod 
        def get_name(self):
                pass

        def avalue(self, Ncritical, i, c=10.0):
                return np.exp(-((Ncritical-i)*(Ncritical-i))/(2*c*c))

        def A1matrix(self, traj, Ncritical, W):
                [Ndim, Nwaypoints] = traj.getWaypointDim(W)
                A1 = np.zeros(Nwaypoints)

                assert(Ncritical<Nwaypoints)

                i = Nwaypoints
                while i > 0:
                        #if abs(i-Ncritical)<M and i<Nwaypoints:
                        if i<Nwaypoints-1:
                                #A1[i] = self.avalue(Ncritical, i,15.0)
                                A1[i] = self.avalue(Ncritical, i,50.0)
                        i -= 1
                return A1