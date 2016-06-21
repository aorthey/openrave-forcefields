import numpy as np
import abc

class DeformationModule():
        __metaclass__ = abc.ABCMeta
        DeformInfo = {}
        COLLISION_ENABLED = False
        lambda_coeff = None

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
                                print "## [",self.get_name(),"] update: ok, lambda:",lambda_coeff
                        else:
                                print "## [",self.get_name(),"] >>>>> collision <<<<<"
                                ### make lambda smaller ?
                else:
                        dU += dUtmp

                return dU

        @abc.abstractmethod 
        def get_gradient(self, lambda_coeff):
                pass

        @abc.abstractmethod 
        def get_name(self):
                pass

        def SmoothVector(self, traj, Ncritical, W, smoothing_factor=None):
                if smoothing_factor is None:
                        smoothing_factor = self.DeformInfo['smoothing_factor']
                [Ndim, Nwaypoints] = traj.getWaypointDim(W)
                A = np.zeros(Nwaypoints)

                assert(Ncritical<Nwaypoints)

                i = Nwaypoints-1
                while i >= 0:
                        A[i] = self.avalue(Ncritical, i, smoothing_factor)
                        i -= 1
                return A

        def avalue(self, Ncritical, i, c):
                return np.exp(-((Ncritical-i)*(Ncritical-i))/(2*c*c))
