import abc
from deformation_module import *

class DeformationModuleEndPointProjection(DeformationModule):

        def get_gradient(self, lambda_coeff):

                traj = self.DeformInfo['traj']
                env = self.DeformInfo['env']
                Wori = self.DeformInfo['Wori']
                Ndim = self.DeformInfo['Ndim']
                Nwaypoints = self.DeformInfo['Nwaypoints']
                dU = self.DeformInfo['dU']

                Wmove = np.zeros((Ndim,Nwaypoints))
                dUtmp = np.zeros((Ndim,Nwaypoints))

                dEnd = dU[:,-1]
                dStart = dU[:,0]

                AEnd = self.SmoothVector(traj, Nwaypoints-1, Wori)
                AStart = self.SmoothVector(traj, 0, Wori)
                for i in range(0,Nwaypoints):
                        dU[:,i] += -AEnd[i]*dEnd
                        dU[:,i] += -AStart[i]*dStart

                return dU

        def get_name(self):
                return "endpoint subspace projection"
