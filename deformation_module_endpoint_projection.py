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

                ##specify how many points at the start/end are static
                Mstatic = int(0.05*Nwaypoints)

                print Mstatic
                dEnd = dU[:,-Mstatic]
                dStart = dU[:,Mstatic]

                sf = 10.0
                #AEnd = self.SmoothVector(traj, Nwaypoints-1, Wori, smoothing_factor=sf)
                #AStart = self.SmoothVector(traj, 0, Wori, smoothing_factor=sf)
                AEnd = self.SmoothVector(traj, Nwaypoints-Mstatic, Wori, smoothing_factor=sf)
                AStart = self.SmoothVector(traj, Mstatic, Wori, smoothing_factor=sf)

                for i in range(0,Nwaypoints):
                        dUtmp[:,i] += -AEnd[i]*dEnd
                        dUtmp[:,i] += -AStart[i]*dStart

                ## set all static points to zero
                for i in range(0,Mstatic):
                        dUtmp[:,i] = -dU[:,i]
                        k = Nwaypoints-i-1
                        dUtmp[:,k] = -dU[:,k]

                return dUtmp


        def get_name(self):
                return "endpoint subspace projection"
