import abc
from deformation_module import *

class DeformationModuleCounterWrench(DeformationModule):
        SMOOTHING_FACTOR = 30.0

        def get_gradient(self, lambda_coeff):

                traj = self.DeformInfo['traj']
                env = self.DeformInfo['env']
                eta = self.DeformInfo['eta']
                Wori = self.DeformInfo['Wori']
                dWori = self.DeformInfo['dWori']
                Ndim = self.DeformInfo['Ndim']
                Nwaypoints = self.DeformInfo['Nwaypoints']
                F = self.DeformInfo['F']
                FN = self.DeformInfo['FN']

                dUtmp = np.zeros((Ndim,Nwaypoints))
                for i in range(0,Nwaypoints):
                        A = self.SmoothVector(traj,i,Wori,smoothing_factor=self.SMOOTHING_FACTOR)
                        dUtmp[:,i] += np.dot(A,( -lambda_coeff * FN.T))
                return dUtmp

        def get_name(self):
                return "counter-wrench"

        #def SmoothVector(self, traj, Ncritical, W, smoothing_factor=None):
        #        if smoothing_factor is None:
        #                smoothing_factor = self.DeformInfo['smoothing_factor']
        #        [Ndim, Nwaypoints] = traj.getWaypointDim(W)
        #        A = np.zeros(Nwaypoints)

        #        assert(Ncritical<Nwaypoints)

        #        i = Nwaypoints-1
        #        while i >= 0:
        #                if i > Ncritical:
        #                        A[i] = 0.0
        #                else:
        #                        A[i] = self.avalue(Ncritical, i, smoothing_factor)
        #                i -= 1
        #        return A

