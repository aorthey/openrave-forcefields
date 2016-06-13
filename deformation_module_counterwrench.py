import abc
from deformation_module import *

class DeformationModuleCounterWrench(DeformationModule):

        def get_gradient(self, lambda_coeff):

                traj = self.DeformInfo['traj']
                env = self.DeformInfo['env']
                eta = self.DeformInfo['eta']
                Wori = self.DeformInfo['Wori']
                dWori = self.DeformInfo['dWori']
                Ndim = self.DeformInfo['Ndim']
                Nwaypoints = self.DeformInfo['Nwaypoints']
                F = self.DeformInfo['F']

                FNxy = self.getForceNormalComponent(F, dWori)
                FNtorque = self.getTorqueNormalComponent(F, Wori, dWori)

                dUtmp = np.zeros((Ndim,Nwaypoints))
                for i in range(0,Nwaypoints):
                        A = self.SmoothVector(traj,i,Wori)
                        dUtmp[:,i] += np.dot(A,( -lambda_coeff * FNxy.T))
                        #dUtmp[:,i] += np.dot(A,( -lambda_coeff * FNtorque.T))
                return dUtmp

        def get_name(self):
                return "counter-wrench"

        def getForceNormalComponent(self, F, dWori):
                Nwaypoints = dWori.shape[1]
                FN = np.zeros((F.shape))
                for i in range(0,Nwaypoints):
                        dWn = dWori[:,i]/np.linalg.norm(dWori[:,i])
                        FN[:,i] = F[:,i]-np.dot(F[:,i],dWn)*dWn
                return FN

        def getTorqueNormalComponent(self, F, Wori, dWori):
                Nwaypoints = Wori.shape[1]
                FNtorque = np.zeros((F.shape))
                for i in range(0,Nwaypoints):
                        FNtorque[3,i] = F[3,i]
                return FNtorque
