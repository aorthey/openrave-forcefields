import abc
from deformation_module import *

class DeformationModuleProjectionReachableSet(DeformationModule):

        def get_gradient(self, lambda_coeff):

                traj = self.DeformInfo['traj']
                env = self.DeformInfo['env']
                eta = self.DeformInfo['eta']
                Wori = self.DeformInfo['Wori']
                dWori = self.DeformInfo['dWori']
                Ndim = self.DeformInfo['Ndim']
                Nwaypoints = self.DeformInfo['Nwaypoints']
                F = self.DeformInfo['F']
                dpmin = self.DeformInfo['dpmin']
                dpmax = self.DeformInfo['dpmax']

                Wmove = np.zeros((Ndim,Nwaypoints))
                dUtmp = np.zeros((Ndim,Nwaypoints))

                ds = traj.DISCRETIZATION_TIME_STEP
                #### compute the intersection of a forward simulation and the
                #### ball with radius ds

                for i in range(0,Nwaypoints-1):
                        p = Wori[:,i]
                        dp = dWori[:,i]
                        pnext = Wori[:,i+1]
                        smax = dpmax[:,i]/4
                        #smax = 0

                        if np.linalg.norm(F[:,i])>1e-3:
                                qnext = np.zeros((Ndim))
                                qnext = p
                                tstep = 1e-4
                                dt = 0.0

                                dold = 1e5
                                dnew = abs(np.linalg.norm(p-qnext) - ds)
                                while dnew < dold:
                                        dold = dnew
                                        dt += tstep
                                        dt2 = dt*dt/2
                                        qnext = p + dt*smax*dp + dt2*F[:,i]
                                        dnew = abs(np.linalg.norm(p-qnext) - ds)

                                dpq = np.dot(np.linalg.norm(qnext),np.linalg.norm(pnext))
                                if dpq < 0.1:
                                        print "WaRNING: dpq:",dpq
                                        sys.exit(0)

                                ## project onto orthogonal direction to path
                                dp = qnext - pnext
                                dp = dp - np.dot(dp,pnext-p)*(pnext-p)
                                Wmove[:,i] = dp

                for i in range(0,Nwaypoints):
                        A = self.SmoothVector(traj,i,Wori)
                        dUtmp[:,i] = np.dot(A, (lambda_coeff * Wmove.T))
                return dUtmp

        def get_name(self):
                return "reachable-set projection"

        def SmoothVector(self, traj, Ncritical, W):
                [Ndim, Nwaypoints] = traj.getWaypointDim(W)
                A = np.zeros(Nwaypoints)

                assert(Ncritical<Nwaypoints)
                i = Nwaypoints-1
                while i > 0:
                        if i>Ncritical:
                                #A[i] = avalue(Ncritical, i, k*20.0)
                                A[i] = 0.0
                        else:
                                A[i] = (Ncritical-i+1)*traj.DISCRETIZATION_TIME_STEP
                                #A[i] = 2.0**(Ncritical-i)
                        i -= 1
                return A
