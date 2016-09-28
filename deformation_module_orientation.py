import abc
from deformation_module import *
from util import *

class DeformationModuleOrientation(DeformationModule):

        def get_gradient(self, lambda_coeff):

                traj = self.DeformInfo['traj']
                env = self.DeformInfo['env']
                eta = self.DeformInfo['eta']
                Wori = self.DeformInfo['Wori']
                dWori = self.DeformInfo['dWori']
                FN = self.DeformInfo['FN']
                Ndim = self.DeformInfo['Ndim']
                Nwaypoints = self.DeformInfo['Nwaypoints']
                F = self.DeformInfo['F']
                dpmin = self.DeformInfo['dpmin']
                dpmax = self.DeformInfo['dpmax']
                amin = self.DeformInfo['amin']
                amax = self.DeformInfo['amax']

                ds = traj.DISCRETIZATION_TIME_STEP
                Tdir = np.zeros((1,Nwaypoints))
                dUtmp = np.zeros((Ndim,Nwaypoints))

                for i in range(0,Nwaypoints-1):
                        p = Wori[:,i]
                        theta = p[3]
                        dp = dWori[:,i]
                        pnext = Wori[:,i+1]
                        smax = dpmax[:,i]/2

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


                                pori = np.zeros((Ndim))
                                dline = np.zeros((Ndim))
                                pori[0:3] = p[0:3]

                                dline[0:3] = np.dot(Rz(pi/2),(dWori[:,i])[0:3])
                                theta_step = 1e-2

                                pori[3] = p[3]+theta_step
                                [R,atmp,a2tmp] = traj.getControlMatrix(pori)
                                R = R[:,:,0]
                                #qnext1 = pori + dt*smax*dp + dt2*F[:,i]
                                qproj_f1 = np.dot(np.dot(Rz(pori[3]),ex),FN[0:3,i])
                                vol1 = self.GetReachableSetProjectedVolume( dt, ds, qnext, dline, pori, smax, dp, F[:,i], R,amin, amax)

                                pori[3] = p[3]-theta_step
                                [R,atmp,a2tmp] = traj.getControlMatrix(pori)
                                R = R[:,:,0]
                                #qnext2 = pori + dt*smax*dp + dt2*F[:,i]
                                qproj_f2 = np.dot(np.dot(Rz(pori[3]),ex),FN[0:3,i])
                                vol2 = self.GetReachableSetProjectedVolume( dt, ds, qnext, dline, pori, smax, dp, F[:,i], R,amin, amax)

                                pori[3] = p[3]
                                [R,atmp,a2tmp] = traj.getControlMatrix(pori)
                                R = R[:,:,0]
                                vol3 = self.GetReachableSetProjectedVolume( dt, ds, qnext, dline, pori, smax, dp, F[:,i], R,amin, amax)

                                if vol1 <= vol3 and vol2 <= vol3:
                                        Tdir[:,i]=0
                                else:
                                        if vol1 > vol3 and vol2 > vol3:
                                                ## both directions are good, so
                                                ## go against force field to
                                                ## increase reachability
                                                if qproj_f1 > qproj_f2:
                                                        ## move to v2
                                                        Tdir[:,i] = -1
                                                else:
                                                        Tdir[:,i] = 1

                                        elif vol1 > vol2:
                                                ##move into + direction
                                                Tdir[:,i] = 1
                                        else:
                                                Tdir[:,i] = -1

                for i in range(0,Nwaypoints):
                        A = self.SmoothVector(traj,i,Wori)
                        dUtmp[3,i] = np.dot(A, (lambda_coeff*Tdir.T))

                return dUtmp

        def get_name(self):
                return "orientation alignment"

        def GetReachableSetVerticesSE2(self, dt, p, s, dp, F, R, amin, amax):
                Ndim = p.shape[0]
                dt2 = dt*dt*0.5
                a = np.zeros((amin.shape))

                rs_vertices = np.zeros((Ndim,4))

                ## create ordered set (clockwise)
                a[0] = amin[0]
                a[1] = amin[1]
                rs_vertices[:,0] = p + dt*s*dp + dt2*F + dt2*np.dot(R,a)
                a[0] = amin[0]
                a[1] = amax[1]
                rs_vertices[:,1] = p + dt*s*dp + dt2*F + dt2*np.dot(R,a)
                a[0] = amax[0]
                a[1] = amax[1]
                rs_vertices[:,2] = p + dt*s*dp + dt2*F + dt2*np.dot(R,a)
                a[0] = amax[0]
                a[1] = amin[1]
                rs_vertices[:,3] = p + dt*s*dp + dt2*F + dt2*np.dot(R,a)

                return rs_vertices

        def GetReachableSetProjectedVolume( self, dt, ds, pline, dline, p, smax, dp, F, R, amin, amax):

                rs_vertices = self.GetReachableSetVerticesSE2( dt, p, smax, dp, F, R, amin, amax)
                Nvertices = rs_vertices.shape[1]

                ##project onto ds-ball around p

                rs_vertices_proj = np.zeros(rs_vertices.shape)

                dline = dline/np.linalg.norm(dline)

                dp = np.zeros(Nvertices)

                for i in range(0,Nvertices):
                        pvertex = rs_vertices[:,i]-pline
                        dp[i] = np.dot(pvertex[0:2], dline[0:2])

                mindp = np.min(dp)
                maxdp = np.max(dp)
                return maxdp - mindp
