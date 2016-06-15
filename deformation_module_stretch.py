import abc
import sys
from deformation_module import *
from util import *

class DeformationModuleStretch(DeformationModule):

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
                critical_pt = self.DeformInfo['critical_pt']

                #dUtmp = np.zeros((Ndim,Nwaypoints))
                tangent = dWori[:,critical_pt]
                tangent /= np.linalg.norm(tangent)
                normal = np.dot(Rz(pi/2),np.array((tangent[0],tangent[1],1)))
                normal = np.array((normal[0], normal[1], 0.0, 0.0))

                ictr = 0
                epsilon = 1e-2
                dc = 0.0

                ### move backwards, identify the space of allowed points
                while dc < epsilon:
                        if critical_pt-ictr >= 0:
                                [dc,Wnext] = self.GetClosestDistanceLineToPoint( Wori[:,critical_pt], Wori[:,critical_pt-ictr], tangent )
                                print "critical_pt:",critical_pt,"wp:",critical_pt-ictr,"dist:",dc
                                ictr+=1
                        else: 
                                print "done"
                                break
                first_pt = critical_pt
                last_pt = critical_pt-ictr+1
                print "first pt:",first_pt,"last pt:",last_pt

                if ContainsReidemeisterTwist( tangent, Wori, first_pt, last_pt):

                        [idxW1,idxW2,idxW3] = self.IdentifyReidemeisterSubpaths( tangent, Wori, first_pt, last_pt)
                        Wdir = np.zeros((Ndim,Nwaypoints))
                        dUtmp = np.zeros((Ndim,Nwaypoints))

                        for i in idxW1:
                                Wdir[:,i] = -tangent
                        #for i in idxW2:
                                #Wdir[:,i] = 0 
                        for i in idxW3:
                                Wdir[:,i] = tangent+0.1*normal

                        Wdir[:,idxW3[0]] = normal

                        for i in range(0,Nwaypoints):
                                A = self.SmoothVector(traj,i,Wori)
                                #dUtmp[:,i] += np.dot(A,( lambda_coeff * Wdir.T))
                                dUtmp[:,i] += lambda_coeff * Wdir[:,i]

                else:
                        ### identify the three reidemeister segments, or create them

                        print "Creating NEW infinitesimal reidemeister twist"
                        dUtmp = self.InfinitesimalReidemeisterType1Twist( tangent, Wori, idxW2)

                        sys.exit(0)


                return dUtmp

                #sys.exit(0)

        def get_name(self):
                return "reidemeister stretch"

        def InfinitesimalReidemeisterType1Twist( self, tangent, W, idx):
                print len(idx)

        def ContainsReidemeisterTwist( self, tangent, W, first_pt, last_pt ):
                i = first_pt 
                dW = W[:,i]-W[:,i-1]
                direction = np.dot(dW,tangent)
                sec1_start = i
                while direction > 0:
                        if i <= last_pt:
                                return False
                        dW = W[:,i]-W[:,i-1]
                        direction = np.dot(dW,tangent)
                        i-=1

                return True

        def IdentifyReidemeisterSubpaths( self, tangent, W, first_pt, last_pt ):
                i = first_pt 
                print W.shape

                dW = W[:,i]-W[:,i-1]
                direction = np.dot(dW,tangent)
                sec1_start = i
                while direction > 0:
                        if i <= last_pt:
                                print "end section 1"
                                break
                        dW = W[:,i]-W[:,i-1]
                        direction = np.dot(dW,tangent)
                        i-=1
                sec1_end = i

                if i <= last_pt:
                        print "reidemeistering subpath"
                        print "create new reidemeister twist"

                        ### divide section 1 into three parts
                        N1 = sec1_start - sec1_end
                        Nd = int(N1/3)
                        sec1_end = N1-Nd
                        sec2_start = sec1_end
                        sec2_end = N1-2*Nd
                        sec3_start = sec2_end
                        sec3_end = 0
                else:
                        sec2_start = i
                        print "identify sections 2/3"
                        while direction <= 0:
                                if i < last_pt:
                                        print "end section 2"
                                        break
                                dW = W[:,i]-W[:,i-1]
                                direction = np.dot(dW,tangent)
                                i-=1
                        sec2_end = i

                print "section1:",sec1_start,"->",sec1_end
                print "section2:",sec2_start,"->",sec2_end
                print "section3:",sec3_start,"->",sec3_end
                idx1 = np.arange(sec1_end,sec1_start)
                idx2 = np.arange(sec2_end,sec2_start)
                idx3 = np.arange(sec3_end,sec3_start)
                return [idx1,idx2,idx3]

        def GetClosestDistanceLineToPoint( self, p1, p2, dp1 ):
                ### start at p1, move along dp1 until minimum to p2 reached
                dp1 /= np.linalg.norm(dp1)

                dold = 1e5
                dnew = 1e3
                gamma = 0.0
                gamma_step = 1e-3
                while dnew < dold:
                        dold = dnew
                        pn = p1 - gamma*dp1
                        dnew = np.linalg.norm(p2-pn)
                        gamma += gamma_step

                pn = p1 - (gamma-gamma_step)*dp1
                return [np.linalg.norm(pn-p2),pn]

                

