import abc
import pylab as plt
import sys
from deformation_module import *
from util import *

class DeformationModuleStretch(DeformationModule):

        DEBUG = True
        handler = []
        def get_gradient(self, lambda_coeff):

                traj = self.DeformInfo['traj']
                env = self.DeformInfo['env']
                Wori = self.DeformInfo['Wori']
                dWori = self.DeformInfo['dWori']
                Ndim = self.DeformInfo['Ndim']
                Nwaypoints = self.DeformInfo['Nwaypoints']
                critical_pt = self.DeformInfo['critical_pt']

                dUtmp = np.zeros((Ndim,Nwaypoints))
                tangent = dWori[:,critical_pt]
                tangent /= np.linalg.norm(tangent)
                normal = np.dot(Rz(pi/2),np.array((tangent[0],tangent[1],1)))
                normal = np.array((normal[0], normal[1], 0.0, 0.0))

                ictr = 0
                epsilon = 1e-2
                dc = 0.0

                ### move backwards, identify the space of allowed points
                #while dc < epsilon:
                #        if critical_pt-ictr >= 0:
                #                ### check if point is STLC
                #                [dc,Wnext] = self.GetClosestDistanceLineToPoint( Wori[:,critical_pt], Wori[:,critical_pt-ictr], tangent )
                #                print "critical_pt:",critical_pt,"wp:",critical_pt-ictr,"dist:",dc
                #                ictr+=1
                #        else: 
                #                print "done"
                #                break

                first_pt = critical_pt
                #last_pt = critical_pt-ictr+1
                last_pt = 0
                print "first pt:",first_pt,"last pt:",last_pt

                if self.ContainsReidemeisterTwist( tangent, Wori, first_pt, last_pt):

                        [idxW1,idxW2,idxW3] = self.IdentifyReidemeisterSubpaths(
                                        tangent, Wori, dWori, first_pt, last_pt)
                        Wdir = np.zeros((Ndim,Nwaypoints))

                        #### apply deform potential onto reidemeister twist

                        print idxW1,idxW2,idxW3
                        for idx in idxW1:
                                Wdir[:,idx] = -tangent

                        for idx in idxW2:
                                Wdir[:,idx] = 0.5*normal

                        for i in range(0,Nwaypoints):
                                A = self.SmoothVector(traj,i,Wori,smoothing_factor=30.0)
                                dUtmp[:,i] += np.dot(A,( lambda_coeff * Wdir.T))

                        #### COUNTERACT any deformation on the left part of the
                        #### twist
                        #for idx in idxW2:

                        #idx = idxW3[0]
                        #A = self.SmoothVector(traj,idx,Wori,smoothing_factor=15.0)
                        #for i in range(0,Nwaypoints):
                        #        dUtmp[:,i] += -A[i]*dUtmp[:,idx]

                        #Wdirtmp = np.zeros((Ndim))
                        #Wdirtmp[:] = dUtmp[:,idx2]

                        #A = self.SmoothVector(traj,idx2,Wori,smoothing_factor=15.0)
                        #for i in range(0,Nwaypoints):
                                #dUtmp[:,i] += -A[i]*Wdirtmp

                        Wnext = Wori + dUtmp
                        Ic = traj.GetFirstCollisionPointIdx(env, Wnext)

                        ### CHECK IF IN COLLISION
                        return dUtmp
                        #while Ic is not None:
                        #print "[Stretch]",Ic
                        #while Ic is not None:
                        #        dW = Wnext[:,Ic] - Wnext[:,Ic-1]
                        #        ## slide along contact
                        #        print "[Stretch]: sliding along contact, standby"
                        #        A = self.SmoothVector(traj,Ic,Wori,smoothing_factor=15.0)
                        #        for i in range(0,Nwaypoints):
                        #                dUtmp[:,i] += -A[i]*dW
                        #        Wnext = Wori + dUtmp
                        #        Ic = traj.GetFirstCollisionPointIdx(env, Wnext)


                else:
                        ### identify the three reidemeister segments, or create them
                        print "Creating NEW infinitesimal reidemeister twist"
                        dUtmp = self.InfinitesimalReidemeisterType1Twist(
                                        traj,
                                        tangent, Wori, first_pt, last_pt,
                                        traj.DISCRETIZATION_TIME_STEP)

                return dUtmp

        def get_name(self):
                return "reidemeister stretch"


        def ScaleFactorBall(self, i, Nprev, ds_ball):
                d = ds_ball * np.sqrt((i)/(Nprev-1.0))
                #d = ds_ball * np.sqrt((i+1)/(Nprev-1.0))
                return d

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
                print "Found reidemeister twist at",i
                return True

        def IdentifyReidemeisterSubpaths( self, tangent, W, dW, first_pt, last_pt ):
                i = first_pt 

                #dW = W[:,i]-W[:,i-1]
                direction = np.dot(dW[:,i],tangent)
                sec1_start = i

                while direction > 0:
                        if i <= last_pt:
                                print "end section 1"
                                break
                        #dW = W[:,i]-W[:,i-1]
                        direction = np.dot(dW[:,i],tangent)
                        i-=1

                sec1_end = i

                if i <= last_pt:
                        print "reidemeister subpath:",
                        print "only one subpath exist"
                        sys.exit(0)
                else:
                        sec2_start = i
                        print "identify section 2/3"
                        epsilon = 0.01
                        while direction <= epsilon:
                                if i <= last_pt:
                                        print "end section 2"
                                        break
                                #dW = W[:,i]-W[:,i-1]
                                direction = np.dot(dW[:,i],tangent)
                                print direction
                                i-=1
                        sec2_end = i
                        sec3_start = sec2_end

                        while direction > 0:
                                if i <= last_pt:
                                        print "end section 3"
                                        break
                                direction = np.dot(dW[:,i],tangent)
                                i-=1
                        sec3_end = i

                Wpos = W[:,last_pt:first_pt]

                Wdir = 0.2*dW[:,last_pt:first_pt]

                for j in range(0,Wdir.shape[1]):
                        Wdir[:,j]/=np.linalg.norm(Wdir[:,j])
                        Wdir[:,j]*=0.1

                if self.DEBUG:
                        plt.plot([Wpos[0,:],Wpos[0,:]+Wdir[0,:]],[Wpos[1,:],Wpos[1,:]+Wdir[1,:]],'-om')
                        plt.plot(Wpos[0,:],Wpos[1,:],'-k',linewidth=2)
                        #Ws1 = Wpos[0:2,sec1_start:sec1_end]
                        Ws1 = Wpos[0:2,sec1_end:sec1_start]
                        Ws2 = Wpos[0:2,sec2_end:sec2_start]
                        Ws3 = Wpos[0:2,sec3_end:sec3_start]
                        plt.plot(Ws1[0,:],Ws1[1,:],'-or',linewidth=6)
                        plt.plot(Ws2[0,:],Ws2[1,:],'-ob',linewidth=6)
                        plt.plot(Ws3[0,:],Ws3[1,:],'-og',linewidth=6)
                        #plt.plot([Wpos[0,-1],Wpos[0,-1]+0.1*tangent[0]],[Wpos[1,-1],Wpos[1,-1]+0.1*tangent[1]],'-ob',linewidth=6)

                        plt.show()

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

        def InfinitesimalReidemeisterType1Twist( self, traj, tangent, W, first_pt, last_pt, ds_ball):
                Ndim = W.shape[0]
                Nwaypoints = W.shape[1]
                dUtmp = np.zeros((Ndim,Nwaypoints))
                Mclip = 0.4
                ### clip the first Mclip percent
                N = first_pt - last_pt
                Nclip = int(Mclip*N)
                print N,Mclip,Nclip

                N = N - 2*Nclip
                ds_ball = ds_ball*int(N/2)

                #N = int(0.66*(first_pt - last_pt))
                if N < 3:
                        print "ERROR: cannot twist path with only",N,"samples"
                        return dUtmp
                else:
                        if N%2 == 0:
                                ## even
                                Nmid = N/2 + Nclip + last_pt
                                ## draw circle around mid point
                                Wmid = W[:,Nmid-1] + 0.5*(W[:,Nmid] - W[:,Nmid-1])
                        else:
                                ## odd
                                Nmid = (N-1)/2
                                Wmid = 0.5*(W[:,Nmid] - W[:,Nmid-1])

                        print "Nmid:",Nmid,"last:",last_pt,"first:",first_pt,"N:",N
                        print "Nlast:",last_pt+Nclip,"Nfirst:",first_pt-Nclip,"Nclip:",Nclip

                ## project points onto ds_ball-ball

                ## left side points
                #Nnext = W[:,Nmid:first_pt-Nclip].shape[1]
                Nprev = W[:,last_pt+Nclip:Nmid].shape[1]

                piclip = pi/32
                pistep = (pi/2-2*piclip)/float(Nprev-1.0)

                tangentSE2 = np.zeros(3)
                tangentSE2[0:2] = tangent[0:2]
                tangentSE2[2] = tangent[3]

                Wcirc = np.zeros((Ndim, 2*Nprev))
                Wupdate = np.zeros((W.shape))
                for i in range(0,Nprev):
                        picur = piclip+pistep*i
                        picur2 = pi/2 + piclip+pistep*i

                        Wcircdir = np.dot(Rz(picur),tangentSE2)
                        Wcircdir2 = np.dot(Rz(picur2),tangentSE2)

                        Wtmp = np.zeros((Ndim))
                        Wtmp[0:2] = Wcircdir[0:2]
                        Wtmp[3] = Wcircdir[2]

                        ds_scale = self.ScaleFactorBall(i, Nprev+1, ds_ball)
                        print ds_scale
                        Wcirc[:,i] = Wmid + ds_scale*Wtmp

                        Wtmp[0:2] = Wcircdir2[0:2]
                        Wtmp[3] = Wcircdir2[2]

                        ds_scale = self.ScaleFactorBall(Nprev-i-1, Nprev+1, ds_ball)
                        Wcirc[:,Nprev+i] = Wmid + ds_scale*Wtmp

                        #plt.plot(Wcirc[0,i],Wcirc[1,i],'-og',markersize=10)
                        #plt.plot(W[0,last_pt+Nclip+i],W[1,last_pt+Nclip+i],'-og',markersize=10)
                        L1 = Wcirc[:,i]
                        L2 = W[:,last_pt+Nclip+i]
                        plt.plot([L1[0],L2[0]],[L1[1],L2[1]],'-g',markersize=10)
                        L3 = Wcirc[:,Nprev+i]
                        L4 = W[:,Nmid+i]
                        plt.plot([L3[0],L4[0]],[L3[1],L4[1]],'-r',markersize=10)

                        ## left update
                        Wupdate[:,last_pt+Nclip+i] = Wcirc[:,i] - W[:,last_pt+Nclip+i]
                        Wupdate[:,Nmid+i] = Wcirc[:,Nprev+i] - W[:,Nmid+i]

                #print Wupdate

                dUtmp = 0.5*Wupdate

                Wnext = W + dUtmp

                if self.DEBUG:
                        ## original line
                        plt.plot(W[0,last_pt:first_pt],W[1,last_pt:first_pt]+0.01,'-og',linewidth=8)

                        ## infinitesimal twist line
                        for i in range(last_pt,first_pt):
                                Wnext[1,i]+=i*0.001
                        plt.plot(Wnext[0,last_pt:first_pt],Wnext[1,last_pt:first_pt]-0.03,'-or',linewidth=1,markersize=2)

                        #plt.plot(W[0,last_pt+Nclip:first_pt-Nclip],W[1,last_pt+Nclip:first_pt-Nclip],'-or')
                        plt.plot(Wmid[0],Wmid[1],'-og',markersize=10)

                        circle = plt.Circle((Wmid[0],Wmid[1]),ds_ball,color='r',fill=False)
                        plt.gca().add_artist(circle)
                        plt.axis('equal')
                        plt.show()
                print "INSIDE:",Wnext
                print dUtmp
                return dUtmp
