import abc
import pylab as plt
from pylab import *
import sys
from deformation_module import *
from trajectory import *
from util import *

class DeformationModuleStretch(DeformationModule):

        DEBUG = False
        handler = []
        def get_gradient(self, lambda_coeff):

                traj = self.DeformInfo['traj']
                env = self.DeformInfo['env']
                Wori = self.DeformInfo['Wori']
                dWori = self.DeformInfo['dWori']
                Ndim = self.DeformInfo['Ndim']
                Nwaypoints = self.DeformInfo['Nwaypoints']
                critical_pt = self.DeformInfo['critical_pt']

                first_pt = critical_pt-1
                last_pt = 0
                if self.DEBUG:
                        print "first pt:",first_pt,"last pt:",last_pt

                tangent = dWori[:,first_pt]
                tangent /= np.linalg.norm(tangent)
                normal = np.dot(Rz(pi/2),np.array((tangent[0],tangent[1],1)))
                normal = np.array((normal[0], normal[1], 0.0, 0.0))

                ictr = 0
                epsilon = 1e-2
                dc = 0.0


                dUtmp = np.zeros((Ndim,Nwaypoints))
                #[cp_semin, cp_semax] = traj.GetSpeedIntervalAtCriticalPoint(env, Wori, dWori, critical_pt)
                #if cp_semax > 0.5:
                        #return dUtmp

                if self.ContainsReidemeisterTwist( tangent, Wori, first_pt, last_pt):

                        [idxW1,idxW2,idxW3] = self.IdentifyReidemeisterSubpaths(
                                        tangent, Wori, dWori, first_pt, last_pt)

                        sf = 30.0
                        ########################### idx1
                        Wdir = np.zeros((Ndim,Nwaypoints))
                        idxW1
                        idxW1_start = idxW1[0]
                        idxW1_end = idxW1[-1]

                        for idx in idxW1:
                                Wdir[:,idx] = -tangent-0.1*normal

                        dUtmp1 = np.zeros((Ndim,Nwaypoints))
                        for i in range(0,Nwaypoints):
                                A = self.SmoothVectorStretch(traj,critical_pt,i,Wori,smoothing_factor=sf)
                                dUtmp1[:,i] += np.dot(A,( lambda_coeff * Wdir.T))

                        dUtmp += dUtmp1
                        #COLLISION_ENABLED=False
                        #Wnext = Wori + dUtmp1
                        #if traj.IsInCollision(env, Wnext):
                        #        print "idx1 collision"
                        #else:
                        #        dUtmp += dUtmp1

                        ########################### idx2
                        Wdir = np.zeros((Ndim,Nwaypoints))
                        for idx in idxW2:
                                Wdir[:,idx] = -tangent + 0.3*normal

                        dUtmp2 = np.zeros((Ndim,Nwaypoints))
                        for i in range(0,Nwaypoints):
                                A = self.SmoothVectorStretch(traj,critical_pt,i,Wori,smoothing_factor=sf)
                                dUtmp2[:,i] += np.dot(A,( lambda_coeff * Wdir.T))
                        dUtmp += dUtmp2
                        #Wnext = Wori + dUtmp2
                        #if traj.IsInCollision(env, Wnext):
                        #        print "idx2 collision"
                        #else:
                        #        dUtmp += dUtmp2
                        ########################### idx3
                        Wdir = np.zeros((Ndim,Nwaypoints))
                        for idx in idxW3:
                                Wdir[:,idx] = tangent - normal
                        dUtmp3 = np.zeros((Ndim,Nwaypoints))

                        for i in range(0,Nwaypoints):
                                A = self.SmoothVectorStretch(traj,critical_pt,i,Wori,smoothing_factor=sf)
                                dUtmp3[:,i] += np.dot(A,( lambda_coeff * Wdir.T))

                        dUtmp += dUtmp3
                        #Wnext = Wori + dUtmp3
                        #if traj.IsInCollision(env, Wnext):
                        #        print "idx3 collision"
                        #else:
                                #dUtmp += dUtmp3
                        ###########################
                        #sys.exit(0)


                else:
                        ### identify the three reidemeister segments, or create them
                        print "Creating NEW infinitesimal reidemeister twist"
                        dUtmp += self.InfinitesimalReidemeisterType1Twist(
                                        traj,
                                        tangent, Wori, first_pt, last_pt,
                                        traj.DISCRETIZATION_TIME_STEP)

                return dUtmp

        def get_name(self):
                return "reidemeister stretch"


        def ScaleFactorBall(self, i, Nprev, ds_ball):
                t=(i+1)/(Nprev-1.0)
                #d = ds_ball * np.sqrt(t)
                #d = ds_ball * np.sqrt((i+1)/(Nprev-1.0))
                #d = ds_ball*(-(1-t)**2+1.0)
                d = ds_ball*(-1.2*(1-t)**2+1.0)
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
                if self.DEBUG:
                        print "Found reidemeister twist at",i
                return True

        def IdentifyReidemeisterSubpaths( self, tangent, W, dW, first_pt, last_pt ):
                i = first_pt 

                #dW = W[:,i]-W[:,i-1]
                direction = np.dot(dW[:,i],tangent)
                sec1_start = i

                while direction > 0:
                        if i <= last_pt:
                                break
                        direction = np.dot(dW[:,i],tangent)
                        i-=1

                sec1_end = i

                if i <= last_pt:
                        print "reidemeister subpath:",
                        print "only one subpath exist"
                        return [[0],[0],[0]]
                        #sys.exit(0)
                else:
                        sec2_start = i
                        epsilon = 0.01
                        while direction <= epsilon:
                                if i <= last_pt:
                                        break
                                direction = np.dot(dW[:,i],tangent)
                                i-=1
                        sec2_end = i
                        sec3_start = sec2_end

                        while direction > 0:
                                if i <= last_pt:
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

                ### for debug drawing only
                offset_old_line = 0.03
                offset_new_line = -0.03
                ### clip the first Mclip percent
                N = first_pt - last_pt
                Nclip = int(Mclip*N)

                N = N - 2*Nclip
                ds_ball = ds_ball*int(N/4)

                #N = int(0.66*(first_pt - last_pt))
                if N < 3:
                        print "ERROR: cannot twist path with only",N,"samples"
                        return dUtmp
                else:
                        if N%2 == 0:
                                ## even
                                Nmid = N/2
                        else:
                                ## odd
                                Nmid = (N-1)/2+1
                        Nmid += Nclip + last_pt
                        Wmid = W[:,Nmid-1] + 0.5*(W[:,Nmid] - W[:,Nmid-1])

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

                        ## left update
                        Wupdate[:,last_pt+Nclip+i] = Wcirc[:,i] - W[:,last_pt+Nclip+i]
                        Wupdate[:,Nmid+i] = Wcirc[:,Nprev+i] - W[:,Nmid+i]

                #print Wupdate

                dUtmp = Wupdate

                if self.DEBUG:

                        first_pt = first_pt - int(Nclip*0.4)
                        Wnext = W + dUtmp
                        #Wnext = self.InsertMissingWaypoints(Wnext, traj.DISCRETIZATION_TIME_STEP)

                        tnext = Trajectory(Wnext)
                        np.savetxt("W",Wnext[0:2,:])
                        [Wnext2,tmp,tmpp] = tnext.get_waypoints_second_order(N=1000)

                        first_pt2=0
                        d = 100
                        while d > traj.DISCRETIZATION_TIME_STEP and first_pt2 < Wnext2.shape[1]-1:
                                d = np.linalg.norm(W[0,first_pt]-Wnext2[0,first_pt2])
                                first_pt2+=1

                        #last_pt = last_pt + int(Nclip*0.5)
                        print first_pt,first_pt2
                        fs = 22
                        ## original line
                        fig=figure(facecolor='white')
                        plt.plot(W[0,last_pt:first_pt],W[1,last_pt:first_pt]+offset_old_line,'-og',linewidth=9,markersize=5)
                        plt.plot(W[0,last_pt:first_pt],W[1,last_pt:first_pt]+offset_old_line,'ok',markersize=5)
                        plt.plot(Wmid[0],Wmid[1]+offset_old_line,'ok',markersize=12)

                        plt.plot(Wnext[0,last_pt:first_pt],Wnext[1,last_pt:first_pt]+offset_new_line,'-or',linewidth=9,markersize=5)
                        plt.plot(Wnext[0,last_pt:first_pt],Wnext[1,last_pt:first_pt]+offset_new_line,'ok',markersize=5)
                        plt.plot(Wmid[0],Wmid[1]+offset_new_line,'ok',markersize=12)

                        plt.plot(Wnext2[0,last_pt:first_pt2],Wnext2[1,last_pt:first_pt2]+2*offset_new_line,'-og',linewidth=4,markersize=5)
                        plt.plot(Wnext2[0,last_pt:first_pt2],Wnext2[1,last_pt:first_pt2]+2*offset_new_line,'ok',markersize=7)
                        plt.plot(Wmid[0],Wmid[1]+2*offset_new_line,'ok',markersize=12)

                        for i in range(last_pt,first_pt):
                                W0 = Wnext[0,i]
                                W1 = Wnext[1,i]+offset_new_line
                                I0 = W[0,i]
                                I1 = W[1,i]+offset_old_line
                                plt.plot([W0,I0],[W1,I1],'-',color=np.array((0.4,0.4,0.4,1)),markersize=10)

                        #plt.plot(W[0,last_pt+Nclip:first_pt-Nclip],W[1,last_pt+Nclip:first_pt-Nclip],'-or')

                        circle = plt.Circle((Wmid[0],Wmid[1]+offset_new_line),ds_ball,color='r',fill=False)
                        plt.gca().add_artist(circle)
                        #circle = plt.Circle((Wmid[0],Wmid[1]+offset_old_line),ds_ball,color='r',fill=False)
                        #plt.gca().add_artist(circle)
                        plt.xlabel('X',fontsize=fs)
                        plt.ylabel('Y',fontsize=fs)
                        #plt.axis('equal')
                        ax = fig.gca()
                        for tick in ax.xaxis.get_major_ticks():
                                tick.label.set_fontsize(fs) 
                        for tick in ax.yaxis.get_major_ticks():
                                tick.label.set_fontsize(fs) 
                        ax.xaxis.labelpad = 20
                        ax.yaxis.labelpad = 20
                        plt.show()

                return dUtmp

        def SmoothVectorStretch(self, traj, Ncritical, Ncur, W, smoothing_factor=None):
                if smoothing_factor is None:
                        smoothing_factor = self.DeformInfo['smoothing_factor']

                [Ndim, Nwaypoints] = traj.getWaypointDim(W)
                A = np.zeros(Nwaypoints)

                assert(Ncritical<Nwaypoints)

                i = Nwaypoints-1
                while i >= 0:
                        if i < Ncritical:
                                A[i] = self.avalue(Ncur, i, smoothing_factor)
                        else:
                                A[i]= 0.0
                        i -= 1
                return A
