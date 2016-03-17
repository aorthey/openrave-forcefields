import string,time
from pylab import *
from numpy import *
from openravepy import *
import TOPP
from TOPP import Utilities
from TOPP import TOPPbindings
from TOPP import TOPPpy
from TOPP import Trajectory
from TOPP import TOPPopenravepy

def GetTrajectoryString(W, dW, ddW):
        ### get trajectory string for any R subspaces
        Ndim = W.shape[0]
        Nwaypoints = W.shape[1]
        durationVector = np.zeros((Nwaypoints-1))

        for i in range(0,Nwaypoints-1):
                ds  = np.linalg.norm(W[:,i+1]-W[:,i])
                dv = np.linalg.norm(dW[:,i])
                duration = ds/dv
                durationVector[i] = duration
                if i==0:
                        trajectorystring = str(duration)
                else:
                        trajectorystring += "\n" + str(duration)

                trajectorystring += "\n" + str(Ndim)
                ## f(s) = a + b s + c s^2
                A = W[:,i]
                B = dW[:,i]
                C = ddW[:,i]*0.5
                for j in range(Ndim):
                        trajectorystring += "\n"
                        trajectorystring += string.join([str(A[j]),str(B[j]),str(C[j])])

        return [trajectorystring, durationVector]

def getSpeedProfileRManifold( F, R, amin, amax, W, dW, ddW ):
        DEBUG = 0

        Ndim = W.shape[0]
        Nwaypoints = W.shape[1]

        [trajstr,durationVector] = GetTrajectoryString(W, dW, ddW)

        L = np.sum(durationVector)

        traj0 = Trajectory.PiecewisePolynomialTrajectory.FromString(trajstr)

        if DEBUG:
                print "FINAL POINT on piecewise C^2 traj:",traj0.Eval(L)
                print "FINAL WAYPOINT                   :",W[:,-1]
                print "###############"

        ### compute a,b,c
        qs = np.zeros((Ndim,Nwaypoints))
        qss = np.zeros((Ndim,Nwaypoints))
        for i in range(0,Nwaypoints):
                duration = np.sum(durationVector[0:i])
                qs[:,i] = traj0.Evald(duration)
                qss[:,i] = traj0.Evaldd(duration)
                print duration,traj0.Eval(duration),qs[:,i],qss[:,i]

        if DEBUG:
                subplot(3,1,1)
                traj0.Plot(0.001)
                subplot(3,1,2)
                traj0.Plotd(0.001)
                subplot(3,1,3)
                traj0.Plotdd(0.001)
                plt.show()

        I = np.identity(Ndim)
        G = np.vstack((I,-I))

        a = np.zeros((Nwaypoints, 2*Ndim))
        b = np.zeros((Nwaypoints, 2*Ndim))
        c = np.zeros((Nwaypoints, 2*Ndim))

        for i in range(0,Nwaypoints):
                Rmax = np.maximum(np.dot(R[:,:,i],amax),np.dot(R[:,:,i],amin))
                Rmin = np.minimum(np.dot(R[:,:,i],amax),np.dot(R[:,:,i],amin))
                H1 = F[:,i] - Rmax
                H2 = -F[:,i] + Rmin
                for j in range(Ndim):
                        #print "qvol[",j,"]=",np.linalg.norm(-H1[j]-H2[j])
                        if H2[j] > -H1[j]:
                                print H2[j],"<= q[",j,"]<=",-H1[j]
                                sys.exit(1)

                c[i,:] = np.hstack((H1,H2)).flatten()

        for i in range(0,Nwaypoints):
                a[i,:] = np.dot(G,qs[:,i]).flatten()
                b[i,:] = np.dot(G,qss[:,i]).flatten()

        vmax = 1000*np.ones(Ndim)
        topp_inst = TOPP.QuadraticConstraints(traj0, durationVector[0], vmax, list(a), list(b), list(c))
        x = topp_inst.solver

        #x.integrationtimestep = 0.001
        #x.reparamtimestep = 0.001
        #x.extrareps = 10

        ret = x.RunComputeProfiles(0,0)
        x.ReparameterizeTrajectory()

        print "TOPP Output:",ret
        if(ret == 4):
                print ret," [ERROR TOPP: MVC hit zero]"

        if(ret == 1):
                x.ReparameterizeTrajectory()
                ion()
                x.WriteResultTrajectory()

                traj1 = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)

                print "Trajectory duration before TOPP: ", traj0.duration
                print "Trajectory duration after TOPP: ", traj1.duration
                if DEBUG:
                        subplot(3,1,1)
                        traj1.Plot(0.01)
                        subplot(3,1,2)
                        traj1.Plotd(0.01)
                        subplot(3,1,3)
                        traj1.Plotdd(0.01)
                        plt.show()

                t = 0
                tstep = traj1.duration/Nwaypoints
                P = []
                while t < traj1.duration:
                        P.append(np.linalg.norm(traj1.Evald(t)))
                        t = t + tstep
                return np.array(P)
        return None

def testMVCgetControlMatrix(W):
        Ndim = W.shape[0]
        Nwaypoints = W.shape[1]

        assert(Ndim == 3)
        R = np.zeros((Ndim,2,Nwaypoints))
        for i in range(0,Nwaypoints):
                t = W[2,i]
                R[0,:,i] = np.array((cos(t),0.0))
                R[1,:,i] = np.array((sin(t),0.0))
                R[2,:,i] = np.array((0.0,1.0))
        return R

from scipy.interpolate import interp1d,splev,splrep,splprep
from scipy.misc import derivative

def testMVCgetDerivWpt(W):
        Ndim = W.shape[0]
        Nwaypoints = W.shape[1]
        dW = np.zeros((W.shape))
        ddW = np.zeros((W.shape))

        traj,tmp = splprep(W,k=5,s=0.01)
        #L = getLengthWpt(W)
        d = 0.0
        for i in range(0,Nwaypoints-1):
                dW[:,i] = splev(d,traj,der=1)
                ddW[:,i] = splev(d,traj,der=2)
                #dW[:,i] = dW[:,i]/np.linalg.norm(dW[:,i])

                ds = np.linalg.norm(W[:,i+1]-W[:,i])
                dv = np.linalg.norm(dW[:,i])
                dt = ds/dv
                #ddW[:,i] = ddW[:,i]/np.linalg.norm(ddW[:,i])
                print d
                d = d + dt

        dW[:,Nwaypoints-1] = splev(d,traj,der=1)
        ddW[:,Nwaypoints-1] = splev(d,traj,der=2)

        return [dW,ddW]

if __name__ == '__main__':
        W = np.array((
                        (0.0,0.1,0.2,0.3,0.4,0.5),
                        (0.0,0.0,0.0,0.0,0.0,0.0),
                        (0.0,0.0,0.0,0.0,0.0,0.0)
                ))
        [dW,ddW] = testMVCgetDerivWpt(W)
        print W
        print dW
        print ddW
        print "######################"
        amin = np.array((-2,-5))
        amax = np.array((3,5))

        R = testMVCgetControlMatrix(W)
        print R
        F = np.zeros((W.shape))
        getSpeedProfileRManifold(F, R, amin, amax, W, dW, ddW)

