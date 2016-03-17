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

def Extractabc(abc):
    lista = [float(x) for x in abc[0].split()]
    listb = [float(x) for x in abc[1].split()]
    listc = [float(x) for x in abc[2].split()]
    n= len(lista)/6
    a = zeros((n,6))
    b = zeros((n,6))
    c = zeros((n,6))
    for i in range(n):
        a[i,:] = lista[i*6:i*6+6]
        b[i,:] = listb[i*6:i*6+6]
        c[i,:] = listc[i*6:i*6+6]
    return a, b, c

def getLengthWpt(W):
        Ndim = W.shape[0]
        Nwaypoints = W.shape[1]
        d = 0.0
        for i in range(0,Nwaypoints-1):
                d = d+ np.linalg.norm(W[:,i+1]-W[:,i])
        return d
def GetTrajectoryStringR(RdimIdx, W, dW, ddW):
        ### get trajectory string for any R subspaces
        Ndim = W.shape[0]
        Nwaypoints = W.shape[1]

        L = getLengthWpt(W[RdimIdx,:])
        #duration = 1.0/float(Nwaypoints-1)
        #trajectorystring = str(duration)
        Rdim = len(RdimIdx)

        durationVector = np.zeros((Nwaypoints-1))
        for i in range(0,Nwaypoints-1):
                ds = np.linalg.norm(W[RdimIdx,i+1]-W[RdimIdx,i])
                dv = np.linalg.norm(dW[RdimIdx,i])
                duration = ds/dv
                durationVector[i] = duration
                #duration = d/L
                if i==0:
                        trajectorystring = str(duration)
                else:
                        trajectorystring += "\n" + str(duration)

                trajectorystring += "\n" + str(Rdim)
                ## f(s) = a + b s + c s^2
                A = W[RdimIdx,i]
                B = dW[RdimIdx,i]
                C = ddW[RdimIdx,i]*0.5
                for j in range(Rdim):
                        idx = RdimIdx[j]
                        trajectorystring += "\n"
                        trajectorystring += string.join([str(A[idx]),str(B[idx]),str(C[idx])])
        return [trajectorystring,durationVector]

def GetTrajectoryStringSO3(SO3dimIdx, W, dW, ddW):
        ### get trajectory string for any SO3 subspaces
        Ndim = W.shape[0]
        Nwaypoints = W.shape[1]
        SO3dim = len(SO3dimIdx)
        assert(SO3dim==3)
        SO3Valid = SO3dimIdx[SO3dimIdx>=0]

        SO3Valid = []
        for j in range(SO3dim):
                idx = SO3dimIdx[j]
                if idx is not None:
                        SO3Valid.append(idx)
        SO3Valid = np.array(SO3Valid)

        L = getLengthWpt(W[SO3Valid,:])
        durationVector = np.zeros((Nwaypoints-1))

        for i in range(0,Nwaypoints-1):
                duration = np.linalg.norm(W[SO3Valid,i+1]-W[SO3Valid,i])
                durationVector[i] = duration
                if i==0:
                        trajectorystring = str(duration)
                else:
                        trajectorystring += "\n" + str(duration)
                trajectorystring += "\n" + str(3)
                ## f(s) = a + b s + c s^2

                for j in range(SO3dim):
                        trajectorystring += "\n"

                        idx = SO3dimIdx[j]
                        if idx is None:
                                trajectorystring += string.join([str(0.0),str(0.0),str(0.0)])
                        else:
                                A = W[idx,i]
                                B = dW[idx,i]
                                C = ddW[idx,i]*0.5
                                trajectorystring += string.join([str(A),str(B),str(C)])

        
        print L
        if L < 1e-5:
                return [None,None]
        print trajectorystring,durationVector
        return trajectorystring

def GetTrajectoryStringComplete(W, dW, ddW):
        ### get trajectory string for any R subspaces
        Ndim = W.shape[0]
        Nwaypoints = W.shape[1]
        L = getLengthWpt(W)
        durationVector = np.zeros((Nwaypoints-1))

        for i in range(0,Nwaypoints-1):
                duration = np.linalg.norm(W[:,i+1]-W[:,i])
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

        print trajectorystring
        return [trajectorystring, durationVector]


def GetABCFromSO3Manifold(SO3dimIdx, amin, amax, W, dW, ddW):
        Ndim = W.shape[0]
        Nwaypoints = W.shape[1]
        duration = 1.0/float(Nwaypoints-1)

        [SO3trajstr,durationVector] = GetTrajectoryStringSO3(SO3dimIdx, W, dW, ddW)

        if SO3trajstr is None:
                a = np.zeros((Nwaypoints,1))
                b = np.zeros((Nwaypoints,1))
                c = np.zeros((Nwaypoints,1))
                return [a,b,c]


        #SO3traj = Trajectory.PiecewisePolynomialTrajectory.FromString(SO3trajstr)
        #SO3traj.Plot(0.05)
        #plt.show()
        constraintsstr = str(duration)

        SO3dim = len(SO3dimIdx)
        assert(SO3dim==3)
        #for j in range(SO3dim):
        ### TODO: make generic amax[1] -> amax[SO3]
        constraintsstr += "\n" + ' '.join([str(amax[1]),str(0.0),str(0.0)]) 
        constraintsstr += "\n" + ' '.join([str(amin[1]),str(0.0),str(0.0)]) 
        #print constraintsstr
        #print "###########################"

        inertia = eye(3)
        for v in inertia:
            constraintsstr += "\n" + ' '.join([str(i) for i in v])

        abc = TOPPbindings.RunComputeSO3Constraints(SO3trajstr,constraintsstr)
        a,b,c = Extractabc(abc)
        return [a,b,c]

def getSpeedProfile( RdimIdx, SO3dimIdx, F, R, amin, amax, W, dW, ddW ):

        Ndim = W.shape[0]
        Nwaypoints = W.shape[1]

        ### ManifoldIDX: 0 -> R-subspace. 1-> SO(2) subspace. 2-> SO(3) subspace
        ### (ordered three consecutive numbers as yaw/pitch/roll

        Rdim = RdimIdx.shape[0]

        [Rtrajstr,durationVector] = GetTrajectoryStringR(RdimIdx, W, dW, ddW)
        L = np.sum(durationVector)

        Rtraj = Trajectory.PiecewisePolynomialTrajectory.FromString(Rtrajstr)

        print "FINAL POINT on piecewise C^2 traj:",Rtraj.Eval(L)
        print "###############"

        ### compute a,b,c
        qs = np.zeros((Rdim,Nwaypoints))
        qss = np.zeros((Rdim,Nwaypoints))
        for i in range(0,Nwaypoints):
                duration = np.sum(durationVector[0:i])
                qs[:,i] = Rtraj.Evald(duration)
                qss[:,i] = Rtraj.Evaldd(duration)
                print duration,Rtraj.Eval(duration),qs[:,i],qss[:,i]

        print "###############"
        subplot(3,1,1)
        Rtraj.Plot(0.001)
        subplot(3,1,2)
        Rtraj.Plotd(0.001)
        subplot(3,1,3)
        Rtraj.Plotdd(0.001)
        plt.show()
        #sys.exit(0)

        I = np.identity(Rdim)
        G = np.vstack((I,-I))

        ar = np.zeros((Nwaypoints, 2*Rdim))
        br = np.zeros((Nwaypoints, 2*Rdim))
        cr = np.zeros((Nwaypoints, 2*Rdim))

        for i in range(0,Nwaypoints):
                Rmax = np.maximum(np.dot(R[RdimIdx,:,i],amax),np.dot(R[RdimIdx,:,i],amin))
                Rmin = np.minimum(np.dot(R[RdimIdx,:,i],amax),np.dot(R[RdimIdx,:,i],amin))
                H1 = F[RdimIdx,i] - Rmax
                H2 = -F[RdimIdx,i] + Rmin
                #for j in RdimIdx:
                #        print "qvol[",j,"]=",np.linalg.norm(-H1[j]-H2[j])
                #        if H2[j] > -H1[j]:
                #                print H2[j],"<= q[",j,"]<=",-H1[j]
                #                sys.exit(1)

                cr[i,:] = np.hstack((H1,H2)).flatten()

        for i in range(0,Nwaypoints):
                ar[i,:] = np.dot(G,qs[:,i]).flatten()
                br[i,:] = np.dot(G,qss[:,i]).flatten()

        # Constraints
        t0 = time.time()

        [aso3,bso3,cso3] = GetABCFromSO3Manifold(SO3dimIdx, amin, amax, W, dW, ddW)

        a = np.hstack((ar,aso3))
        b = np.hstack((br,bso3))
        c = np.hstack((cr,cso3))
        traj0 = GetTrajectoryStringComplete(W, dW, ddW)
        vmax = 1000*np.ones(Ndim)

        ############################################################################
        ### DEBUG: only retime on R manifold subspace
        ############################################################################
        #a = list(ar)
        #b = list(br)
        #c = list(cr)
        #traj0 = Rtraj
        #vmax = 1000*np.ones(Rdim)
        ############################################################################

        #topp_inst = TOPP.QuadraticConstraints(traj0, duration, vmax, list(a), list(b), list(c))
        print "###########################"

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
                # Display results
                ion()
                #x.WriteProfilesList()
                #x.WriteSwitchPointsList()
                #profileslist = list(TOPPpy.ProfilesFromString(x.resprofilesliststring))
                #switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
                #TOPPpy.PlotProfiles(profileslist,switchpointslist,4)
                x.WriteResultTrajectory()

                print x.restrajectorystring
                traj1 = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
                print "Trajectory duration before TOPP: ", traj0.duration
                print "Trajectory duration after TOPP: ", traj1.duration
                P = []
                tstep = traj1.duration/Nwaypoints
                subplot(3,1,1)
                traj1.Plot(0.001)
                subplot(3,1,2)
                traj1.Plotd(0.001)
                subplot(3,1,3)
                traj1.Plotdd(0.001)
                plt.show()

                t = 0
                while t < traj1.duration:
                        P.append(np.linalg.norm(traj1.Evald(t)))
                        t = t + tstep
                print "profile computed"
                #mvcbobrow = profileslist.pop(0)
                #mvc = np.array(mvcbobrow[3])
                #print mvc
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
        amin = np.array((-5,-5))
        amax = np.array((5,5))

        R = testMVCgetControlMatrix(W)
        print R
        F = np.zeros((W.shape))
        RdimIdx = np.array((0,1))
        SO3dimIdx = np.array((2,None,None))
        getSpeedProfile(RdimIdx, SO3dimIdx, F, R, amin, amax, W, dW, ddW)

