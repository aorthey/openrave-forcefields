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

def GetTrajectoryStringR(RdimIdx, W, dW, ddW):
        ### get trajectory string for any R subspaces
        Ndim = W.shape[0]
        Nwaypoints = W.shape[1]
        duration = 1.0/float(Nwaypoints-1)
        trajectorystring = str(duration)
        Rdim = len(RdimIdx)

        for i in range(0,Nwaypoints-1):
                if i > 0:
                        trajectorystring += "\n" + str(duration)
                trajectorystring += "\n" + str(Rdim)
                ## f(s) = a + b s + c s^2
                A = W[:,i]
                B = dW[:,i]
                C = ddW[:,i]*0.5
                for j in range(Rdim):
                        idx = RdimIdx[j]
                        trajectorystring += "\n"
                        trajectorystring += string.join([str(A[idx]),str(B[idx]),str(C[idx])])
        return trajectorystring

def GetTrajectoryStringSO3(SO3dimIdx, W, dW, ddW):
        ### get trajectory string for any SO3 subspaces
        Ndim = W.shape[0]
        Nwaypoints = W.shape[1]
        duration = 1.0/float(Nwaypoints-1)
        trajectorystring = str(duration)
        SO3dim = len(SO3dimIdx)
        assert(SO3dim==3)

        for i in range(0,Nwaypoints-1):
                if i > 0:
                        trajectorystring += "\n" + str(duration)
                trajectorystring += "\n" + str(3)
                ## f(s) = a + b s + c s^2
                A = W[:,i]
                B = dW[:,i]
                C = ddW[:,i]*0.5

                for j in range(SO3dim):
                        trajectorystring += "\n"

                        idx = SO3dimIdx[j]
                        if idx is None:
                                trajectorystring += string.join([str(0.0),str(0.0),str(0.0)])
                        else:
                                trajectorystring += string.join([str(A[idx]),str(B[idx]),str(C[idx])])
        return trajectorystring

def GetTrajectoryStringComplete(W, dW, ddW):
        ### get trajectory string for any R subspaces
        Ndim = W.shape[0]
        Nwaypoints = W.shape[1]
        duration = 1.0/float(Nwaypoints-1)
        trajectorystring = str(duration)

        for i in range(0,Nwaypoints-1):
                if i > 0:
                        trajectorystring += "\n" + str(duration)
                trajectorystring += "\n" + str(Ndim)
                ## f(s) = a + b s + c s^2
                A = W[:,i]
                B = dW[:,i]
                C = ddW[:,i]*0.5
                for j in range(Ndim):
                        trajectorystring += "\n"
                        trajectorystring += string.join([str(A[j]),str(B[j]),str(C[j])])

        return trajectorystring


def GetABCFromSO3Manifold(SO3dimIdx, amin, amax, W, dW, ddW):
        Ndim = W.shape[0]
        Nwaypoints = W.shape[1]
        duration = 1.0/float(Nwaypoints-1)

        SO3trajstr = GetTrajectoryStringSO3(SO3dimIdx, W, dW, ddW)

        #SO3traj = Trajectory.PiecewisePolynomialTrajectory.FromString(SO3trajstr)
        #SO3traj.Plot(0.05)
        #plt.show()
        inertia = eye(3)
        constraintsstr = str(duration)

        SO3dim = len(SO3dimIdx)
        assert(SO3dim==3)
        #for j in range(SO3dim):
        ### TODO: make generic amax[1] -> amax[SO3]
        constraintsstr += "\n" + ' '.join([str(amax[1]),str(0.0),str(0.0)]) 
        #print constraintsstr
        #print "###########################"

        for v in inertia:
            constraintsstr += "\n" + ' '.join([str(i) for i in v])

        abc = TOPPbindings.RunComputeSO3Constraints(SO3trajstr,constraintsstr)
        a,b,c = Extractabc(abc)
        return [a,b,c]

def getSpeedProfile( F, R, amin, amax, W, dW, ddW ):
        RdimIdx = [0,1,2]
        Rdim = 3
        SO3dim = 1
        SO3dimIdx = [3,None,None]

        Ndim = W.shape[0]
        Nwaypoints = W.shape[1]

        duration = 1.0/float(Nwaypoints-1)
        Rtrajstr = GetTrajectoryStringR(RdimIdx, W, dW, ddW)

        Rtraj = Trajectory.PiecewisePolynomialTrajectory.FromString(Rtrajstr)

        print Rtraj.Eval(1.0)
        print "###############"
        Rtraj.Plot(0.1)
        plt.show()
        #sys.exit(0)

        ### compute a,b,c
        qs = np.zeros((Rdim,Nwaypoints))
        qss = np.zeros((Rdim,Nwaypoints))
        for i in range(0,Nwaypoints):
                qs[:,i] = Rtraj.Evald(i*duration)
                qss[:,i] = Rtraj.Evaldd(i*duration)
                print i*duration,qs[:,i],qss[:,i]
        print "###############"

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
                for i in RdimIdx:
                        if H2[i] > -H1[i]:
                                print H2[i],"<= q[",i,"]<=",-H1[i]
                                sys.exit(1)

                cr[i,:] = np.hstack((H1,H2)).flatten()

        for i in range(0,Nwaypoints):
                ar[i,:] = np.dot(G,qs[:,i]).flatten()
                br[i,:] = np.dot(G,qss[:,i]).flatten()

        print ar[0,:],br[0,:],cr[0,:]
        print qs[:,1]
        print qss[:,0]
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

        topp_inst = TOPP.QuadraticConstraints(traj0, duration, vmax, list(a), list(b), list(c))
        x = topp_inst.solver

        x.integrationtimestep = 0.001
        x.reparamtimestep = 0.001
        x.extrareps = 10

        ret = x.RunComputeProfiles(0.0,0.0)
        x.ReparameterizeTrajectory()

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

                traj1 = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
                print "Trajectory duration before TOPP: ", traj0.duration
                print "Trajectory duration after TOPP: ", traj1.duration
                t = 0
                P = []
                tstep = traj1.duration/Nwaypoints
                while t < traj1.duration:
                        P.append(np.linalg.norm(traj1.Evald(t)))
                        t = t + tstep
                print "profile computed"
                #mvcbobrow = profileslist.pop(0)
                #mvc = np.array(mvcbobrow[3])
                #print mvc
                return np.array(P)
        return None



def piecewisePolynomialTrajectoryFromWaypoints(W):
        Ndim = W.shape[0]
        Nwaypoints = W.shape[1]
        Trajectory.PiecewisePolynomialTrajectory(chunkslist)

def GetMVC(W, dW, u1min, u1max, u2min, u2max, dF):
        Ndim = W.shape[0]
        Nwaypoints = W.shape[1]

        trajectorystring = str(1.0)
        trajectorystring += "\n" + str(Ndim)
        for i in range(Ndim):
                trajectorystring += "\n" + string.join([str(w) for w in W[i,:]])

        ## Trajectory
        #ndof = 5
        #trajectorystring = """1.0
        #5
        #-0.495010 1.748820 -2.857899 1.965396
        #0.008319 0.004494 1.357524 -1.289918
        #-0.354081 1.801074 -1.820616 0.560944
        #0.221734 -1.320792 3.297177 -2.669786
        #-0.137741 0.105246 0.118968 -0.051712
        #1.0
        #5
        #0.361307 1.929207 -4.349490 2.275776
        #0.080419 -1.150212 2.511645 -1.835906
        #0.187321 -0.157326 -0.355785 0.111770
        #-0.471667 -2.735793 7.490559 -4.501124
        #0.034761 0.188049 -1.298730 1.553443"""

        traj0 = Trajectory.PiecewisePolynomialTrajectory.FromString(trajectorystring)

        # Constraints
        vmax = 200*ones(Ndim)
        vmax[2] = 0.0 ##z-axis
        amax = 200*ones(Ndim)

        dt = 0.01
        t0 = time.time()

        #constraintstring = str(dt)
        #constraintstring += "\n" + string.join([str(v) for v in vmax])
        #constraintstring += "\n" + string.join([str(a) for a in amax])
        #x = TOPPbindings.TOPPInstance(None,"KinematicLimits",constraintstring,trajectorystring);
        #else: #Using the general QuadraticConstraints (fully supported)

        constraintstring = str(dt)
        constraintstring += "\n" + string.join([str(v) for v in vmax])
        constraintstring += TOPPpy.ComputeKinematicConstraints(traj0, amax, dt) 
        x = TOPPbindings.TOPPInstance(None,"QuadraticConstraints",constraintstring,trajectorystring);

        x.integrationtimestep = 0.005
        x.reparamtimestep = 0.02
        x.extrareps = 10
        #x.passswitchpointnsteps = 5

        # Run TOPP
        t1 = time.time()
        ret = x.RunComputeProfiles(0,0)
        x.ReparameterizeTrajectory()
        t2 = time.time()

        print "Setup TOPP:", t1-t0
        print "Run TOPP:", t2-t1
        print "Total:", t2-t0

        # Display results
        ion()
        #x.WriteProfilesList()
        #x.WriteSwitchPointsList()
        #profileslist = list(TOPPpy.ProfilesFromString(x.resprofilesliststring))
        #switchpointslist = TOPPpy.SwitchPointsFromString(x.switchpointsliststring)
        #TOPPpy.PlotProfiles(profileslist,switchpointslist,4)
        x.WriteResultTrajectory()

        traj1 = Trajectory.PiecewisePolynomialTrajectory.FromString(x.restrajectorystring)
        print "Trajectory duration before TOPP: ", traj0.duration
        print "Trajectory duration after TOPP: ", traj1.duration
        t = 0
        P = []
        tstep = traj1.duration/Nwaypoints
        while t < traj1.duration:
                P.append(np.linalg.norm(traj1.Evald(t)))
                t = t + tstep
        print "profile computed"
        #mvcbobrow = profileslist.pop(0)
        #mvc = np.array(mvcbobrow[3])
        #print mvc
        return np.array(P)

