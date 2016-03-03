import string,time
from pylab import *
from numpy import *
from openravepy import *
from TOPP import Utilities
from TOPP import TOPPbindings
from TOPP import TOPPpy
from TOPP import Trajectory
from TOPP import TOPPopenravepy

#def FromRaveTraj(robot, traj):
#    N = traj.GetNumWaypoints()
#    if N < 2:
#        print "RAVE trajectory has less than 2 waypoints"
#        return None
#    timespec = openravepy.ConfigurationSpecification()
#    timespec.AddGroup('deltatime', 1, 'linear')
#    posspec = robot.GetActiveConfigurationSpecification()
#    velspec = posspec.ConvertToVelocitySpecification()
#    chunkslist = []
#    vel = traj.GetWaypoint(0, velspec)
#    ndof = len(vel)
#    for i in range(N - 1):
#        pos = traj.GetWaypoint(i, posspec)
#        deltatime = traj.GetWaypoint(i + 1, timespec)[0]
#        if deltatime < 1e-5:
#            continue
#        nextvel = traj.GetWaypoint(i + 1, velspec)
#        polynomialslist = []
#        for j in range(ndof):
#            x = pos[j]
#            v = vel[j]
#            a = (nextvel[j] - vel[j]) / deltatime
#            polynomialslist.append(Trajectory.Polynomial([x, v, a / 2]))
#        chunkslist.append(Trajectory.Chunk(deltatime, polynomialslist))
#        vel = nextvel
#    return Trajectory.PiecewisePolynomialTrajectory(chunkslist)


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

        print trajectorystring

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
        vmax = 2*ones(Ndim)
        vmax[2] = 0.0
        vmax[3] = 0.0
        amax = 0.1*ones(Ndim)

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
        tstep = traj1.duration/100
        while t < traj1.duration:
                P.append(np.linalg.norm(traj1.Evald(t)))
                t = t + tstep
        print "profile computed"
        #mvcbobrow = profileslist.pop(0)
        #mvc = np.array(mvcbobrow[3])
        #print mvc
        return np.array(P)

