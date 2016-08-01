from TOPP import Utilities
import TOPP
import numpy as np

def TOPPInterface(Ndim, trajectorystring, discrtimestep, a, b, c):
        vmax = 1e5*np.ones(Ndim) ## no velocity constraints
        topp_inst = TOPP.QuadraticConstraints(trajectorystring, discrtimestep, vmax, list(a), list(b), list(c))
        x = topp_inst.solver

        ##### compute all possible speeds at end of traj
        ret = x.RunVIP(0,0)
        print "discrtimestep",discrtimestep,"TOPP RunVIP             (0,0) code:",ret,
        if ret == 1:
                print "(Success)",
                semin = x.sdendmin
                semax = x.sdendmax
                print "sd_end: [",semin,",",semax,"]"

                ##### compute speed profile for some speeds inside sd_end
                for speed in np.linspace(semin,semax,5):
                        ret = x.RunComputeProfiles(0,speed)
                        print "discrtimestep",discrtimestep,"TOPP RunComputeProfiles (0,",speed,") code:",ret,

                        if ret==1:
                                print "(Success)"
                        else:
                                print "(Failure)"


Ndim = 4
discrtimestep= 1e-3
a = np.loadtxt('topp/a')
b = np.loadtxt('topp/b')
c = np.loadtxt('topp/c')
with open("topp/traj0", "r") as fh:
        trajectorystring = "%s" % fh.read()
#print trajectorystring
TOPPInterface(Ndim, trajectorystring, discrtimestep, a, b, c)
