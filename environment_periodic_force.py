import abc
import time
import openravepy
if not __openravepy_build_doc__:
    from openravepy import *
    from numpy import *
from environment_force import *

class PeriodicForceEnvironment(ForceEnvironment):
        forces_max = None
        forces_min = None
        forces_period = None
        time_init = None

        force_handle = None

        def __init__(self):
                ForceEnvironment.__init__(self)
                self.time_init = time.time()

        @abc.abstractmethod
        def GetCells(self):
                pass

        @abc.abstractmethod
        def GetForces(self):
                pass

        @abc.abstractmethod
        def RobotGetInitialPosition(self):
                pass

        @abc.abstractmethod
        def RobotGetGoalPosition(self):
                pass

        def ComputeForceFromTime(self, t, Fmax, Fmin, Fperiod):
                Ft = numpy.array((0.0,0.0,0.0))
                for i in range(0,3):
                        if Fperiod > 0.001:
                                Ft[i]=0.5*(Fmax[i]-Fmin[i])*cos(-2*pi*t/Fperiod)+Fmin[i]
                        else:
                                Ft[i]=0.5*(Fmax[i]-Fmin[i])+Fmin[i]
                        #print Fmax[i],Fmin[i],t,Ft[i],Fperiod
                return Ft

        def DisplayForces(self):
                if self.forces_max is None:
                        self.forces = self.GetForces()
                if self.cells is None:
                        self.cells = self.GetCells()

                time_cur = time.time()
                t = time_cur - self.time_init

                assert(len(self.forces_max)==len(self.GetCells()))

                tmp_force_handle = []

                for i in range(0,len(self.cells)):
                        Ft = self.ComputeForceFromTime(t, self.forces_max[i], self.forces_min[i], self.forces_period[i])
                        h = self.DrawForceArrowsInCell(self.cells[i], Ft)
                        tmp_force_handle.append(h)
                        h = self.DrawBorderAroundCell(self.cells[i])
                        tmp_force_handle.append(h)
                self.force_handle = tmp_force_handle


if __name__ == "__main__":
        env = PeriodicForceEnvironment()
        xmlenv='environments/the_hide_out.env.xml'
        xmlrobot='robots/pointrobot.robot.xml'
        env.setrobotenv(xmlrobot,xmlenv)

