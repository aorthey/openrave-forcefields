from environment_force_humanoid import *
from environment_force_humanoid_contact_world import *


class EnvironmentHumanoidContactWorldSurfboard(EnvironmentHumanoidContactWorld):
        def __init__(self, xmlenv='environments/surfboard.env.xml'):
                EnvironmentHumanoid.__init__(self, xmlenv)

        def GetCells(self):
                C = self.GetCellsAll(verbose=False)
                self.cells = C[0:3]
                return self.cells

        def GetForces(self):
                ##
                self.forces = np.array((0.0,0.0,0.0))
                self.forces = np.vstack([self.forces,(-0.5,0.0,0.0)])
                self.forces = np.vstack([self.forces,(0.0,0.0,0.0)])
                return self.forces


if __name__ == "__main__":
        env = EnvironmentHumanoidContactWorldSurfboard()
