from environment_force import *
import numpy as np
import openravepy

class EnvironmentRadial(ForceEnvironment):
        def __init__(self):
                ForceEnvironment.__init__(self)
                xmlenv='environments/the_radial.env.xml'
                xmlrobot='robots/virus.robot.xml'
                self.setrobotenv(xmlrobot,xmlenv)

        def GetCells(self):
                C = self.GetCellsAll()
                self.cells = C[0:3]
                return self.cells

        def GetForces(self):
                self.forces = np.array((2.5,2.5,0.0))
                self.forces = np.vstack([self.forces,(0.0,0.0,0.0)])
                self.forces = np.vstack([self.forces,(0.0,0.0,0.0)])
                return self.forces

        def RobotGetInitialPosition(self):
                return [-4.0,0.0,0.1,0,0,0,0,0]

        def RobotGetGoalPosition(self):
                return [4.0,0.0,0.1,0,0,0,0,0]

        def GetForceAtX(self, X):
                self.robot.SetDOFValues(X)

                F = self.GetZeroForce()
                for i in range(0,len(self.cells)):
                        robotIsInsideCell = self.env.CheckCollision(self.robot, self.cells[i])
                        if robotIsInsideCell:
                                #F = F + self.forces[i]
                                F += self.GetForceAtXInVectorField(X, i)

                return F

        def GetForceAtXInVectorField(self, X, i):
                F = np.zeros((3))
                if np.linalg.norm(self.forces[i])>0.001:
                        x=X[0]
                        y=X[1]
                        c = np.linalg.norm(X[0:2])**2
                        F[0] = -self.forces[i][0]*x/c
                        F[1] = -self.forces[i][1]*y/c

                return F

        def DisplayForces(self):
                if self.forces is None:
                        self.forces = self.GetForces()
                if self.cells is None:
                        self.cells = self.GetCells()

                #with self.env:
                assert(len(self.forces)==len(self.cells))

                self.ResetForceHandles()
                i = 0
                h = self.DrawForceArrowsInCell(self.cells[i], self.forces[i])
                self.AddForceHandles(h)

        def DrawForceArrowsInBox(self, T, B, F):
                ## force in box
                Fx = F[0]
                Fy = F[1]
                Fz = F[2]

                ## middle point of box
                mean_x = T[0,3]
                mean_y = T[1,3]
                mean_z = T[2,3]
                ## extend of box
                ext_x = B[0]
                ext_y = B[1]
                ext_z = B[2]

                Msamples = 15
                Ndim = 4
                handles=[]
                for x in np.linspace(-ext_x,ext_x,Msamples):
                        for y in np.linspace(-ext_y,ext_y,Msamples):
                                for z in np.linspace(ext_z,ext_z,Msamples):
                                        X=np.zeros((Ndim))
                                        X[0] = x
                                        X[1] = y
                                        X[2] = z

                                        if np.linalg.norm(X)>0.1:
                                                #X[0] = x+mean_x
                                                #X[1] = y+mean_y
                                                #X[2] = z+mean_z
                                                #F = self.GetForceAtX(X)
                                                F = self.GetForceAtXInVectorField(X, 0)

                                                X[0] = x#+mean_x
                                                X[1] = y-0.06#+mean_y
                                                X[2] = z/2#+mean_z

                                                #nF = 2*np.linalg.norm(F)
                                                nF=10
                                                Fx=F[0]/nF
                                                Fy=F[1]/nF
                                                Fz=F[2]/nF

                                                lx = Fx
                                                ly = Fy
                                                lz = Fz

                                                A = self.env.drawarrow(p1=array((X[0],X[1],X[2])),p2=array((X[0]+lx,X[1]+ly,X[2]+lz)),linewidth=self.FORCE_FIELD_ARROW_SIZE,color=self.FORCE_FIELD_COLOR)
                                                handles.append(A)

                return handles

if __name__ == "__main__":
        env = EnvironmentRadial()
        env.MakeRobotInvisible()
        env.DisplayForces()
