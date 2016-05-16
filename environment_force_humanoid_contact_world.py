from environment_force_humanoid import *

class EnvironmentHumanoidContactWorld(EnvironmentHumanoid):
        def __init__(self):
                xmlenv='environments/contactworld.env.xml'
                EnvironmentHumanoid.__init__(self, xmlenv)

        def GetCells(self):
                C = self.GetCellsAll()
                self.cells = C[0:6]
                return self.cells

        def GetForces(self):
                ##
                self.forces = np.array((0.0,0.0,0.0))
                self.forces = np.vstack([self.forces,(0.0,2.0,0.0)])
                self.forces = np.vstack([self.forces,(0.0,0.0,0.0)])
                self.forces = np.vstack([self.forces,(0.0,0.0,0.0)])
                self.forces = np.vstack([self.forces,(0.0,0.0,0.0)])
                self.forces = np.vstack([self.forces,(0.0,0.0,0.0)])
                return self.forces

        def GetSurfaces(self):
                if self.cells is None:
                        self.cells = self.GetCells()

                S = []
                for i in range(0,len(self.cells)):
                        print i
                        C = self.cells[i]
                        G1 = C.GetGeometries()[0]
                        B = G1.GetBoxExtents()
                        T = G1.GetTransform()
                        #self.DrawFrameFromTransform(T)

                        bx = B[0]
                        by = B[1]
                        bz = B[2]

                        ## position of center of plane
                        pN1 = np.dot(T,np.array((bx,0,0,1)))[0:3]
                        pN2 = np.dot(T,np.array((-bx,0,0,1)))[0:3]
                        pN3 = np.dot(T,np.array((0,by,0,1)))[0:3]
                        pN4 = np.dot(T,np.array((0,-by,0,1)))[0:3]
                        pN5 = np.dot(T,np.array((0,0,bz,1)))[0:3]
                        pN6 = np.dot(T,np.array((0,0,-bz,1)))[0:3]

                        s = 0.1
                        nN1 = np.dot(T,np.array((s+bx,0,0,1)))[0:3]
                        nN2 = np.dot(T,np.array((-(s+bx),0,0,1)))[0:3]
                        nN3 = np.dot(T,np.array((0,s+by,0,1)))[0:3]
                        nN4 = np.dot(T,np.array((0,-(s+by),0,1)))[0:3]
                        nN5 = np.dot(T,np.array((0,0,s+bz,1)))[0:3]
                        nN6 = np.dot(T,np.array((0,0,-(s+bz),1)))[0:3]

                        ## normed direction of normal of plane
                        center_box = np.dot(T,np.array((0,0,0,1)))[0:3]
                        
                        dN1 = np.dot(T,np.array((1,0,0,1)))[0:3]-center_box
                        dN2 = np.dot(T,np.array((-(1),0,0,1)))[0:3]-center_box
                        dN3 = np.dot(T,np.array((0,1,0,1)))[0:3]-center_box
                        dN4 = np.dot(T,np.array((0,-(1),0,1)))[0:3]-center_box
                        dN5 = np.dot(T,np.array((0,0,1,1)))[0:3]-center_box
                        dN6 = np.dot(T,np.array((0,0,-(1),1)))[0:3]-center_box

                        dN1 = dN1/np.linalg.norm(dN1)
                        dN2 = dN2/np.linalg.norm(dN2)
                        dN3 = dN3/np.linalg.norm(dN3)
                        dN4 = dN4/np.linalg.norm(dN4)
                        dN5 = dN5/np.linalg.norm(dN5)
                        dN6 = dN6/np.linalg.norm(dN6)

                        posN = np.array((pN1,pN2,pN3,pN4,pN5,pN6))
                        dirN = np.array((dN1,dN2,dN3,dN4,dN5,dN6))


                        ## tangential,and orthogonal tangential vector of plane
                        ## make sure they have right-handed orientation
                        tN1 = dN3
                        tN2 = dN3
                        tN3 = dN5
                        tN4 = dN5
                        tN5 = dN1
                        tN6 = dN1

                        oN1 = np.cross(dN1,tN1)
                        oN2 = np.cross(dN2,tN2)
                        oN3 = np.cross(dN3,tN3)
                        oN4 = np.cross(dN4,tN4)
                        oN5 = np.cross(dN5,tN5)
                        oN6 = np.cross(dN6,tN6)

                        tN = np.array((tN1,tN2,tN3,tN4,tN5,tN6))
                        oN = np.array((oN1,oN2,oN3,oN4,oN5,oN6))

                        #self.DrawNormalVectors(posN,dirN)
                        #self.DrawTangentialVectors(posN,tN,oN)

                        self.DrawCoordinateFrame(posN, dirN, tN, oN)

                        ## extension into tN1 and oN1
                        etN1 = by
                        etN2 = by
                        etN3 = bz
                        etN4 = bz
                        etN5 = bx
                        etN6 = bx

                        eoN1 = bz
                        eoN2 = bz
                        eoN3 = bx
                        eoN4 = bx
                        eoN5 = by
                        eoN6 = by

                        ### pos of center, dir of normal, dir of tangential, dir
                        ### of binormal, extension in tangential dir, extension in
                        ### binormal dir
                        Si = np.zeros((6,6,3))
                        Si[0,:,:] = np.array((pN1,dN1,tN1,oN1,np.array((etN1,0,0)),np.array((eoN1,0,0))))
                        Si[1,:,:] = np.array((pN2,dN2,tN2,oN2,np.array((etN2,0,0)),np.array((eoN2,0,0))))
                        Si[2,:,:] = np.array((pN3,dN3,tN3,oN3,np.array((etN3,0,0)),np.array((eoN3,0,0))))
                        Si[3,:,:] = np.array((pN4,dN4,tN4,oN4,np.array((etN4,0,0)),np.array((eoN4,0,0))))
                        Si[4,:,:] = np.array((pN5,dN5,tN5,oN5,np.array((etN5,0,0)),np.array((eoN5,0,0))))
                        Si[5,:,:] = np.array((pN6,dN6,tN6,oN6,np.array((etN6,0,0)),np.array((eoN6,0,0))))

                        if i==0:
                                S = Si
                        else:
                                S = np.vstack((S,Si))


                return S



        def DisplayForces(self):
                if self.forces is None:
                        self.forces = self.GetForces()
                if self.cells is None:
                        self.cells = self.GetCells()

                #with self.env:
                assert(len(self.forces)==len(self.cells))

                self.ResetForceHandles()
                for i in range(0,len(self.cells)):
                        C = self.cells[i]
                        G1 = C.GetGeometries()[0]
                        B = G1.GetBoxExtents()
                        T = G1.GetTransform()
                        T[2,3] = 1.0
                        B[2] = 1.0
                        ######
                        h = self.DrawBoxMesh(T,B)
                        self.AddForceHandles([h])
                        ######
                        F = self.forces[i]
                        h = self.DrawForceArrowsInBox(T, B, F)
                        self.AddForceHandles(h)

if __name__ == "__main__":
        env = EnvironmentHumanoidContactWorld()
