class EnvironmentForces():
        def GetCells(self,env):
                B = env.GetBodies()[1]
                W = B.GetLinks()
                for b in W:
                        print b.GetName()
                return W

        def GetForceArrowsInsideBox(self,W1,Fx,Fy):
            G1 = W1.GetGeometries()[0]
            B = G1.GetBoxExtents()
            T = G1.GetTransform()
            print B,T
            bx = T[0,3]
            by = T[1,3]
            N = int(B[0])
            M = int(B[1])
            dxspacing = B[0]/N
            dyspacing = B[1]/M

            dx = min(dxspacing-dxspacing/10,Fx)
            dy = min(dyspacing-dyspacing/10,Fy)
            if (dx == 0) & (dy==0):
                    dx=dx+0.001

            l = sqrt(dx*dx+dy*dy)
            print dx,dy,N,M

            print dxspacing,dyspacing
            handles=[]
            xstart = -B[0]+bx
            ystart = -B[1]+by
            zval=0.1
            for i in range(0,2*N):
                for j in range(0,2*M):
                    x = xstart+i*dxspacing+0.3*dxspacing
                    y = ystart+j*dyspacing+0.3*dyspacing
                    A = env.drawarrow(array((x,y,zval)),array((x+dx,y+dy,zval)),linewidth=0.08*l,color=array((1,0,0,0.5)))
                    handles.append(A)
            return handles


