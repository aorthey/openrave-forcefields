

from util import Rz

### return the disturbance force after we counteract them with the best
### acceleration available (from our dt reachable set)
def get_minimal_disturbance_forces(self, dt, W, F, u1min, u1max, u2min, u2max):
        dt2 = dt*dt/2
        Ndim = W.shape[0]
        Nsamples = W.shape[1]
        Fp = np.zeros((Ndim,Nsamples))

        for i in range(0,Nsamples):
                theta = W[3,i]
                Fp[0:3,i] = np.dot(Rz(theta),F[0:3,i])[0:3].flatten()

        ### question (1) : can we produce an acceleration ddW, such that it counteracts F?
        nearest_ddq = ForceCanBeCounteractedByAccelerationVector(dt, Fp, u1min, u1max, u2min, u2max)

        Ndim = Fp.shape[0]
        Nsamples = Fp.shape[1]
        dF = np.zeros((Ndim,Nsamples))

        for i in range(0,Nsamples):
                for j in range(0,Ndim):
                        dF[j,i] = (Fp[j,i] - nearest_ddq[j,i])

        for i in range(0,Nsamples):
                theta = W[3,i]
                dF[0:3,i] = np.dot(Rz(theta).T,dF[0:3,i])[0:3].flatten()
        return dF

def get_minimum_feasible_velocity_vector(self, dt, epsilon, W, dW, F):

        u1min = 0.0
        u1max = 2.0
        u2min = -2.0
        u2max = 2.0
        dt2 = dt*dt/2

        ## assume that W contains [x,y,z,theta]
        ## note: a force along x,y,z corresponds to a linear acc field, a force along theta would be a twisted acc field 

        ## project F onto the identity element, here the identity is at (0,0,0,0), with x pointing towards the 'right'
        Ndim = W.shape[0]
        Nsamples = W.shape[1]
        Fp = np.zeros((2,Nsamples))
        for i in range(0,Nsamples):
                theta = W[3,i]
                Fp[:,i] = np.dot(Rz(theta).T,F[:,i])[0:2].flatten()

        ### question (1) : can we produce an acceleration ddW, such that it counteracts F?
        nearest_ddq = ForceCanBeCounteractedByAccelerationVector(dt, Fp, u1min, u1max, u2min, u2max)

        ### question (2) : if ddW cannot counteract F, what is the minimal speed dW, such that we still follow tau?
        ## check that we are in an epsilon neighborhood of next waypoint!
        Ndim = Fp.shape[0]
        Nsamples = Fp.shape[1]
        dF = np.zeros((Ndim,Nsamples))

        for i in range(0,Nsamples):
                for j in range(0,Ndim):
                        dF[j,i] = dt2*Fp[j,i] - nearest_ddq[j,i]

        vmin = GetMinimalSpeedToReachEpsilonNeighbordhoodVector(dt, epsilon, W, dW, dF)

        return vmin

def get_minimum_feasible_velocity(self, dt, W0, W1, dW0, F):
        epsilon = 0.01

        ## assume that W contains [x,y,z,theta]
        ## note: a force along x,y,z corresponds to a linear acc field, a force along theta would be a twisted acc field 

        ## project F onto the identity element, here the identity is at (0,0,0,0), with x pointing towards the 'right'
        theta = W0[3]
        R = Rz(theta)
        Fp = np.dot(R.T,F)
        Fp = Fp[0:2].flatten()

        u1min = 0.0
        u1max = 2.0
        u2min = -2.0
        u2max = 2.0

        dt2 = dt*dt/2

        ### question (1) : can we produce an acceleration ddW, such that it counteracts F?
        [d,nearest_ddq] = ForceCanBeCounteractedByAcceleration(dt, Fp, u1min, u1max, u2min, u2max)

        if d <= 0.0001:
                #print "force can be counteracted by system dynamics"
                #print "no need for having a fixed velocity"
                return 0.0
        
        ### question (2) : if ddW cannot counteract F, what is the minimal speed dW, such that we still follow tau?
        ## check that we are in an epsilon neighborhood of next waypoint!
        dF = np.zeros((2,1))
        dF[0] = dt2*Fp[0] - nearest_ddq[0]
        dF[1] = dt2*Fp[1] - nearest_ddq[1]


        vmin = GetMinimalSpeedToReachEpsilonNeighbordhood(dt, epsilon, W0, W1, dW0, dF)

        return vmin


def normalize_velocity_vector(self, dW):
        Ndim = dW.shape[1]
        for i in range(0,Ndim):
                dW[:,i] = dW[:,i]/np.linalg.norm(dW[:,i])
        return dW

def plot_speed_profile(self, env):

        L = self.get_length()
        dt = 0.05
        Nwaypoints = int(L/dt)

        [W,dW] = self.get_waypoints(N=Nwaypoints)
        F = self.get_forces_at_waypoints(W, env)
        dW = self.normalize_velocity_vector(dW)
        ### we need F, W, dW(normalized)

        #V = np.zeros((Nwaypoints-1))
        #for i in range(0,Nwaypoints-1):
                #V[i] = self.get_minimum_feasible_velocity(dt, W[:,i],W[:,i+1],dW[:,i],F[:,i])

        epsilon = 0.01
        T = np.linspace(0.0,1.0,num=Nwaypoints-1)
        V = self.get_minimum_feasible_velocity_vector(dt,epsilon, W,dW,F)

        plot(T,V,'or',linewidth=3,markersize=2)        
        title('Speed Profile Trajectory', fontsize=self.FONT_SIZE)
        xlabel('$\theta$', fontsize=self.FONT_SIZE)
        ylabel('$s$', fontsize=self.FONT_SIZE)
        plt.show()

def plot_speed_profile(self, P):
        Nwaypoints = P.shape[0]
        S = np.linspace(0,1,Nwaypoints)

        fig=figure(facecolor='white')
        plot(S,P,'-r',linewidth=5,markersize=8)        
        plt.tick_params(labelsize=self.FONT_SIZE)
        title('Speed Profile Trajectory', fontsize=self.FONT_SIZE)
        xlabel('$s$', fontsize=2*self.FONT_SIZE)
        ylabel('$\dot s$', fontsize=2*self.FONT_SIZE)

        plt.show()
def computeDelta(self, W, dW, Qoff, Rmin, Rmax):
        Ndim = Rmin.shape[0]
        eta = cvx.Variable(1)
        xvar = cvx.Variable(Ndim)
        yvar = cvx.Variable(Ndim)

        constraints = []
        constraints.append( Qoff + Rmin <= xvar )
        constraints.append( xvar <= Qoff + Rmax ) ## Q
        constraints.append( W + eta*dW == yvar ) ## L
        constraints.append( eta > 0)
        objfunc = cvx.norm(xvar-yvar)

        objective = cvx.Minimize(objfunc)
        prob = cvx.Problem(objective, constraints)
        result = prob.solve(solver=cvx.SCS, eps=1e-5)
        #result = prob.solve(solver=SCS, use_indirect=True, eps=1e-2, verbose=True)
        #prob.solve(verbose=True, abstol_inacc=1e-2,reltol_inacc=1e-2,max_iters= 300, reltol=1e-2)

        delta = np.linalg.norm(xvar.value-yvar.value)

        Qnearest = np.array(xvar.value).flatten()
        Lnearest = np.array(yvar.value).flatten()

        dL = dW/np.linalg.norm(dW)
        dQ = (Qnearest-W)/np.linalg.norm(Qnearest-W)
        angle = acos(np.dot(dL,dQ))

        return [delta, angle, Lnearest, Qnearest]

def drawReachableSet(self, env, dt, speed, W, dW, amin, amax, Lnearest, Qnearest, R, F, angle, delta):
        dt2 = dt*dt*0.5
        from scipy.spatial import ConvexHull
        import itertools
        Ndim = W.shape[0]

        witv = dt
        plot([W[0],W[0]+witv*dW[0]],[W[1],W[1]+witv*dW[1]],'-k',linewidth=3)
        plot([W[0]+witv*dW[0],W[0]+(witv+witv/3.0)*dW[0]],[W[1]+witv*dW[1],W[1]+(witv+witv/3.0)*dW[1]],'--k',linewidth=3)
        plot(W[0],W[1],'or',markersize=10)
        plot([Lnearest[0],Qnearest[0]],[Lnearest[1],Qnearest[1]],'--og',markersize=10,linewidth=4)
        plot([Lnearest[0]],[Lnearest[1]],'ok',markersize=20)
        Mlq = Lnearest+0.5*(Qnearest-Lnearest)
        strdelta = ' $\\delta$='+str(np.around(delta,decimals=2))
        plt.text(Mlq[0], Mlq[1], strdelta, fontsize=22)
        dQ = Qnearest-W
        plot([W[0],W[0]+dQ[0]],[W[1],W[1]+dQ[1]],'--g',linewidth=3)
        Mlq = W+0.2*(Qnearest-W)
        strangle = ' $\\theta$='+str(np.around(angle,decimals=2))
        plt.text(Mlq[0], Mlq[1], strangle, fontsize=22)


        ctr = 0
        Adim = amin.shape[0]
        Mf = 2**Adim
        Qt = np.zeros((Mf,Ndim))
        for val in itertools.product([0,1], repeat=Adim):
                ## if 0, then choose amin, if 1, then choose amax
                a = val*amax + abs(1-np.array(val))*amin
                Qt[ctr,:] = W + speed*dt*dW + dt2*F.flatten() + dt2*np.dot(R,a) 

                ctr=ctr+1

        #plt.fill_betweenx(Qt[:,0], Qt[:,1])
        plot(Qt[:,0],Qt[:,1],'ob',markersize=5)
        hull = ConvexHull(Qt[:,0:2])
        plt.fill(Qt[hull.vertices,0], Qt[hull.vertices,1], facecolor='blue', alpha=0.2, lw=2)
        plt.show()

def visualizeReachableSetsAtT(self, env, t0, dt = 0.1):
        speed = 0.0
        dt2 = dt*dt*0.5
        [W,dW,ddW] = self.evaluate_at(t0,der=2)
        print "REACHABLE SET AT POSITION:",W
        [Ndim, Nwaypoints] = self.getWaypointDim(W)
        F = self.get_forces_at_waypoints(W, env)
        #F = np.zeros((F.shape))
        [R,amin,amax] = self.getControlMatrix(W)
        Rmin = np.zeros((Ndim))
        Rmax = np.zeros((Ndim))
        Rmin = dt2*np.minimum(np.dot(R[:,:,0],amax),np.dot(R[:,:,0],amin))
        Rmax = dt2*np.maximum(np.dot(R[:,:,0],amax),np.dot(R[:,:,0],amin))

        Qoff = W + speed*dt*dW + dt2*F.flatten()
        [delta0, angle0, Lnearest, Qnearest] = self.computeDelta(W, dW, Qoff, Rmin, Rmax)
        self.drawReachableSet(env, dt, speed, W, dW, amin, amax, Lnearest, Qnearest, R[:,:,0], F, angle0, delta0)
        
        M = 50
        lambda1 = np.linspace(-0.02,0.02,M)
        lambda2 = np.linspace(0.0,1,M)

        delta = np.zeros((M,M))
        angle = np.zeros((M,M))
        ctri = 0
        for l1 in lambda1:
                ctrj = 0
                for l2 in lambda2:
                        speed = sqrt(l2*amax[0]*2)
                        trajNormal = np.array((-W[1],W[0],0.0,0.0))
                        trajNormal = trajNormal/np.linalg.norm(trajNormal)

                        Fproj = np.dot(F.flatten().T, trajNormal)
                        Qoff = W + speed*dt*dW + dt2*F.flatten() - l1*Fproj*trajNormal
                        [delta[ctri,ctrj], angle[ctri,ctrj], Lnearest, Qnearest] = self.computeDelta(W, dW, Qoff, Rmin, Rmax)
                        ctrj = ctrj+1
                ctri = ctri+1

        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib import cm
        from matplotlib.ticker import LinearLocator, FormatStrFormatter
        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        X, Y = np.meshgrid(lambda2, lambda1)
        surf = ax.plot_surface(X, Y, angle, rstride=1, cstride=1, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)

        ax.set_xlabel('$\\lambda_2$',fontsize=22)
        ax.set_ylabel('$\\lambda_1$',fontsize=22)
        ax.set_zlabel('$\\theta$',fontsize=22)

        dmax = np.amax(delta)
        ax.plot([0],[0],[angle0],'ob',markersize=10)
        #ax.set_zlim(-1.01, 1.01)

        #fig.colorbar(surf, shrink=0.5, aspect=5)

        plt.show()
        #print delta

def computeReachableSets(self, dt, env):
        [W,dW,ddW] = self.get_waypoints_second_order()

        [Ndim, Nwaypoints] = self.getWaypointDim(W)

        F = self.get_forces_at_waypoints(W, env)
        [R,amin,amax] = self.getControlMatrix(W)

        Rmin = np.zeros((Nwaypoints, Ndim))
        Rmax = np.zeros((Nwaypoints, Ndim))
        for i in range(0,Nwaypoints):
                Rmin[i,:] = np.minimum(np.dot(R[:,:,i],amax),np.dot(R[:,:,i],amin))
                Rmax[i,:] = np.maximum(np.dot(R[:,:,i],amax),np.dot(R[:,:,i],amin))

        dt2 = dt*dt*0.5

        delta = np.zeros((Nwaypoints))
        eta = cvx.Variable(Nwaypoints)
        xvar = cvx.Variable(Ndim,Nwaypoints)
        yvar = cvx.Variable(Ndim,Nwaypoints)
        constraints = []
        objfunc = 0.0
        for i in range(0,Nwaypoints):
                Qoff = W[:,i] + dt*dW[:,i] + dt2*F[:,i]
                constraints.append( Qoff + Rmin[i,:] <= xvar[:,i] )
                constraints.append( xvar[:,i] <= Qoff + Rmax[i,:] ) ## Q
                constraints.append( W[:,i] + eta[i]*dW[:,i] == yvar[:,i] ) ## L
                constraints.append( eta[i] >= 0)
                objfunc += cvx.norm(xvar[:,i]-yvar[:,i])

        objective = cvx.Minimize(objfunc)
        prob = cvx.Problem(objective, constraints)

        print "solve minimal speed"
        result = prob.solve(solver=cvx.SCS)
        for i in range(0,Nwaypoints):
                delta[i] = np.linalg.norm(xvar[:,i].value-yvar[:,i].value)

        plt.clf()
        T = np.linspace(0,1.0,Nwaypoints)
        subplot(2,1,1)
        title('Delta and Force Profiles', fontsize=self.FONT_SIZE)
        plot(T,delta,'-r',linewidth=4)
        ylabel('Delta $\delta$', fontsize=self.FONT_SIZE)
        subplot(2,1,2)
        plot(T,F[0,:],'-b',linewidth=4,label="$F_x$")
        plot(T,F[1,:],'-g',linewidth=4,label="$F_y$")
        xlabel('Path parameter $s$',fontsize=self.FONT_SIZE)
        ylabel('Force $F$', fontsize=self.FONT_SIZE)
        legend = plt.legend(loc='upper center', shadow=True, fontsize=self.FONT_SIZE)
        plt.show()
