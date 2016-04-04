

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
