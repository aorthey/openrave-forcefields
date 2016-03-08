from deformation_factory import *
from util_mvc import *

import sys

DEBUG = 0

#def compute_lambda_1(i, W, dW, dF, Vmax, Amax_linear, dt, lambda_2):
#        dt2 = dt[i]*dt[i]/2
#        A = -dF[:,i]
#        B = W[:,i]-W[:,i+1]+dt2*dF[:,i]+dW[:,i]*dt[i]*(Vmax[i]+sqrt(2*lambda_2[i]*Amax_linear[i]))
#        A2 = np.dot(A.T,A)
#        if A2 < 0.001:
#                lambda_1 = 0.0
#        else:
#                lambda_1 = - np.dot(A.T,B)/A2
#        return lambda_1

def eval_norm(i, W, dW, dF, Vmax, Amax_linear, dt, lambda_1, lambda_2):
        dt2 = dt[i]*dt[i]/2
        #print dt2,dF[:,i],dW[:,i],dt[i],Vmax[i],Amax_linear[i],lambda_2[i]
        #print Vmax[i]

        
        ### project dF onto dW to obtain dWf. then subtract dWf from dF to
        ### obtain the component orthogonal to dW
        dWf = np.dot(dF[:,i].T,dW[:,i])*dW[:,i]
        dFcomponent = dF[:,i] - dWf

        Wst = dt2*dF[:,i]+dW[:,i]*dt[i]*(Vmax[i]+sqrt(2*lambda_2[i]*Amax_linear[i]))+W[:,i]-lambda_1[i]*dFcomponent
        d = np.linalg.norm(Wst-W[:,i+1])
        return [d,Wst]

def compute_minimum_distance(i, W, dW, dF, Vmax, Amax_linear, dt, lambda_1, lambda_2):
        dold = 2000.0
        d = 1000.0
        dt[i]=0.0
        tstep = 0.005
        Wstarr = W[:,i]
        while d < dold:
                dold = d
                dt[i] = dt[i]+tstep
                [d,Wst] = eval_norm(i, W, dW, dF, Vmax, Amax_linear, dt, lambda_1, lambda_2)
        return dold

def eval_norm2(i, W, dW, dF, Vmax, Amax_linear, dt, lambda_2):
        dt2 = dt[i]*dt[i]/2
        Wst = dt2*dF[:,i]+dW[:,i]*dt[i]*(Vmax[i-1]+lambda_2[i]*sqrt(2*Amax_linear[i]))+W[:,i]
        d = np.linalg.norm(Wst-W[:,i+1])
        return [d,Wst]

def compute_minimum_distance2(i, W, dW, dF, Vmax, Amax_linear, dt, lambda_2):
        dold = 2000.0
        d = 1000.0
        dt[i]=0.0
        tstep = 0.005
        Wstarr = W[:,i]
        while d < dold:
                dold = d
                dt[i] = dt[i]+tstep
                [d,Wst] = eval_norm2(i, W, dW, dF, Vmax, Amax_linear, dt, lambda_2)
        return dold

def compute_lambda_updates(epsilon, W, dW, dF, Vmax, Amax_linear):
        #[lambda_1, lambda_2] = compute_lambda_updates(W, dW, dF, Vmax, Amax_linear)
        DEBUG = 1
        if W.ndim <= 1:
                print "Cannot compute with only one waypoint. Abort"
                sys.exit(1)

        Ndim = W.shape[0]
        Nwaypoints = W.shape[1]

        lambda_1 = np.zeros(Nwaypoints)
        lambda_2 = np.zeros(Nwaypoints)
        lambda_zero = np.zeros(Nwaypoints)
        dt = np.zeros((Nwaypoints))

        lstep = 0.01
        for i in range(0,Nwaypoints-1):

                ddf = np.linalg.norm(dF[:,i])
                if ddf > 0:
                        dold = 2000.0
                        d_lambda = compute_minimum_distance(i, W, dW, dF, Vmax, Amax_linear, dt, lambda_1, lambda_zero)
                        while d_lambda < dold:
                                dold = d_lambda
                                lambda_1[i] += lstep 
                                d_lambda = compute_minimum_distance(i, W, dW, dF, Vmax, Amax_linear, dt, lambda_1, lambda_zero)

        for i in range(1,Nwaypoints-1):
                print i,"/",Nwaypoints
                d_lambda = compute_minimum_distance2(i, W, dW, dF, Vmax, Amax_linear, dt, lambda_2)
                if d_lambda < epsilon:
                        ## step down
                        dold = d_lambda+1.0
                        while d_lambda < dold:
                                lambda_2[i] -= lstep
                                dold = d_lambda
                                d_lambda = compute_minimum_distance2(i, W, dW, dF, Vmax, Amax_linear, dt, lambda_2)
                        lambda_2[i] += lstep
                else:
                        while d_lambda > epsilon:
                                lambda_2[i] += lstep 
                                d_lambda = compute_minimum_distance2(i, W, dW, dF, Vmax, Amax_linear, dt, lambda_2)

        if DEBUG:
                plt.subplot(2,1,1)
                plot(lambda_1,'-r')
                plt.ylabel('$\lambda_1$',fontsize=25)
                plt.subplot(2,1,2)
                plot(lambda_2,'-r')
                plt.ylabel('$\lambda_2$',fontsize=25)
                plt.xlabel('time (s)',fontsize=25)
                plt.show()


        return [lambda_1,lambda_2]


if __name__ == '__main__':
        plt.clf()
        dF = np.array(((0.3,-0.02),(0,0))).T
        W = np.array(((0,0),(2,0))).T
        dW = np.array(((1,0),(1,0))).T

        Amax_linear = np.array(((0.5),(0.5)))
        Vmax = np.array(((0.0),(0.0)))
        [l1,l2]=compute_lambda_updates(W, dW, dF, Vmax, Amax_linear)

        plot(W[0,:],W[1,:],'ok',markersize=10)
        plot([W[0,:], W[0,:]+dW[0,:]],[W[1,:], W[1,:]+dW[1,:]],'-r',linewidth=3)
        plot([W[0,:], W[0,:]+dF[0,:]],[W[1,:], W[1,:]+dF[1,:]],'-g',linewidth=3)
        plt.show()
        print l1,l2


class DeformationPotentials(Deformation):


        ## change only traj_deformed here
        def deform_onestep(self):
                u1min = 0.0
                u1max = 20.0
                u2min = -3.0
                u2max = 3.0
                dt = 0.01
                epsilon = 0.005

                traj = self.traj_deformed
                L = traj.get_length()
                Nwaypoints = int(L/dt)
                print "WAYPOINTS:",Nwaypoints,"LENGTH:",L

                [Wori,dWori] = traj.get_waypoints(N=Nwaypoints)#(self.traj_ori, N = 100)
                Ndim = Wori.shape[0]
                Nwaypoints = Wori.shape[1]

                ###############################################################
                ### compute dF from F (minimum disturbance force at each
                ### waypoint,  projected into the world frame)
                ###############################################################

                F = self.GetForcesAtWaypoints(Wori)
                dF = traj.get_minimal_disturbance_forces(dt, Wori, F, u1min, u1max, u2min, u2max)
                #self.draw_forces_at_waypoints(Wori, dF)

                ###############################################################
                ### compute the maximum velocity curve (MVC)
                ###############################################################

                #Vmax = GetMinimalSpeedToReachEpsilonNeighbordhoodVector(dt, epsilon, Wori, dWori, dF)
                #Vmax = GetMVC(Wori, dWori, u1min, u1max, u2min, u2max, dF)

                #Vmax = GetEpsilonVC(epsilon, Wori, dWori, u1min, u1max, u2min, u2max, dF)
                Vmax = get_minimal_epsilon_velocity_curve(epsilon, Wori, dWori, dF)
                traj.plot_speed_profile(Vmax)
                #sys.exit(0)

                ###############################################################
                ### compute orthogonal component of the disturbance force onto
                ### the direction vector
                ###############################################################

                dFcomponent = np.zeros(dF.shape)
                for i in range(0,Nwaypoints-1):
                        dWf = np.dot(dF[:,i].T,dWori[:,i])*dWori[:,i]
                        dFcomponent[:,i] = dF[:,i] - dWf

                ###############################################################
                ## for all points where dF = 0, we can say that lambda_1 = lambda_2 = 0
                ###############################################################

                Amax_linear = u1max*np.ones(Nwaypoints)
                [lambda_1, lambda_2] = compute_lambda_updates(epsilon, Wori, dWori, dF, Vmax, Amax_linear)

                ###############################################################
                ## update trajectory into lambda directions
                ###############################################################

                dU = np.zeros((Ndim,Nwaypoints))
                for i in range(0,Nwaypoints):
                        dU[:,i] = -lambda_1[i] * dFcomponent[:,i] - 3*lambda_2[i] * dWori[:,i]

                eta = 10
                W = Wori + eta*dU
                self.traj_deformed.new_from_waypoints(W)

