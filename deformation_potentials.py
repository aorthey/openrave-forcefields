from deformation_factory import *

import sys

DEBUG = 0

def compute_lambda_1(i, W, dW, dF, Vmax, Amax_linear, dt, lambda_2):
        dt2 = dt[i]*dt[i]/2
        A = -dF[:,i]
        B = W[:,i]-W[:,i+1]+dt2*dF[:,i]+dW[:,i]*dt[i]*(Vmax[i]+sqrt(2*lambda_2[i]*Amax_linear[i]))
        A2 = np.dot(A.T,A)
        if A2 < 0.001:
                lambda_1 = 0.0
        else:
                lambda_1 = - np.dot(A.T,B)/A2
        return lambda_1

def eval_norm(i, W, dW, dF, Vmax, Amax_linear, dt, lambda_1, lambda_2):
        dt2 = dt[i]*dt[i]/2
        #print dt2,dF[:,i],dW[:,i],dt[i],Vmax[i],Amax_linear[i],lambda_2[i]
        #print Vmax[i]
        Wst = dt2*dF[:,i]+dW[:,i]*dt[i]*(Vmax[i]+sqrt(2*lambda_2[i]*Amax_linear[i]))+W[:,i]-lambda_1[i]*dF[:,i]
        d = np.linalg.norm(Wst-W[:,i+1])
        return [d,Wst]

def compute_lambda_updates(W, dW, dF, Vmax, Amax_linear):
        #[lambda_1, lambda_2] = compute_lambda_updates(W, dW, dF, Vmax, Amax_linear)
        Ndim = W.shape[0]
        lstep = 0.01
        tstep = 0.001
        if W.ndim <= 1:
                print "Cannot compute with only one waypoint. Abort"
                sys.exit(1)

        Nwaypoints = W.shape[1]


        lambda_1 = np.zeros(Nwaypoints)
        lambda_2 = np.zeros(Nwaypoints)
        dt = np.zeros((Nwaypoints))
        for i in range(0,Nwaypoints-1):
                dold_lambda = 2000.0
                d_lambda = 1000.0

                while d_lambda < dold_lambda:
                        dold_lambda = d_lambda
                        dold = 2000.0
                        d = 1000.0
                        dt[i]=0.0
                        Wstarr = W[:,i]
                        while d < dold:
                                dold = d
                                dt[i] = dt[i]+tstep
                                #sys.stdout.write(str(np.around(dt[i],decimals=2))+"("+str(np.around(d,decimals=2))+')|')
                                #print np.around(dt[i],decimals=2),np.around(d,decimals=2)
                                [d,Wst] = eval_norm(i, W, dW, dF, Vmax, Amax_linear, dt, lambda_1, lambda_2)
                                Wstarr = np.vstack((Wstarr,Wst))

                        lambda_1[i]=compute_lambda_1(i, W, dW, dF, Vmax, Amax_linear, dt, lambda_2)
                        [d_lambda,Wst] = eval_norm(i, W, dW, dF, Vmax, Amax_linear, dt, lambda_1, lambda_2)
                        if DEBUG:
                                if d_lambda > dold_lambda:
                                        plot(Wstarr[:,0],Wstarr[:,1],'-',color=[1,0,1],linewidth=3)
                                else: 
                                        plot(Wstarr[:,0],Wstarr[:,1],'-',color=[0,0,1],linewidth=3)
                        lambda_2[i] += lstep 

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
                u1max = 2.0
                u2min = -2.0
                u2max = 2.0
                dt = 0.02
                epsilon = 0.001

                traj = self.traj_deformed
                L = traj.get_length()
                Nwaypoints = int(L/dt)

                [Wori,dWori] = traj.get_waypoints(N=100)#(self.traj_ori, N = 100)
                Ndim = Wori.shape[0]
                Nwaypoints = Wori.shape[1]

                ###############################################################
                ### compute dF from F (minimum disturbance force at each
                ### waypoint,  projected into the world frame)
                ###############################################################

                F = self.GetForcesAtWaypoints(Wori)
                dF = traj.get_minimal_disturbance_forces(dt, Wori, F, u1min, u1max, u2min, u2max)

                ###############################################################
                ### compute the maximum velocity curve (MVC)
                ###############################################################

                ##TODO!! integrate mvc from TOPP
                Vmax = GetMinimalSpeedToReachEpsilonNeighbordhoodVector(dt, epsilon, Wori, dWori, dF)
                self.Vmax=Vmax
                print Vmax

                ###############################################################
                ## for all points where dF = 0, we can say that lambda_1 = lambda_2 = 0
                ###############################################################

                Amax_linear = u1max*np.ones(Nwaypoints)
                [lambda_1, lambda_2] = compute_lambda_updates(Wori, dWori, dF, Vmax, Amax_linear)

                ###############################################################
                ## update trajectory into lambda directions
                ###############################################################

                dU = np.zeros((Ndim,Nwaypoints))
                print Wori.shape
                print dU.shape
                print lambda_1
                print lambda_2
                for i in range(0,Nwaypoints):
                        dU[:,i] = -lambda_1[i] * dF[:,i] - lambda_2[i] * dWori[:,i]

                eta = 100
                W = Wori + eta*dU
                self.traj_deformed.new_from_waypoints(W)


