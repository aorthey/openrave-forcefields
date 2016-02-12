import cvxpy as cvx
from cvxpy import Variable, Problem, Minimize, norm, SCS, OPTIMAL
from cvxpy import SCS, CVXOPT
from numpy import sqrt,sin,cos,pi
import numpy as np
from trajectory import *
import sys


## create a set of K orthonormal polynomial basis functions
def Fpoly(x,K): 
        ## assume x in [0,1]
        t = 2*x - 1.0
        #t = x
        return np.array(map(lambda e: map(lambda t: pow(t,e), t),np.arange(0,K)))

def dFpoly(x,K): 
        t = 2*x - 1.0
        A= np.array(map(lambda e: map(lambda t: pow(t,e-1)/e, t),np.arange(1,K)))
        O = np.zeros((1,len(t)))
        return np.vstack((O,A))

def Ffourier(x, K): 
        N = K/2
        ##sqrt(2)*sin(2*pi*N*x)|N\in |N} \cup sqrt(2)*cos(2*pi*N*x) \cup 1
        Fsin = map(lambda N: map(lambda x: sqrt(2)*sin(2*pi*N*x),x),np.arange(0,N))
        Fcos = map(lambda N: map(lambda x: sqrt(2)*cos(2*pi*N*x),x),np.arange(0,N))
        Fone = pow(x,0)
        F = np.vstack((Fsin,Fcos,Fone))
        #print N,len(x),F.shape
        return F


def find_nearest(array, value):
        idx = (np.abs(array-value)).argmin()
        return array[idx]

def find_nearest_idx(array, value):
        idx = (np.abs(array-value)).argmin()
        return idx

def WaypointsToWeights(waypts):
        ## fit a polynomial from a set of basis functions to estimate a N-dim curve

        #if waypts.shape[1] > 20:
                #waypts=np.delete(waypts, list(range(0, waypts.shape[1], 2)), axis=1)
        #if waypts.shape[1] > 20:
        #        waypts=np.delete(waypts, list(range(0, waypts.shape[1], 2)), axis=1)
        #if waypts.shape[1] > 10:
        #        waypts=np.delete(waypts, list(range(0, waypts.shape[1], 2)), axis=1)
        #if waypts.shape[1] > 10:
        #        waypts=np.delete(waypts, list(range(0, waypts.shape[1], 2)), axis=1)

        Ndim = waypts.shape[0]
        Nsamples = waypts.shape[1]

        #######################################################################
        ## discretization of trajectory
        #######################################################################
        M = 300 ## points on precomputed functions
        K = 500  ## number of precomputed basis functions
        plotFunctionalSpace = True
        #######################################################################

        if M < Nsamples:
                print "ERROR: more waypoints than discretization, abord"
                sys.exit(0)

        constraints = []
        print np.around(waypts,2)

        ##### FUNC SPACE CONSTRAINTS
        T = np.linspace(0.0,1.0,M)
        F = Fpoly(T,K)
        #F = Ffourier(T,K)
        dF = dFpoly(T,K)
        Weights = Variable(K,Ndim)

        #print np.around(F,decimals=2)
        #sys.exit(0)

        dw = 1.0/float(Nsamples-1)

        ctr=0
        Twpt = np.zeros((Nsamples,1))

        
        for i in range(0,Nsamples):
                tidx = find_nearest_idx(T,i*dw)
                Twpt[ctr]=tidx
                ctr=ctr+1
                Ftmp = np.reshape(F[:,tidx],(K,1))
                constraints.append(norm(waypts[:,i] - Weights.T*Ftmp) <= 0.05)
                #constraints.append(waypts[:,i] == Weights.T*Ftmp)

        ## add smoothing condition
        for t in T[1:]:
                tidx = find_nearest_idx(T,t)
                Ftmp0 = np.reshape(F[:,tidx-1],(K,1))
                Ftmp1 = np.reshape(F[:,tidx],(K,1))
                constraints.append(norm(Weights.T*Ftmp0 - Weights.T*Ftmp1) <= 0.05)


        if plotFunctionalSpace:
                plt.title('Waypoints')
                plt.subplot(3, 1, 1)
                plot(Twpt,waypts[0,:].flatten(),'ok',markersize=10)        
                plt.ylabel('X')
                plt.subplot(3, 1, 2)
                plot(Twpt,waypts[1,:].flatten(),'ok',linewidth=3,markersize=10)        
                plt.ylabel('Y')
                plt.subplot(3, 1, 3)
                plot(Twpt,waypts[2,:].flatten(),'ok',linewidth=3,markersize=10)
                plt.ylabel('Z')
                plt.show()

        objective = Minimize(norm(Weights,1))
        prob = Problem(objective, constraints)

        #ECOS, ECOS_BB, CVXOPT, SCS
        #result = prob.solve(solver=SCS, use_indirect=True, eps=1e-2, verbose=True)
        #prob.solve(verbose=True, abstol_inacc=1e-2,reltol_inacc=1e-2,max_iters= 300, reltol=1e-2)
        result = prob.solve(solver=SCS)

        if plotFunctionalSpace:
                Y = np.zeros((M,Ndim))
                ctr=0
                for t in T:
                        tidx = find_nearest_idx(T,t)
                        Ftmp = np.reshape(F[:,tidx],(K,1))
                        WF = Weights.T.value*Ftmp
                        Y[ctr,0] = WF[0]
                        Y[ctr,1] = WF[1]
                        Y[ctr,2] = WF[2]
                        ctr=ctr+1
                plt.title('Waypoints')
                plt.subplot(3, 1, 1)
                plot(Twpt,waypts[0,:].flatten(),'ok',markersize=10)        
                plot(Y[:,0].flatten(),'or',markersize=3)        
                plt.ylabel('X')
                plt.subplot(3, 1, 2)
                plot(Twpt,waypts[1,:].flatten(),'ok',linewidth=3,markersize=10)        
                plot(Y[:,1].flatten(),'or',linewidth=3,markersize=3)        
                plt.ylabel('Y')
                plt.subplot(3, 1, 3)
                plot(Twpt,waypts[2,:].flatten(),'ok',linewidth=3,markersize=10)
                plot(Y[:,2].flatten(),'or',linewidth=3,markersize=3)        
                plt.ylabel('Z')

                plt.show()

        if not (prob.status == OPTIMAL):
                print "ERROR: infeasible cvx program"
                sys.exit(0)

        return [Weights.value,T,F,dF]

class TrajectoryPolynomial(Trajectory):
        t_param = []
        F_basis = []
        dF_basis = []


        @classmethod
        def from_waypoints(cls, W):
                [W,t,F,dF] = WaypointsToWeights(W)
                cls.traj = W
                cls.t_param = t
                cls.F_basis = F
                cls.dF_basis = dF

                return cls(W)

        def new_from_waypoints(self, W):
                [W,t,F,dF] = WaypointsToWeights(W)
                self.traj = W
                self.t_param = t
                self.F_basis = F
                self.dF_basis = dF

        def evaluate_at(self, t):
                tidx = find_nearest_idx(self.t_param,t)
                W = self.traj
                K = self.F_basis.shape[0]
                Ftmp = self.F_basis[:,tidx]
                dFtmp = self.dF_basis[:,tidx]

                Ftmp = np.reshape(Ftmp,(K,1))
                dFtmp = np.reshape(dFtmp,(K,1))
                
                f0 = np.array((W.T[0,:]*Ftmp,W.T[1,:]*Ftmp,W.T[2,:]*Ftmp))
                df0 = np.array((W.T[0,:]*dFtmp,W.T[1,:]*dFtmp,W.T[2,:]*dFtmp))
                f0 = f0.flatten()
                df0 = df0.flatten()

                return [f0,df0]
