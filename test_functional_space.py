import cvxpy as cvx
from cvxpy import Variable, Problem, Minimize, norm, SCS, OPTIMAL
from cvxpy import SCS, CVXOPT
from numpy import sqrt,sin,cos,pi,arccos
import numpy as np
import sys
from pylab import plot,title,xlabel,ylabel
import pylab as plt

## create a set of K orthonormal polynomial basis functions
def Fpoly(x,N): 
        ## assume x in [0,1]
        t = 2*x - 1.0
        ## Tn(x) = t^n
        return np.array(map(lambda e: map(lambda t: pow(t,e), t),np.arange(0,N)))

def dFpoly(x,N): 
        t = 2*x - 1.0
        A= np.array(map(lambda e: map(lambda t: pow(t,e-1)/e, t),np.arange(1,N)))
        O = np.zeros((1,len(t)))
        return np.vstack((O,A))

def Fchebyshev(x,K):
        t = 2*x - 1.0
        ## Tn(t) = cos( n*arccos(t) )
        return np.array( \
                        map(lambda n:  \
                                map(lambda t: cos(n*arccos(t)), t), \
                        np.arange(0,K)))

def dFchebyshev(x,K):
        t = 2*x - 1.0
        ## dTn/dt = n*Un-1(t)
        return np.array( \
                        map(lambda n:  \
                                map(lambda t: 
                                        n*((t+sqrt(t*t-1))**n - (t-sqrt(t*t-1))**n)/(2*sqrt(t*t-1)), \
                                t), \
                        np.arange(0,K))\
                )



def find_nearest(array, value):
        idx = (np.abs(array-value)).argmin()
        return array[idx]

def find_nearest_idx(array, value):
        idx = (np.abs(array-value)).argmin()
        return idx

def WaypointsToWeights(waypts):
        ## fit a polynomial from a set of basis functions to estimate a N-dim curve

        Ndim = waypts.shape[0]
        Nsamples = waypts.shape[1]

        #######################################################################
        ## discretization of trajectory
        #######################################################################
        M = 500 ## points on precomputed functions
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
        dF = dFpoly(T,K)
        #F = Fchebyshev(T,K)
        #dF = dFchebyshev(T,K)

        print np.around(F,decimals=2)
        Weights = Variable(K,Ndim)

        if plotFunctionalSpace:
                plt.title('Basis Functions')
                Kp = min(10,K)
                print T.shape,F.shape
                for i in range(0,Kp):
                        plt.subplot(Kp, 1, i)
                        plot(T,F[i,:],'-r',markersize=5)        
                        plt.ylabel(i)
                plt.show()
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
                constraints.append(norm(waypts[:,i] - Weights.T*Ftmp) <= 0.01)
                #constraints.append(waypts[:,i] == Weights.T*Ftmp)

        ## add smoothing condition
        for t in T[1:]:
                tidx = find_nearest_idx(T,t)
                Ftmp0 = np.reshape(F[:,tidx-1],(K,1))
                Ftmp1 = np.reshape(F[:,tidx],(K,1))
                constraints.append(norm(Weights.T*Ftmp0 - Weights.T*Ftmp1) <= 0.01)

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
        result = prob.solve(solver=SCS, verbose=True)

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

if __name__ == "__main__":
        M = 10
        W = np.random.rand(3,M)
        [W,t,F,dF] = WaypointsToWeights(W)
