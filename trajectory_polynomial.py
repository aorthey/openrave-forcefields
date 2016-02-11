import cvxpy as cvx
from cvxpy import Variable, Problem, Minimize, norm, SCS
import numpy as np
from trajectory import *

def Fpoly(x,K): 
        return np.array(map(lambda e: map(lambda x: pow(x,e), x),np.arange(0,K)))

def dFpoly(x,K): 
        A= np.array(map(lambda e: map(lambda x: pow(x,e-1)/e, x),np.arange(1,K)))
        O = np.zeros((1,len(x)))
        return np.vstack((O,A))

def find_nearest(array, value):
        idx = (np.abs(array-value)).argmin()
        return array[idx]

def WaypointsToWeights(waypts):
        ## fit a polynomial from a set of basis functions to estimate a N-dim curve
        Ndim = waypts.shape[0]
        Nsamples = waypts.shape[1]

        ## discretization of trajectory
        M = 1000

        constraints = []
        objective = []
        K=Nsamples+100

        print waypts
        ##### FUNC SPACE CONSTRAINTS
        t = np.linspace(0.0,1.0,M)
        F = Fpoly(t,K)
        dF = dFpoly(t,K)

        Weights = Variable(K,Ndim)

        dw = 1.0/float(Nsamples-1)

        for i in range(0,Nsamples):
                tt = find_nearest(t,i*dw)
                print tt,waypts[:,i]
                constraints.append(waypts[:,i] == Weights.T*F[:,tt])

        objective = Minimize(norm(Weights,1))
        prob = Problem(objective, constraints)

        #prob.solve(solver=SCS, use_indirect=True, eps=1e-2, verbose=True)
        prob.solve(solver=SCS, eps=1e-3)
        return [Weights,t,F,dF]

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
                time_latticed = find_nearest(self.t_param,t)
                W = self.traj

                f0 = W.T*self.F_basis[:,time_latticed]
                df0 = W.T*self.dF_basis[:,time_latticed]

                f0 = np.array(f0.value).flatten()
                df0 = np.array(df0.value).flatten()

                print t,time_latticed,f0,df0

                return [f0,df0]
