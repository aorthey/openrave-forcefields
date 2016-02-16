import numpy as np
from util import *
from cvxpy import Variable, Problem, Minimize, norm, SCS, OPTIMAL
from cvxpy import SCS, CVXOPT
from numpy import sqrt,sin,cos,pi
from pylab import plot,title,xlabel,ylabel,subplot
import pylab as plt


def PlotReachableSetForceDistance(dt, u1min, u1max, u2min, u2max, F, p):

        dt2 = dt*dt/2
        M = 200
        X = np.zeros((M,2))
        u2 = np.linspace(u2min,u2max,M-1)
        for i in range(0,M-1):
                r = u1max*dt2
                theta = u2[i]*dt2
                X[i,0] = r*cos(theta)
                X[i,1] = r*sin(theta)

        plot(X[:,0],X[:,1],'ob')
        plot([0,F[0]],[0,F[1]],'-g',linewidth=10,markersize=10)
        plot([p[0]],[p[1]],'or',markersize=10)
        plt.fill(X[:,0], X[:,1])
        plt.show()

def GetMinimalSpeedToReachEpsilonNeighbordhood(dt, epsilon, W0, W1, dW0, dF):
        dw = np.linalg.norm(W1-W0)

        p = Variable(1)
        sM = Variable(1)

        constraints = []
        constraints.append( norm(p*dt*dW0[0:2] + dF +np.dot(dw,np.array((1,0))) ) < epsilon )
        constraints.append( sM > p )
        constraints.append( sM > 0 )
        constraints.append( p > 0 )

        objective = Minimize(norm(sM))

        prob = Problem(objective, constraints)
        result = prob.solve(solver=SCS)

        return sM.value


def ForceCanBeCounteractedByAcceleration(dt, Fp, u1min, u1max, u2min, u2max, plot=False) :

        ### question (1) : can we produce an acceleration ddW, such that it counteracts F?

        ## dynamics projected onto identity element, it becomes obvious that in an infinitesimal neighborhood, 
        ## we can only counteract forces along the x and the theta axes due to non-holonomicity

        X1 = np.array((1,0,0,0))
        X2 = np.array((0,0,0,1))

        dt2 = dt*dt/2

        ## span dt2-hyperball in Ndim
        F = dt2*Fp
        thetamin = dt2*u2min
        thetamax = dt2*u2max
        xmin = 0.0
        xmax = dt2*u1max

        Xlow = np.dot(np.dot(Rz(-pi/2),Rz(thetamin)),np.array((1,0,0)))
        Xhigh = np.dot(np.dot(Rz(pi/2),Rz(thetamax)),np.array((1,0,0)))

        Ndim = 2
        p = Variable(Ndim,1)

        constraints = []
        constraints.append( norm(p) <= xmax )
        constraints.append( np.matrix(Xlow[0:2])*p <= 0 )
        constraints.append( np.matrix(Xhigh[0:2])*p <= 0 )

        objective = Minimize(norm(p - F))

        prob = Problem(objective, constraints)
        result = prob.solve(solver=SCS)
        nearest_ddq = np.array(p.value)

        d = np.linalg.norm(nearest_ddq)

        if plot:
                PlotReachableSetForceDistance(dt, u1min, u1max, u2min, u2max, F, nearest_ddq)
        return [d,nearest_ddq]

if __name__ == '__main__':
        F = np.array((2.5,0.01))

        u1min = 0.0
        u1max = 2.0
        u2min = -5.0
        u2max = 5.0

        dt = 0.05
        ForceCanBeCounteractedByAcceleration(dt, F, u1min, u1max, u2min, u2max, plot=True)
