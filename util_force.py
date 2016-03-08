import numpy as np
from util import *
from cvxpy import Variable, Problem, Minimize, norm, SCS, OPTIMAL
from cvxpy import SCS, CVXOPT
from numpy import sqrt,sin,cos,pi
from pylab import plot,title,xlabel,ylabel,subplot
import pylab as plt
import sys

inf = float('inf')

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
        plot([0,F[0]],[0,F[1]],'-g',linewidth=8,markersize=10)

        l = 0.1*sqrt(F[0]*F[0]+F[1]*F[1])
        hw = min(F[0]*l,F[1]*l)

        plot([p[0]],[p[1]],'or',markersize=10)
        plt.fill(X[:,0], X[:,1])
        plt.show()

## dF(q), W(q), W(q+1), dW(q)
def eval_norm(speed, dt, W, Wnext, dW, dF):
        dt2 = dt*dt/2
        Wstep = dt2*dF+speed*dW*dt+W
        d = np.linalg.norm(Wnext - Wstep)
        return [d,Wstep]

#dist_curve_wpt = get_minimal_distance_curve_point(speed, W, Wnext, dW, dF)
def get_minimal_distance_curve_point(speed, W, Wnext, dW, dF):
        DEBUG_TEST = 0

        dtold = 2000.0
        dt = 1000.0
        t = 0.0
        t_step = 0.005

        if DEBUG_TEST:
                S = []
                D = []

        Wstarr = np.zeros((4))
        while dt < dtold:
                dtold = dt
                [dt, Wstep] = eval_norm(speed, t, W, Wnext, dW, dF)
                t = t+t_step
                #Wstarr = np.vstack((Wstarr,Wstep))
                if DEBUG_TEST:
                        S.append(t)
                        D.append(dt)
                        print t,dt,dF


        #plot(Wstarr[:,0],Wstarr[:,1],'-',color=[0,0,1],linewidth=3)
        #plot([0,dF[0]],[0,dF[1]],'-',color=[1,0,0],linewidth=5)
        #plt.show()
        #sys.exit(0)

        if DEBUG_TEST:
                plot(S,D,'-r')
                plt.show()
                sys.exit(0)

        return dt

#rho[i] = get_minimal_epsilon_velocity_at_wpt(epsilon, W[:,i], W[:,i+1], dW[:,i], dist_w[i], dF[:,i])
def get_minimal_epsilon_velocity_at_wpt(epsilon, W, Wnext, dW, dist_next_wpt, dF):
        DEBUG = 0

        #dW = np.zeros(W.shape)
        #dW[0] = 1
        #Wnext = dist_next_wpt*dW
        dist_next_wpt = np.linalg.norm(W-Wnext)

        dold_speed = 2000.0
        dist_curve_wpt = 1000.0
        speed = 0.001
        speed_step = 0.01
        print "############"
        print np.linalg.norm(dF)
        if DEBUG:
                S = []
                D = []

        while dist_curve_wpt > epsilon:#d_speed < dold_speed:
                dold_speed = dist_curve_wpt
                Wstarr = W

                dist_curve_wpt = get_minimal_distance_curve_point(speed, W, Wnext, dW, dF)

                if dist_curve_wpt >= dist_next_wpt:
                        ### no progress, meaning the disturbance force moved us
                        ### into a direction opposite to the waypoint => still
                        ### increase speed to overcome the force
                        dold_speed = dist_curve_wpt + 1

                speed += speed_step
                if DEBUG:
                        S.append(speed)
                        D.append(dist_curve_wpt)

        if DEBUG:
                plot(S,D,'-r')
                plt.show()
        return speed-speed_step

def get_minimal_epsilon_velocity_curve(epsilon, W, dW, dF):
        Ndim = W.shape[0]
        Nsamples = W.shape[1]
        dist_w = np.zeros((Nsamples-1))
        rho = np.zeros((Nsamples))

        for i in range(0,Nsamples-1):
                dist_w[i] = np.linalg.norm(W[:,i]-W[:,i+1])

        for i in range(0,Nsamples-1):
                rho[i] = get_minimal_epsilon_velocity_at_wpt(epsilon, W[:,i], W[:,i+1], dW[:,i], dist_w[i], dF[:,i])

        rho[Nsamples-1]=0.0
        return rho



def GetMinimalSpeedToReachEpsilonNeighbordhoodVector(dt, epsilon, W, dW, dF):
        Ndim = W.shape[0]
        Nsamples = W.shape[1]

        dist_w = np.zeros((Nsamples-1))
        for i in range(0,Nsamples-1):
                dist_w[i] = np.linalg.norm(W[:,i]-W[:,i+1])

        p = Variable(Nsamples-1)
        sM = Variable(Nsamples-1)

        constraints = []
        objfunc = 0.0
        for i in range(0,Nsamples-1):
                #constraints.append( norm(p*dt*dW0[0:2] + dF +np.dot(dw,np.array((1,0))) ) < epsilon )
                constraints.append( norm(p[i]*dt*dW[:,i] + dt*dt/2*dF[:,i] + np.dot(dist_w[i],np.array((1,0,0,0))) ) < epsilon )
                constraints.append( sM[i] >= p[i] )
                constraints.append( sM[i] >= 0.0)
                constraints.append( p[i] >= 0.0 )
                objfunc += norm(sM[i])

        objective = Minimize(objfunc)

        prob = Problem(objective, constraints)

        print "solve minimal speed"
        result = prob.solve(solver=SCS)
        print "done.(",prob.value,"|",np.min(sM.value),")"

        if prob.value < inf:
                return np.array(sM.value).flatten()
        else:
                return np.array(sM.value).flatten()

def ForceCanBeCounteractedByAccelerationVector(dt, Fp, u1min, u1max, u2min, u2max, plot=False) :

        ### question (1) : can we produce an acceleration ddW, such that it counteracts F?

        ## dynamics projected onto identity element, it becomes obvious that in an infinitesimal neighborhood, 
        ## we can only counteract forces along the x and the theta axes due to non-holonomicity

        dt2 = dt*dt/2

        ## span dt2-hyperball in Ndim
        F = dt2*Fp
        thetamin = dt2*u2min
        thetamax = dt2*u2max
        xmin = 0.0
        xmax = dt2*u1max

        Xlow = np.dot(np.dot(Rz(-pi/2),Rz(thetamin)),np.array((1,0,0)))
        Xhigh = np.dot(np.dot(Rz(pi/2),Rz(thetamax)),np.array((1,0,0)))

        Ndim = Fp.shape[0]
        if Fp.ndim <= 1:
                Nsamples = 1
        else:
                Nsamples = Fp.shape[1]
        p = Variable(3,Nsamples)

        constraints = []
        objfunc = 0.0
        for i in range(0,Nsamples):
                constraints.append( norm(p[:,i]) <= xmax )
                constraints.append( np.matrix(Xlow[0:3])*p[:,i] <= 0 )
                constraints.append( np.matrix(Xhigh[0:3])*p[:,i] <= 0 )
                if Fp.ndim <= 1:
                        objfunc += norm(p[:,i]-F[0:3])
                else:
                        objfunc += norm(p[:,i]-F[0:3,i])
                #objfunc.append(norm(p[:,i]-F[:,i]))

        objective = Minimize(objfunc)
        prob = Problem(objective, constraints)

        result = prob.solve(solver=SCS, eps=1e-7)

        #nearest_ddq = np.array(p.value)
        nearest_ddq = np.array(p.value/dt2)

        codimension = Ndim-nearest_ddq.shape[0]

        #print Ndim, nearest_ddq.shape
        #print codimension
        zero_rows = np.zeros((codimension,Nsamples))

        if nearest_ddq.shape[0] < Ndim:
                nearest_ddq = np.vstack((nearest_ddq,zero_rows))

        if plot:
                PlotReachableSetForceDistance(dt, u1min, u1max, u2min, u2max, -F, dt2*nearest_ddq)

        return nearest_ddq

if __name__ == '__main__':
        F = np.array((1,-0.05,0.0))

        u1min = 0.0
        u1max = 5.0
        u2min = -5.0
        u2max = 5.0

        dt = 0.05
        ForceCanBeCounteractedByAccelerationVector(dt, F, u1min, u1max, u2min, u2max, plot=True)
