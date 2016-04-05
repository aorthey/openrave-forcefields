from trajectory import *
from scipy.interpolate import PPoly
import sys

DEBUG=1

class TrajectoryBSpline(Trajectory):
        SMOOTH_CONSTANT=0
        POLYNOMIAL_DEGREE=3
        MIN_NUMBER_WAYPOINTS = 5

        @classmethod
        def addMinimalWaypoints(cls, W):
                Nwaypoints = W.shape[1]
                if Nwaypoints >= cls.MIN_NUMBER_WAYPOINTS:
                        return W

                Wtmp = W
                for i in range(0,Nwaypoints-1):
                        Wnew = W[:,i]+0.5*(W[:,i+1]-W[:,i])
                        Wtmp = np.insert(Wtmp, 2*i+1, values=Wnew, axis=1)
                W = Wtmp
                return cls.addMinimalWaypoints(W)

        @classmethod
        def prettify(cls, W):
                Nwaypoints = W.shape[1]
                if Nwaypoints <= 1:
                        print "cannot create trajectory with only one waypoint"
                        sys.exit(1)
                if Nwaypoints <= 3:
                        W = cls.addMinimalWaypoints(W)
                        Nwaypoints = W.shape[1]
                return W


        @classmethod
        def from_waypoints(cls, W):
                print W
                W = cls.prettify(W)
                print W
                #trajectory,tmp = splprep(W,k=cls.POLYNOMIAL_DEGREE,s=cls.SMOOTH_CONSTANT)
                [T,D] = cls.computeTrajectoryStringForTOPP(W)

                self.trajectorystring_topp = T
                self.durationvector_topp = D
                print T,D
                sys.exit(0)

                return cls(trajectory)

        @classmethod
        def computeTrajectoryStringForTOPP(cls, W):
                Ndim = W.shape[0]
                Nwaypoints = W.shape[1]


                ###############################################################
                ##### get number of intervals between breakpoints
                tvec = np.linspace(0,1,W.shape[1])
                trajectory = splrep(tvec,W[0,:],k=cls.POLYNOMIAL_DEGREE,s=cls.SMOOTH_CONSTANT)
                poly= PPoly.from_spline(trajectory)
                [B,idx]=np.unique(poly.x,return_index=True)
                Ninterval = B.shape[0]-1
                print Ninterval,Ndim
                ###############################################################
                P = np.zeros((Ninterval, Ndim, 4))

                durationVector = np.zeros((Ninterval))
                for j in range(1,B.shape[0]):
                        d = B[j]-B[j-1]
                        durationVector[j-1] = d
                print durationVector

                for i in range(0,Ndim):

                        tvec = np.linspace(0,1,Nwaypoints)
                        trajectory = splrep(tvec,W[i,:],k=cls.POLYNOMIAL_DEGREE,s=cls.SMOOTH_CONSTANT)

                        poly= PPoly.from_spline(trajectory)
                        dpoly = poly.derivative(1)
                        ddpoly = poly.derivative(2)
                        [B,idx]=np.unique(poly.x,return_index=True)
                        coeff = poly.c[:,idx]

                        def eval_coeff(t,coeff):
                                return t*t*t*coeff[0] + t*t*coeff[1] + t*coeff[2] + coeff[3]
                        def deval_coeff(t,coeff):
                                return 3*t*t*coeff[0] + 2*t*coeff[1] + coeff[2]
                        def ddeval_coeff(t,coeff):
                                return 6*t*coeff[0] + 2*coeff[1]
                        def g(dt,A,B,C,D):
                                return A + B*dt + C*dt**2 + D*dt**3
                        def dg(dt,A,B,C,D):
                                return B + 2*C*dt + 3*D*dt**2
                        def ddg(dt,A,B,C,D):
                                return 2*C + 6*D*dt

                        print coeff[3,:-1]
                        P[:,i,0] = coeff[3,:-1]
                        P[:,i,1] = coeff[2,:-1]
                        P[:,i,2] = coeff[1,:-1]
                        P[:,i,3] = coeff[0,:-1]

                        #Ninterval = B.shape[0]-1
                        #Ai = np.zeros((Ninterval))
                        #Bi = np.zeros((Ninterval))
                        #Ci = np.zeros((Ninterval))
                        #Di = np.zeros((Ninterval))
                        #for j in range(Ninterval):
                        #        t = B[j]
                        #        print "poly in:", B[j],B[j+1]
                        #        tstep = 0.1
                        #        Ai[j] = coeff[3,j]
                        #        Bi[j] = coeff[2,j]
                        #        Ci[j] = coeff[1,j]
                        #        Di[j] = coeff[0,j]
                        #        while t <= B[j+1]+1e-5:
                        #                tin = t-B[j]
                        #                print t,"/",B[j+1],eval_coeff(tin, coeff[:,j]), deval_coeff(tin, coeff[:,j]), ddeval_coeff(tin, coeff[:,j])
                        #                t=t+tstep
                        #                print g(tin,Ai[j],Bi[j],Ci[j],Di[j])
                        #        print coeff[:,j]
                for i in range(0,durationVector.shape[0]):
                        duration = durationVector[i]
                        if i==0:
                                trajectorystring = str(duration)
                        else:
                                trajectorystring += "\n" + str(duration)
                        trajectorystring += "\n" + str(Ndim)

                        for j in range(Ndim):
                                trajectorystring += "\n"
                                trajectorystring += string.join([str(P[i,j,0]),str(P[i,j,1]),str(P[i,j,2]),str(P[i,j,3])])

                return [trajectorystring, durationVector]


        def new_from_waypoints(self, W):
                if W.shape[1]<5:
                        print "number of waypoints diminished below",self.MIN_NUMBER_WAYPOINTS
                        print "wpts:",W.shape[1]
                        print W
                        sys.exit(1)
                self.traj,tmp = splprep(W,k=self.POLYNOMIAL_DEGREE,s=self.SMOOTH_CONSTANT)
                [W0, tmp] = self.evaluate_at(0)
                [W1, tmp] = self.evaluate_at(1)
                #print "AFTER WAYPOINTS:"
                #print np.around(W0,decimals=2), np.around(W1,decimals=2)

        def evaluate_at(self, t, der=1):
                f = splev(t,self.traj)
                df = splev(t,self.traj,der=1)
                f = np.array(f)
                df = np.array(df)
                df[2]=0.0

                #ndf = np.linalg.norm(df)
                #df=df/ndf

                from util import Rz
                if der>1:
                        ddf = splev(t,self.traj,der=2)
                        ddf = np.array(ddf)
                        ddf[2]=0.0
                        #ddf = np.array(ddf0)
                        #ddf = np.zeros((f.shape))
                        #ddf[0] = df[1]
                        #ddf[1] = -df[0]
                        ### remove component along df0 from ddf0
                        #ddf[2]=0.0
                        #a = np.dot(ddf,df)
                        #ddf = ddf - a*df

                        #nddf = np.linalg.norm(ddf)
                        #ddf=ddf/nddf
                        return [f,df,ddf]
                return [f,df]
