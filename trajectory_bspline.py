from trajectory import *
import sys

DEBUG=1

class TrajectoryBSpline(Trajectory):
        SMOOTH_CONSTANT=0.0
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
                trajectory,tmp = splprep(W,k=cls.POLYNOMIAL_DEGREE,s=cls.SMOOTH_CONSTANT)
                return cls(trajectory)

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
