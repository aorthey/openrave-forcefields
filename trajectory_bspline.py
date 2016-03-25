from trajectory import *

DEBUG=1


class TrajectoryBSpline(Trajectory):
        SMOOTH_CONSTANT=0.01
        @classmethod
        def from_waypoints(cls, W):
                Nwaypoints = W.shape[1]
                if Nwaypoints<=4:
                        degree=2
                else:
                        degree=5
                print Nwaypoints,degree
                trajectory,tmp = splprep(W,k=degree,s=cls.SMOOTH_CONSTANT)
                return cls(trajectory)

        def new_from_waypoints(self, W):
                #print "BEFORE WAYPOINTS:"
                #print np.around(W[:,0],decimals=2), np.around(W[:,-1],decimals=2)
                Nwaypoints = W.shape[1]
                if Nwaypoints<=4:
                        degree=2
                else:
                        degree=5
                self.traj,tmp = splprep(W,k=degree,s=self.SMOOTH_CONSTANT)
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
