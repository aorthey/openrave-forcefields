from trajectory import *

DEBUG=1


class TrajectoryBSpline(Trajectory):
        SMOOTH_CONSTANT=0.2
        @classmethod
        def from_waypoints(cls, W):
                Nwaypoints = W.shape[1]
                if Nwaypoints<=4:
                        degree=1
                else:
                        degree=5
                trajectory,tmp = splprep(W,k=degree,s=cls.SMOOTH_CONSTANT)
                return cls(trajectory)

        def new_from_waypoints(self, W):
                #print "BEFORE WAYPOINTS:"
                #print np.around(W[:,0],decimals=2), np.around(W[:,-1],decimals=2)
                Nwaypoints = W.shape[1]
                if Nwaypoints<=4:
                        degree=1
                else:
                        degree=5
                self.traj,tmp = splprep(W,k=degree,s=self.SMOOTH_CONSTANT)
                [W0, tmp] = self.evaluate_at(0)
                [W1, tmp] = self.evaluate_at(1)
                #print "AFTER WAYPOINTS:"
                #print np.around(W0,decimals=2), np.around(W1,decimals=2)

        def evaluate_at(self, t):
                f0 = splev(t,self.traj)
                df0 = splev(t,self.traj,der=1)
                f0 = np.array(f0)
                df0 = np.array(df0)
                return [f0,df0]
