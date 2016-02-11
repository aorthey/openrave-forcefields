from trajectory import *

class TrajectoryBSpline(Trajectory):
        @classmethod
        def from_waypoints(cls, W):
                Nwaypoints = W.shape[1]
                if Nwaypoints<=4:
                        degree=1
                else:
                        degree=3
                trajectory,tmp = splprep(W,k=degree)
                return cls(trajectory)

        def new_from_waypoints(self, W):
                Nwaypoints = W.shape[1]
                if Nwaypoints<=4:
                        degree=1
                else:
                        degree=3
                self.traj,tmp = splprep(W,k=degree)

        def evaluate_at(self, t):
                f0 = splev(t,self.traj)
                df0 = splev(t,self.traj,der=1)
                f0 = np.array(f0)
                df0 = np.array(df0)
                return [f0,df0]
