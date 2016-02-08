
import numpy as np
from scipy.interpolate import interp1d,splev,splrep,splprep
from scipy.misc import derivative

class TrajectoryAnalyzer():
        tau = []
        handles = []
        degree = 3

        def __init__(self, W):
                if W.shape[1]<=4:
                        self.degree=1
                else:
                        self.degree=3

                self.tau,tmp = splprep(W,k=self.degree)

        def draw(self,env):
                N = 100
                for t in np.linspace(0.0, 1.0, num=N):
                        [f0,df0] = self.funcEval(t)
                        if f0.shape[0] == 2:
                                pts = np.array(((f0[0],f0[1],0.15)))
                        else:
                                pts = f0

        ptsize = 0.05
        def DrawRedPoint(self,env,X):
                self.handles.append(env.env.plot3(points=X,
                                   pointsize=self.ptsize,
                                   colors=np.array(((0.8,0.0,0.0,0.9))),
                                   drawstyle=1))
        def DrawGreenPoint(self,env,X):
                self.handles.append(env.env.plot3(points=X,
                                   pointsize=self.ptsize,
                                   colors=np.array(((0.0,0.8,0.0,0.9))),
                                   drawstyle=1))

        def funcEval(self, t):
                f0 = splev(t,self.tau)
                df0 = splev(t,self.tau,der=1)

                f0 = np.array(f0)
                df0 = np.array(df0)
                return [f0,df0]

        
        def analyze(self,env):
                N = 100
                for t in np.linspace(0.0, 1.0, num=N):
                        [f0,df0] = self.funcEval(t)
                        pts = np.array(((f0[0],f0[1],-0.1,0.001)))
                        F = env.GetForceAtX(pts)

                        d = np.linalg.norm(F)

                        pts = np.array(((f0[0],f0[1],0.1)))
                        if d > 0.2:
                                self.DrawRedPoint(env, pts)
                        else:
                                self.DrawGreenPoint(env, pts)




