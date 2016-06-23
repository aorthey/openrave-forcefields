import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PPoly
from scipy.interpolate import splev, splrep

def RepairTrajectory(W,max_dist):
        ictr=0
        i=0
        while i < W.shape[1]-1:

                d = np.linalg.norm(W[:,i+1]-W[:,i])

                istep = 1
                if d > max_dist:

                        dstep = max_dist
                        Nnew = int(d/dstep)

                        dcur = 0.0
                        print i,d,"/",max_dist,"d/dstep",d/dstep,Nnew
                        for k in range(1,Nnew):
                                Wstep = (W[:,i+k]-W[:,i])
                                Wstep /= np.linalg.norm(Wstep)
                                Wnew = W[:,i] + k*dstep*Wstep
                                W= np.insert(W, i+k, Wnew, axis=1)
                                ictr+=1 
                        istep = Nnew

                i = i+istep
        print "new points:",ictr
        return W

Win = np.loadtxt('../W')

W = RepairTrajectory(Win, 1e-4)

kp = 3
tvec = np.linspace(0,1,W.shape[1])
trajx = splrep(tvec, W[0,:], s=0, k=kp)
trajy = splrep(tvec, W[1,:], s=0, k=kp)

t = np.linspace(0.0, 1.0, 10000)
x = splev(t, trajx)
y = splev(t, trajy)

#M=100
#W = W[:,0:M]

plt.plot(Win[0,:], Win[1,:], '-or')
plt.plot(W[0,:], W[1,:], 'ok')
plt.plot(x,y,'-b')
plt.show()

