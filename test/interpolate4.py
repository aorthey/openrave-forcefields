import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PPoly
from scipy.interpolate import splev, splrep

def SimpleInterpolate(q0,q1,qd0,qd1,T):
        a=((qd1-qd0)*T-2*(q1-q0-qd0*T))/T**3
        b=(3*(q1-q0-qd0*T)-(qd1-qd0)*T)/T**2
        c=qd0
        d=q0
        return [d,c,b,a]

T = 0.01
q0 = np.array([-2,0])
q1 = np.array([-2.0294,-2.3622731e-05])
qd0 = np.array([-1,-8.03e-04])
qd1 = np.array([-1,-8.03e-04])
[a,b,c,d] = SimpleInterpolate(q0,q1,qd0,qd1,T)

M = 100
tvec = np.linspace(0.0,T,M)
X = np.zeros((2,M))
for i in range(0,M):
        t = tvec[i]
        X[:,i] = a + t*b + t*t*c + t*t*t*d

plt.plot(1.0,1.0,'ok')
plt.plot(tvec,X[0,:], '-or')
plt.show()

