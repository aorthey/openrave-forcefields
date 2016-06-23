import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PPoly
from scipy.interpolate import splev, splrep

W = np.array([-2.,-2.5 ,-3.7538003,-5.00760061,-5.50760061])
#W = np.hstack((W[0],W[:],W[-1]))

tvec = np.linspace(0,1,W.shape[0])
#tvec = np.hstack((-0.5,tvec,1.5))

tck = splrep(tvec, W, s=0, k=2)

x2 = np.linspace(0.0, 1.0, 100)
y2 = splev(x2, tck)

print splev(0,tck)
print splev(0,tck,der=1)
print splev(1,tck,der=1)

poly= PPoly.from_spline(tck)
dpoly = poly.derivative(1)
ddpoly = poly.derivative(2)

[B,idx]=np.unique(poly.x,return_index=True)
coeff = poly.c[:,idx]


plt.plot(tvec, W, 'o', x2, y2)
plt.show()



