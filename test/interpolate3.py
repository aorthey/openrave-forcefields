import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PPoly
from scipy.interpolate import splev, splrep

a = -2.0 
b = -0.376331377066 
c = 0.921500133388 
d = -0.553193147758
a = -2.0 
b = -1.0 
c = 2.99140401146 
d = -1.99426934097 1.0 -2.00286532951 -2.00286532951

T = 1.0

M = 100
tvec = np.linspace(0.0,T,M)
X = np.zeros((M))
for i in range(0,M):
        t = tvec[i]
        X[i] = a + t*b + t*t*c + t*t*t*d
plt.plot(tvec,X, '-or')
plt.show()

