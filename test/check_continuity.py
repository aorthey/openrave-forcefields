
import numpy as np
import matplotlib.pyplot as plt

W = np.loadtxt('../W')

M = 100
x = W[0,0:M]
y = W[1,0:M]

for i in range(0,M):
        dx = W[:,i+1]-W[:,i]
        print dx[0]
plt.plot(x, y, '-or')
plt.show()
