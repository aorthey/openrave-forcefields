import numpy as np
import matplotlib.pyplot as plt

x= np.linspace(0,1,100)
y = np.sqrt(x)
plt.plot(x, y, '-or')
y = -1.2*(1-x)**3+1.0
plt.plot(x, y, '-ob')
plt.show()
