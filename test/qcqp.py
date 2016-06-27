import cvxopt
from cvxpy import *
import numpy as np

pnext =np.array( [-2.6181187748246475, -0.20745827989382473, 0.1000000008651387, -2.9842652455082757] )
A =np.array( [[-2000000.0, -0.0, -0.0, -0.0], [-0.0, -2000000.0, -0.0, -0.0], [-0.0, -0.0, -2000000.0, -0.0], [-0.0, -0.0, -0.0, -2000000.0], [2000000.0, 0.0, 0.0, 0.0], [0.0, 2000000.0, 0.0, 0.0], [0.0, 0.0, 2000000.0, 0.0], [0.0, 0.0, 0.0, 2000000.0]] )
b =np.array( [-5399638.220798268, -463140.95008524356, 200000.0017302774, -5919821.330636616, 5399635.7184592495, 463137.52396507084, -200000.0017302774, 5919820.330636616] )
ds= 0.01


Ndim = 4
Id = np.eye(Ndim)
x = Variable(Ndim)

objective = Minimize( norm(x - pnext)) #-c.T*x )
constraints = [ 0.5*quad_form(x,Id) <= ds*ds,
                A*x + b <= 0 ]

prob = Problem(objective, constraints)
res = prob.solve()
qcontrol  =np.array(x.value).flatten()

