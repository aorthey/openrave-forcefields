import sys
import numpy as np
fname = sys.argv[1]
print sys.argv
print fname
A=np.loadtxt(fname)
print "mean",np.mean(A[:,1])
print "std ",np.std(A[:,1])
