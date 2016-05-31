import numpy as np
from math import pi,cos,sin
from pylab import *
from matplotlib import gridspec
import pylab as plt

M=100
radius = 1
xoffset=0.15
yoffset=0.0
fs=22
ylmin=-1.2
ylmax=1.2
xlmin=-1.2
xlmax=1.2
plt.xlim([xlmin,xlmax])
T = np.linspace(-pi/2-0.5,pi/2+0.5,M)
r = np.array((radius*np.cos(T),radius*np.sin(T),np.zeros(M)))
F = np.array((0,1,0))

fig=figure(facecolor='white', figsize=(10,6))

gs = gridspec.GridSpec(2, 3)
###############################################################################
ax1 = plt.subplot(gs[0,:])
###############################################################################
torque = np.array(map( lambda i: np.cross(r[:,i],F), range(0,r.shape[1])))
plt.plot(T,torque[:,2],'r',linewidth=3)
plt.plot(-pi/2,0,'ok',markersize=10)
plt.plot(0,radius,'ok',markersize=10)
plt.plot(pi/2,0,'ok',markersize=10)
ax1.text(-pi/2+xoffset, 0+yoffset, '$1$', fontsize=fs)
ax1.text(0, radius-xoffset, '$2$', fontsize=fs)
ax1.text(pi/2+xoffset, 0+yoffset, '$3$', fontsize=fs)
xlabel('Angle (radians)', fontsize=fs)
ylabel('Torque (Nm)', fontsize=fs)
###############################################################################
ax2 = plt.subplot(gs[1,0])
ms=30
lw=3
xoffset=-0.2
yoffset=-0.05
###############################################################################
ylabel('Y (m)', fontsize=fs)
xlabel('X (m)', fontsize=fs)
ax2.text(xoffset, yoffset, '$1$', fontsize=fs)
ax2.plot(0,0,'ok',markersize=ms)
ax2.plot([0,0],[1,0],'-k',linewidth=lw)
ax2.plot([0,1],[0,0],'--k',linewidth=lw)
plt.ylim([ylmin,ylmax])
plt.xlim([xlmin,xlmax])
###############################################################################
ax3 = plt.subplot(gs[1,1])
ax3.set_yticklabels(())
xlabel('X (m)', fontsize=fs)
plt.text(xoffset, yoffset, '$2$', fontsize=fs)
plt.plot(0,0,'ok',markersize=ms)
plt.plot([0,1],[0,0],'-k',linewidth=lw)
plt.plot([0,1],[0,0],'--k',linewidth=lw)
plt.ylim([ylmin,ylmax])
plt.xlim([xlmin,xlmax])
###############################################################################
ax3.plot(0,0,'ok',markersize=20)
###############################################################################
ax4 = plt.subplot(gs[1,2])
ax4.set_yticklabels(())
xlabel('X (m)', fontsize=fs)
plt.text(xoffset, yoffset, '$3$', fontsize=fs)
plt.plot(0,0,'ok',markersize=ms)
plt.plot([0,0],[0,-1],'-k',linewidth=lw)
plt.plot([0,1],[0,0],'--k',linewidth=lw)
plt.ylim([ylmin,ylmax])
plt.xlim([xlmin,xlmax])
###############################################################################
ax4.plot(0,0,'ok',markersize=20)

plt.show()
