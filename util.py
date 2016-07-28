import numpy as np
import math
from math import pi,cos,sin,acos,asin,atan

def draw_waypoints(W, env, ptsize=0.05, trajectory_color = np.array((0,1,0,1))):
        handle = []
        handle.append(env.env.drawlinestrip(points=W.T,
                                   linewidth=10,
                                   colors=trajectory_color))
        return handle

def draw_ravetraj(ravetraj, env, ptsize=0.05, trajectory_color = np.array((0.1,0.8,0.1,1))):
        N = ravetraj.GetNumWaypoints()
        W=[]
        for i in range(0,N):
                w = np.array((ravetraj.GetWaypoint(i)[0],ravetraj.GetWaypoint(i)[1],ravetraj.GetWaypoint(i)[2]))
                W.append((w))
        W = np.array(W).T
        return draw_waypoints(W,env,ptsize,trajectory_color)

def rand_interval( a, b):
        return rand_interval_NM(a,b,1,1)

def rand_interval_NM( a, b, n, m):
        #return r \in [a,b]
        assert(b>=a)
        r = np.random.rand(n,m)
        rab = r*(b-a)+a
        return rab


inf = float('inf')
black=np.array((0,0,0))
white=np.array((1,1,1))
red=np.array((1,0,0))
green=np.array((0,1,0))
blue=np.array((0,0,1))
xe = np.array((1,0,0))
ye = np.array((0,1,0))
ze = np.array((0,0,1))
ex = np.array((1,0,0))
ey = np.array((0,1,0))
ez = np.array((0,0,1))
exh = np.array((1,0,0,1))
eyh = np.array((0,1,0,1))
ezh = np.array((0,0,1,1))

def Rz(t):
        return np.array([[cos(t),-sin(t),0],[sin(t),cos(t),0],[0,0,1]])
def Ry(t):
        return np.array([[cos(t),0,sin(t)],[0,1,0],[-sin(t),0,cos(t)]])
def Rx(t):
        return np.array([[1,0,0],[0,cos(t),-sin(t)],[0,sin(t),cos(t)]])

def Rax(theta, axis):
#def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    theta = np.asarray(theta)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

#def Rax(theta, u):
#    return [[cos(theta) + u[0]**2 * (1-cos(theta)), 
#             u[0] * u[1] * (1-cos(theta)) - u[2] * sin(theta), 
#             u[0] * u[2] * (1 - cos(theta)) + u[1] * sin(theta)],
#            [u[0] * u[1] * (1-cos(theta)) - u[2] * sin(theta),
#             cos(theta) + u[1]**2 * (1-cos(theta)),
#             u[1] * u[2] * (1 - cos(theta)) + u[0] * sin(theta)],
#            [u[0] * u[2] * (1-cos(theta)) - u[1] * sin(theta),
#             u[1] * u[2] * (1-cos(theta)) - u[0] * sin(theta),
#             cos(theta) + u[2]**2 * (1-cos(theta))]]

def HTT(T):
        I=np.identity(3)
        return HT(I,T)

def HTR(R):
        return HT(R,np.array((0,0,0)))

def HT(R,T):
        H1 = np.vstack((R,[0,0,0]))
        H1 = np.hstack((H1,[[T[0]],[T[1]],[T[2]],[1]]))
        return H1

def getZYsphericalRot(ap, xl):
        pxy = ap-np.dot(ap,ze)*ze
        pxy = pxy/np.linalg.norm(pxy)
        pzx = ap-np.dot(ap,ye)*ye
        pzx = pzx/np.linalg.norm(pzx)
        xln = xl/np.linalg.norm(xl)

        theta = acos( np.dot(pxy,xln) )
        phi = acos( np.dot(pzx,xln) )
        if np.dot(ap,ze) < 0:
                phi = -phi
        if np.dot(ap,ye) > 0:
                theta = -theta
        return [theta,phi]


def getSphericalCoordinatesFromCart(x,y,z):
        ### compute azimuth/inclination from dtau
        r = math.sqrt(x*x + y*y + z*z)
        phi = math.acos(z/r)
        theta = math.atan2(y,x)

        ## compute inverse
        xx = r*sin(phi)*cos(theta)
        yy = r*sin(phi)*sin(theta)
        zz = r*cos(phi)

        return [theta,phi]

def getGlobalTransformation(tau, dtau):
        spheretheta,spherephi = getSphericalCoordinatesFromCart(dtau[0],dtau[1],dtau[2])
        Rglob = np.dot(Rz(spheretheta),Ry(-(math.pi/2-spherephi)))
        Hglob = HT(Rglob,tau)
        return Hglob


def PrintNumpy(name, p, ReturnString=False):
        pstr = str(name+"=np.array("+"".join(str(p.tolist()))+")")
        if ReturnString:
                return pstr
        else:
                print pstr
