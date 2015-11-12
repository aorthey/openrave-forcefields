#!/usr/bin/env python
import time
import openravepy
from math import *

if not __openravepy_build_doc__:
    from openravepy import *
    from numpy import *

def GetForceArrowsInsideBox(W1,Fx,Fy):
    G1 = W1.GetGeometries()[0]
    B = G1.GetBoxExtents()
    T = G1.GetTransform()
    print B,T
    bx = T[0,3]
    by = T[1,3]
    N = int(B[0])
    M = int(B[1])
    dxspacing = B[0]/N
    dyspacing = B[1]/M

    dx = min(dxspacing-dxspacing/10,Fx)
    dy = min(dyspacing-dyspacing/10,Fy)
    if (dx == 0) & (dy==0):
            dx=dx+0.001

    l = sqrt(dx*dx+dy*dy)
    print dx,dy,N,M

    print dxspacing,dyspacing
    handles=[]
    xstart = -B[0]+bx
    ystart = -B[1]+by
    zval=0.1
    for i in range(0,2*N):
        for j in range(0,2*M):
            x = xstart+i*dxspacing+0.3*dxspacing
            y = ystart+j*dyspacing+0.3*dyspacing
            A = env.drawarrow(array((x,y,zval)),array((x+dx,y+dy,zval)),linewidth=0.08*l,color=array((1,0,0,0.5)))
            handles.append(A)
    return handles

def DrawReachableRegion(W1,Fx,Fy,ax):
    zval = 0.02
    G1 = W1.GetGeometries()[0]
    B = G1.GetBoxExtents()
    dx = B[0]
    dy = B[1]
    colorReachableRegion = array(((0,1,0,0.5)))
    if (Fx < ax) & (Fy < ax):
        handles.append(env.drawtrimesh(points=array(((-dx,-dy,zval),(-dx,dy,zval),(dx,dy,zval))),
                                     indices=None,
                                      colors=colorReachableRegion))
        handles.append(env.drawtrimesh(points=array(((-dx,-dy,zval),(dx,dy,zval),(dx,-dy,zval))),
                                     indices=None,
                                      colors=colorReachableRegion))
    else:
        ##compute cone for zero velocity
        print "cone"
        sx = 0.5*dt*dt*Fx
        sy = 0.5*dt*dt*Fy
        r = amax*dt*dt*0.5
        s = sqrt(sx*sx+sy*sy)
        theta = asin(r/s)
        print theta

        H = robot.GetTransform()
        rx = H[0,3]
        ry = H[1,3]
        px = tan(pi*0.5-theta)*dy + rx
        py = dy + ry
        handles.append(env.drawtrimesh(points=array(((rx,ry,zval),(px,py,zval),(dx,dy,zval))),
                                     indices=None,
                                      colors=colorReachableRegion))
        px = tan(pi*0.5-theta)*dy + rx
        py = -dy + ry
        handles.append(env.drawtrimesh(points=array(((rx,ry,zval),(px,py,zval),(dx,-dy,zval))),
                                     indices=None,
                                      colors=colorReachableRegion))
        handles.append(env.drawtrimesh(points=array(((rx,ry,zval),(dx,dy,zval),(dx,-dy,zval))),
                                     indices=None,
                                      colors=colorReachableRegion))
    return handles


if __name__ == "__main__":
    ###########################################################################
    ## environment setter
    env = Environment()
    env.SetViewer('qtcoin')
    env.Reset()

    ###########################################################################
    ## loader
    xmlenv='../../src/data/point_in_forcefield.env.xml'
    xmlrobot='../../src/robots/pointrobot.robot.xml'
    env.Add(env.ReadRobotXMLFile(xmlrobot))
    env.Load(xmlenv)
    ###########################################################################
    ## conventional variables
    robot = env.GetRobots()[0]
    I = eye(4)
    openravepy.misc.DrawAxes(env,I)
    ###########################################################################
    ## create force field
    physics = RaveCreatePhysicsEngine(env,'ode')
    physics.SetGravity(array((0,0,-9.81)))
    env.SetPhysicsEngine(physics)
    ###########################################################################

    handles=[]
    Fx=0.0
    Fy=0.0
    rx=0.0
    ry=0.0

    amax = 4.995

    ## assuming unit mass
    dt = 0.5
    vx = Fx*dt
    vy = Fy*dt
    zval = 0.05

    sx = 0.5*dt*dt*Fx
    sy = 0.5*dt*dt*Fy

    with env:
        #robot.GetLinks()[0].SetStatic(False)
        L = robot.GetLinks()[0]
        #physics.SetBodyForce(L,array((10.0,5.0,0.0)),array((1,1,0)),True)
        physics.SetLinkVelocities(L,(vx,vy,0),(0,0,0))
        r = amax*dt*dt*0.5
        handles.append(env.plot3(points=array(((rx+sx,ry+sy,0.05))),
                                   pointsize=r,
                                   colors=array(((0.2,0.2,0.2,0.2))),
                                   drawstyle=1))

        ## body 0 is robot, body 1 is env
        W = env.GetBodies()[1]
        W1 = W.GetLinks()[0]
        #handles.append(DrawReachableRegion(W1,Fx,Fy,amax))
        handles.append(GetForceArrowsInsideBox(W1,Fx,Fy))
        env.StopSimulation()
        #env.StartSimulation(timestep=0.001)
        starttime = time.time()

    #handles.append(env.plot3(points=array(((0.5,0,1.0),(-0.5,0,1.0))),
    #                   pointsize=45.0,
    #                   colors=array(((0,0,0,0.1),(0,0,0,0.8)))))
    # draw a random texture with alpha channel

    #while True:
    #    time.sleep(0.01)
    #    with env:
    #        curtime = time.time()
    #        T = curtime-starttime
    #        Hr = robot.GetTransform()
    #        xr = Hr[0,3]
    #        yr = Hr[1,3]
    #        Hr[2,3] = zval
    #        handles[0].SetTransform(Hr)
    #        handles[1].SetTransform(Hr)

    #with env:
        #L = robot.GetLinks()[1]
        #physics.SetBodyForce(L,array((-5.0,0.0,0.0)),array((0.0,0.0,0.0)),True)

    raw_input('Enter any key to quit. ')
