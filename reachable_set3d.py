from shapely.ops import cascaded_union, polygonize
import sys
import re
import os
from svglib.svglib import svg2rlg
from reportlab.graphics import renderPDF
from os.path import basename

import subprocess 
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import matplotlib.ticker as mtick
import matplotlib
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon
from scipy.spatial import Delaunay
import shapely.geometry as geometry
import pylab as plt
import numpy as np
import math
from util import *
import parameters_dynamical_system as params
class Arrow3D(FancyArrowPatch):
        def __init__(self, xs, ys, zs, *args, **kwargs):
                FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
                self._verts3d = xs, ys, zs

        def draw(self, renderer):
                xs3d, ys3d, zs3d = self._verts3d
                xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
                self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
                FancyArrowPatch.draw(self, renderer)

class ReachableSet3D():

        pts = None
        poly = []
        rs_boundary_thickness = 4
        loc = 'upper left'
        qhull_options = 'QJ'
        #loc = 'best'

        rs_color = np.array((0.5,0.5,1.0,0.2))
        rs_last_color = np.array((0.2,0.2,1.0,0.3))
        rs_last_edge_color = np.array((0.8,0.8,0.8,0.8))
        rs_alpha = 0.8
        rs_last_alpha = 0.9

        point_size = 200
        path_lw = 3

        force_color = np.array((1.0,0,0))
        tangent_color = np.array((0,0,0))
        orientation_color = np.array((0.9,0.0,0.5))
        fs = 28
        fs_label = 34
        fs_title = 46
        image = None

        def ndmesh(*args):
                args = map(np.asarray,args)
                return np.broadcast_arrays(*[x[(slice(None),)+(None,)*i]
                      for i, x in enumerate(args)])

        def __init__(self, ds, p, s, dp, force, R, amin, amax):
                
                self.pts = None
                self.poly = []
                self.fig = plt.figure(facecolor='white')
                #self.fig.set_size_inches(18.5, 18.5, forward=True)
                from mpl_toolkits.mplot3d import Axes3D
                self.image = self.fig.gca(projection='3d')

                tsamples= 15
                [qnext,dtmp,utmp,tend] = params.ForwardSimulate(p, dp, s, ds, force)
                print "TEND:",tend,qnext
                tstart = 0.0#tend/100

                Ndim = p.shape[0]
                poly = []
                
                ### A all possible control input combination
                A = np.vstack((amin,amax)).T
                A = np.array(np.meshgrid(*A)).T.reshape(-1,amin.shape[0])
                M_time = A.shape[0]

                for dt in self.expspace(tstart,tend,tsamples):
                        dt2 = dt*dt*0.5
                        M_speed = 1
                        speed = np.linspace(0,s,M_speed)
                        speed=[s]
                        q = np.zeros((Ndim,M_speed*M_time))
                        for k in range(0,M_speed):
                                for i in range(0,M_time):
                                        control = np.dot(R,A[i])
                                        q[:,i+k*M_time] = p + dt*speed[k]*dp + dt2 * force + dt2 * control

                        self.add_points( q.T )

                dt = tend
                dt2 = dt*dt*0.5

                qnext = p+dt*s*dp+dt2*force

                pnext = p+ds*dp/np.linalg.norm(dp)

                tstring = 'Reachable Set (T<='+str("%10.2e"%dt)+')'
                self.filename = 'images/reachableset_'+params.FILENAME+'_ori'+str(np.around(p[3],decimals=2))
                self.filename = re.sub('[.]', '-', self.filename)
                self.filename += '.png'

                self.image.set_title(tstring, fontsize=self.fs_title, y=1.1)
                self.image.set_xlabel('\n\nX-Position [m]', fontsize=self.fs_label)
                self.image.set_ylabel('\n\nY-Position [m]', fontsize=self.fs_label)
                self.image.set_zlabel('\n\n$\\theta$-Position [rad]', fontsize=self.fs_label)

                self.arrow_head_size = np.linalg.norm(dt2*force)/5
                self.hw = 0.5*self.arrow_head_size
                self.lw = 0.2*self.arrow_head_size

                self.image.scatter(p[0], p[1], p[3], color='k',
                                s=self.point_size)

                self.image.scatter(pnext[0],pnext[1],pnext[3], 'og', color='g',
                                s=self.point_size)
                self.image.scatter(qnext[0],qnext[1],qnext[3], 'ok',
                                s=self.point_size)

                [qcontrol,qtmp,utmp,dtmp] = params.ForwardSimulate(p, dp, s, ds, force)
                print "TEND:",dtmp,qcontrol

                dq = qcontrol - pnext
                dnp = dp/np.linalg.norm(dp)
                wmove = dq - np.dot(dq,dnp)*dnp

                self.PlotPathSegment( pnext, pnext+wmove)

                ds = np.linalg.norm(p-pnext)

                self.image.scatter(qcontrol[0],qcontrol[1],qcontrol[3], 'or',
                                color='r',
                                s=2*self.point_size)

                text_offset_t = ds/10.0
                text_offset_y = ds/10.0
                text_offset_x = ds/10.0

                self.image.text(p[0]-text_offset_x, p[1]-text_offset_y, p[3]-text_offset_t, 
                                "$p(s)$",
                                color='black',
                                fontsize=self.fs)
                self.image.text(pnext[0]-text_offset_x, pnext[1]-text_offset_y, pnext[3]-text_offset_t,
                                "$p(s+\\Delta s)$", 
                                color='black',
                                fontsize=self.fs)

                dori = np.zeros((Ndim))
                dori[0:2] = np.dot(Rz(p[3]),ex)[0:2]
                dori[0:2] = 0.5*(dori[0:2]/np.linalg.norm(dori[0:2]))
                dori[3] = p[3]

                pathnext = p + 1.5*ds*dp/np.linalg.norm(dp)
                pathlast = p - 0.2*ds*dp/np.linalg.norm(dp)

                self.PlotPathSegment( pathlast, pathnext)

                #arrow0 = self.PlotArrow(p, dt*s*dori, self.orientation_color)
                arrow1 = self.PlotArrow(p+dt2*force, dt*s*dp, self.tangent_color, ls='dashed')
                arrow1 = self.PlotArrow(p, dt*s*dp, self.tangent_color)
                arrow2 = self.PlotArrow(p, dt2*force, self.force_color)

                #arrow0 = self.PlotArrow(qnext, dt*s*dori, self.orientation_color)
                arrow2 = self.PlotArrow(p+dt*s*dp, dt2*force, self.force_color)

                plt.legend([arrow1,arrow2,],
                                ['Velocity Displacement','Wrench Displacement',],
                                fontsize=self.fs,
                                loc=self.loc)
                self.origin = p
                self.tangent = dp

                #draw sphere
                ts1 = acos( np.dot(dp[0:2],ex[0:2])/np.linalg.norm(dp[0:2]))
                ts2 = acos( np.dot((qnext-p)[0:2],ex[0:2])/np.linalg.norm((qnext-p)[0:2]))

                toffset = pi/8
                toffsetz = pi/4
                ori1 = np.cross(ex[0:2],(dp)[0:2])
                ori2 = np.cross(ex[0:2],(qnext-p)[0:2])
                if ori1 < 0:
                        ts1 *= -1
                if ori2 < 0:
                        ts2 *= -1

                if ts1 > ts2:
                        tlimU = ts1+toffset
                        tlimL = ts2-2*toffset
                else:
                        tlimL = ts1-toffset
                        tlimU = ts2+2*toffset

                u, v = np.mgrid[tlimL:tlimU:10j, 0+toffsetz:np.pi/2+toffsetz/4:10j]
                x=ds*np.cos(u)*np.sin(v) + p[0]
                y=ds*np.sin(u)*np.sin(v) + p[1]
                z=ds*np.cos(v) + p[3]
                self.image.plot_surface( x, y, z,  rstride=2, cstride=2,
                                color='c', alpha=0.3, linewidth=2,
                                edgecolors=np.array((0.5,0.5,0.5,0.5)))

                plt.axis('equal')

        def PlotPathSegment(self, plast, pnext):
                self.image.plot([plast[0],pnext[0]],[plast[1],pnext[1]],[plast[3],pnext[3]],'-k',linewidth=self.path_lw)

        def PlotBall(self, ds):
                pass


        def PlotArrow(self, pos, direction, color, ls='solid',arrowstyle="-|>"):

                v = pos + direction
                a = Arrow3D([pos[0], v[0]], [pos[1], v[1]], 
                                [pos[3], v[3]], mutation_scale=20, 
                                lw=3, arrowstyle=arrowstyle, linestyle=ls, color=color)
                self.image.add_artist(a)
                return a
                
        def expspace(self, tstart, tend, tsamples):
                tlin = np.linspace(tstart, tend, tsamples)
                tpow = (tlin-tstart)**2
                dt = tpow[-1]-tpow[0]
                tscale = (tend-tstart)/dt
                texp = tscale * tpow + tstart
                return texp

        def add_points(self, q):
                hull = ConvexHull(q,qhull_options=self.qhull_options)    
                q = q[hull.vertices,:]
                self.poly.append( q )

                #if self.pts is None:
                #        self.pts = q
                #else:
                #        self.pts = np.vstack((self.pts,q))

        def Plot(self):

                from mpl_toolkits.mplot3d import Axes3D
                from mpl_toolkits.mplot3d.art3d import Poly3DCollection
                from matplotlib.tri import Triangulation, TriAnalyzer

                for i in range(0,len(self.poly)):

                        if i < len(self.poly)-1:
                                X = np.vstack((self.poly[i],self.poly[i+1]))
                        else:
                                X = self.poly[i]

                        X = np.vstack((X[:,0],X[:,1],X[:,3])).T
                        hull = ConvexHull(X,qhull_options=self.qhull_options)    


                        x,y,z=X.T
                        tri = Triangulation(x, y, triangles=hull.simplices)

                        triangle_vertices = np.array([np.array([[x[T[0]], y[T[0]], z[T[0]]],
                                [x[T[1]], y[T[1]], z[T[1]]],
                                [x[T[2]], y[T[2]], z[T[2]]]]) for T in tri.triangles])

                        self.tri = tri

                        tri = Poly3DCollection(triangle_vertices)
                        if i == len(self.poly)-1:
                                tri.set_color(self.rs_last_color)
                                tri.set_edgecolor(self.rs_last_edge_color)
                                self.image.scatter(x,y,z, 'ok', color=np.array((0,0,1.0,0.1)),s=30)

                        else:
                                tri.set_color(self.rs_color)
                                tri.set_edgecolor('None')

                        self.image.add_collection3d(tri)


                self.image.xaxis.get_major_formatter().set_powerlimits((0, 2))
                self.image.yaxis.get_major_formatter().set_powerlimits((0, 2))
                self.image.zaxis.get_major_formatter().set_powerlimits((0, 2))

                self.image.tick_params(axis='both', which='major', pad=15)

                self.image.xaxis.get_offset_text().set_size(self.fs)
                self.image.yaxis.get_offset_text().set_size(self.fs)
                self.image.zaxis.get_offset_text().set_size(self.fs)

                for tick in self.image.xaxis.get_major_ticks():
                        tick.label.set_fontsize(self.fs) 

                for tick in self.image.yaxis.get_major_ticks():
                        tick.label.set_fontsize(self.fs) 

                for tick in self.image.zaxis.get_major_ticks():
                        tick.label.set_fontsize(self.fs) 

                #plt.axis('equal')
                #plt.autoscale(enable=True, axis='y', tight=True)
                self.image.relim()
                self.image.autoscale_view(True,False,True)


        def PlotSave(self, filename=None):

                if filename is not None:
                        self.filename = filename
                        
                self.fig.savefig(self.filename,format='svg', dpi=1200)

                svgname = self.filename
                pdfname = os.path.splitext(self.filename)[0]+'.pdf'

                syscmd = 'mogrify -format pdf -trim '+svgname
                os.system(syscmd)
                syscmd = 'pdfcrop '+pdfname+' '+pdfname
                os.system(syscmd)
                os.system('cp /home/`whoami`/git/openrave/sandbox/WPI/images/*.pdf /home/`whoami`/git/papers/images/simulation/')

        def PlotShow(self):
                plt.show()

