from shapely.ops import cascaded_union, polygonize
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
        rs_last_color = np.array((0.5,0.5,1.0,0.2))

        rs_alpha = 0.8
        rs_last_alpha = 0.9

        path_lw = 3

        force_color = np.array((1.0,0,0))
        tangent_color = np.array((0,0,0))
        orientation_color = np.array((0.9,0.0,0.5))
        fs = 28
        fs_label = 34
        fs_title = 46
        image = None

        def __init__(self, p, s, dp, force, R, amin, amax):
                
                self.pts = None
                self.poly = []
                self.fig = plt.figure(facecolor='white')
                from mpl_toolkits.mplot3d import Axes3D
                self.image = self.fig.gca(projection='3d')

                
                M = 10000
                tstart = 0.001
                tend = 0.01
                tsamples= 15

                Ndim = p.shape[0]
                poly = []
                
                for dt in self.expspace(tstart,tend,tsamples):
                        print dt
                        dt2 = dt*dt*0.5
                        q = np.zeros((Ndim,M))
                        for i in range(0,M):
                                ## sample random control
                                ar = np.random.uniform(amin, amax)
                                control = np.dot(R,ar)
                                q[:,i] = p + dt*s*dp + dt2 * force + dt2 * control

                        self.add_points( q.T )

                tstring = 'Reachable Set (<T='+str(dt)+')'
                self.filename = 'images/reachableset_3dcar_ori'+str(np.around(p[3],decimals=2))
                self.filename = re.sub('[.]', '-', self.filename)
                self.filename += '.png'

                self.image.set_title(tstring, fontsize=self.fs_title, y=1.1)
                self.image.set_xlabel('\n\nX-Position [m]', fontsize=self.fs_label)
                self.image.set_ylabel('\n\nY-Position [m]', fontsize=self.fs_label)
                self.image.set_zlabel('\n\n$\\theta$-Position [rad]', fontsize=self.fs_label)
                #self.image.xaxis.labelpad = 100

                self.arrow_head_size = np.linalg.norm(dt2*force)/5
                self.hw = 0.5*self.arrow_head_size
                self.lw = 0.2*self.arrow_head_size

                pnext = p+dt*s*dp+dt2*force

                self.image.scatter(p[0], p[1], p[3], 'ok', s=30)
                self.image.scatter(pnext[0],pnext[1],pnext[3],  'ok', s=10)

                dori = np.zeros((Ndim))
                dori[0:2] = np.dot(Rz(p[3]),ex)[0:2]
                dori[0:2] = 0.5*(dori[0:2]/np.linalg.norm(dori[0:2]))
                dori[3] = p[3]

                if s < 0.001:
                        dnf = np.linalg.norm(force)
                        nforce = force/dnf
                        si = np.linalg.norm(dt2*dnf/(dt*np.linalg.norm(dp)))

                        ## visualize path
                        pathnext = p + 1.5*dt*si*dp
                        pathlast = p - 0.2*dt*si*dp

                        self.PlotPathSegment( pathlast, pathnext)

                        arrow0 = self.PlotArrow(p, dt*si*dori, self.orientation_color)
                        arrow0 = self.PlotArrow(pnext, dt*si*dori, self.orientation_color)
                        arrow2 = self.PlotArrow(p, dt2*force, self.force_color)

                        plt.legend([arrow0,arrow2,],
                                        ['Orientation','Force',],
                                        fontsize=self.fs,
                                        loc=self.loc)
                else:
                        pathnext = p + 1.5*dt*s*dp
                        pathlast = p - 0.2*dt*s*dp

                        self.PlotPathSegment( pathlast, pathnext)

                        arrow0 = self.PlotArrow(p, dt*s*dori, self.orientation_color)
                        arrow1 = self.PlotArrow(p, dt*s*dp, self.tangent_color)
                        arrow2 = self.PlotArrow(dt*s*dp, dt2*force, self.force_color)

                        arrow0 = self.PlotArrow(pnext, dt*s*dori, self.orientation_color)
                        arrow2 = self.PlotArrow(p, dt2*force, self.force_color)

                        plt.legend([arrow0,arrow1,arrow2,],
                                        ['Orientation','Velocity/Tangent Path','Force',],
                                        fontsize=self.fs,
                                        loc=self.loc)
                plt.tight_layout()
                plt.axis('equal')

        def PlotPathSegment(self, plast, pnext):
                self.image.plot([plast[0],pnext[0]],[plast[1],pnext[1]],[plast[3],pnext[3]],'-k',linewidth=self.path_lw)


        def PlotArrow(self, pos, direction, color):

                v = pos + direction
                a = Arrow3D([pos[0], v[0]], [pos[1], v[1]], 
                                [pos[3], v[3]], mutation_scale=20, 
                                lw=3, arrowstyle="-|>", color=color)
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

                for i in range(0,len(self.poly)):

                        if i < len(self.poly)-1:
                                X = np.vstack((self.poly[i],self.poly[i+1]))
                                print i,i+1
                        else:
                                X = self.poly[i]
                                print i

                        X = np.vstack((X[:,0],X[:,1],X[:,3])).T
                        hull = ConvexHull(X,qhull_options=self.qhull_options)    

                        from matplotlib.tri import Triangulation

                        x,y,z=X.T
                        tri = Triangulation(x, y, triangles=hull.simplices)
                        triangle_vertices = np.array([np.array([[x[T[0]], y[T[0]], z[T[0]]],
                                [x[T[1]], y[T[1]], z[T[1]]],
                                [x[T[2]], y[T[2]], z[T[2]]]]) for T in tri.triangles])

                        tri = Poly3DCollection(triangle_vertices)
                        if i == len(self.poly)-1:
                                tri.set_color(self.rs_last_color)
                                tri.set_edgecolor('k')
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

                self.fig.set_size_inches(18.5, 18.5, forward=True)
                plt.savefig(self.filename,format='svg', dpi=1200)


                svgname = self.filename
                pdfname = os.path.splitext(self.filename)[0]+'.pdf'

                syscmd = 'mogrify -format pdf -trim '+svgname
                os.system(syscmd)
                syscmd = 'pdfcrop '+pdfname+' '+pdfname
                os.system(syscmd)
                os.system('cp /home/`whoami`/git/openrave/sandbox/WPI/images/*.pdf /home/`whoami`/git/papers/images/simulation/')

                plt.show()

