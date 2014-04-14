# Copyright (c) 2013 Riccardo Lucchese, riccardo.lucchese at gmail.com
#
# This software is provided 'as-is', without any express or implied
# warranty. In no event will the authors be held liable for any damages
# arising from the use of this software.
#
# Permission is granted to anyone to use this software for any purpose,
# including commercial applications, and to alter it and redistribute it
# freely, subject to the following restrictions:
#
#    1. The origin of this software must not be misrepresented; you must not
#    claim that you wrote the original software. If you use this software
#    in a product, an acknowledgment in the product documentation would be
#    appreciated but is not required.
#
#    2. Altered source versions must be plainly marked as such, and must not be
#    misrepresented as being the original software.
#
#    3. This notice may not be removed or altered from any source
#    distribution.

import math
import numpy

import matplotlib
from matplotlib import rc
rc('font', family='sans-serif')
import matplotlib.pyplot as plt

from math2D import numeric_real, Point2D, Line2D, Rectangle
from geometry import Geometry
from camera import Camera
from target import Target


class Environment():
    def __init__(self, width, height):
        assert numeric_real(width)
        assert width > 0
        assert numeric_real(height)
        assert height > 0

        self._initialized = False
        self._cameras = []
        self._targets = []
        self._geometries = []
        self._walls = []
        self._cur_time = 0
        self._observation = {}

        self._rect = Rectangle(Point2D(0, 0), 0., width, height)
        if not self._rect.area():
            print("warning, environment has zero area")
            return

        self._figure = None

        self._initialized = True

    def _contains_point(self, p):
        return self._rect.contains_point(p)
        
    def _contains_line(self, line):
        return self._rect.contains_line(line)

    def _update_walls(self):
        self._walls = []
        for geom in self._geometries:
            for line in geom.border():
                if line in self._walls:
                    continue
                if not self._contains_line(line):
                    print "warning, discarding wall with both points outside the environment ", line
                    continue

                self._walls.append(line)

    def add_geometry(self, geom):
        assert isinstance(geom, Geometry)
        assert geom not in self._geometries
        self._geometries.append(geom)
        self._update_walls()
    
    def add_camera(self, camera):
        assert isinstance(camera, Camera)
        assert camera not in self._cameras
        o = camera.origin()
        if not self._contains_point(o):
            print "warning, camera with id \"%\" has origin outside the environment ", o
        self._cameras.append(camera)

    def add_cameras(self, cameras):
        for cam in cameras:
            self.add_camera(cam)

    def add_target(self, target):
        assert isinstance(target, Target)
        assert target not in self._targets
        self._targets.append(target)

    def add_targets(self, targets):
        for target in targets:
            self.add_target(target)

    def time(self):
        return self._cur_time

    def step(self, dt):
        assert numeric_real(dt)
        assert dt > 0

        cur_time = self._cur_time
        
        # update target positions
        for target in self._targets:
            target.step(cur_time, dt, self._walls)

        # detect targets in each camera coverage area
        data = []
        for cam in self._cameras:
            # step ptz cameras          
            cam.step(cur_time, dt)

            for target in self._targets:
                cam.detect(target, self._walls)
                
            data.extend(cam.detection_data())
        
        observ = {}
        for entry in data:
             #print entry
             _id = entry.id
             if _id in observ.keys():
                 observ[_id].append(entry.area_id)
             else:
                 observ[_id] = [entry.area_id]
        #for key in sorted(observ.keys()):
        #    print "  observ for target \"%s\": nr. %d, %s" % (key, len(observ[key]), observ[key])

        self._observation = {}
        for targetid,val in observ.iteritems():
            if val == []:
                val = '0'
            else:
                val = '-'.join(sorted(val))
            self._observation[targetid] = val


        # open the building passage    
        #if self._cur_time == 30:
        #    self._geometries[0].set_passage_open(1)
        #    self._update_walls()

        self._cur_time += dt


    def get_observation(self):
        return self._observation
        
    def get_observation_for_target(self, targetid):
	assert targetid

        try:
            entry = self._observation[targetid]
        except:
            entry = '0'        
        return entry

    def plot(self, savepath=None):
        if not self._figure:
            #self._figure = plt.figure(figsize=(6.5,4.5))
            self._figure = plt.figure(figsize=(9,6))
            plt.ion()
            plt.show()

        fig = self._figure
        plt.figure(fig.number)
        # setup figure for this iteration
        fig.clf()
        fig.patch.set_facecolor((1,1,1))

        # setup axis, limits, grid
	gs = matplotlib.gridspec.GridSpec(1, 1, left=0.05, right=0.95, top=0.95, bottom=0.05, wspace=0, hspace=0)
        axis = plt.subplot(gs[0,0], aspect='equal')
        plt.axis('off')
        if 0:
            plt.axis('on')
        #plt.grid('on')
        xmin, xmax = self._rect.xlim()
        plt.xlim(xmin -0.25, xmax)# +0.5)
        ymin, ymax = self._rect.ylim()
        plt.ylim(ymin -0.25, ymax +0.25)

        # plot cameras
        for cam in self._cameras:
            cam.plot(axis)

        # plot occlusion geometries
        for geom in self._geometries:
            geom.plot(axis)
        
        # plot targets
        for target in self._targets:
            target.plot(axis)

        # plot some stats
        #axis.text(0.1, 5, str(self._cur_time), verticalalignment='center',
        #	horizontalalignment='center', family='sans-serif',
        #	color='black', fontsize=15)

        fig.canvas.draw()
        
        if savepath:
            plt.savefig(savepath)
