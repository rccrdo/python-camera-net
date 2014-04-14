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

import matplotlib
import matplotlib.path as mpath
import matplotlib.patches as mpatches

from math2D import *


class Geometry():
    """
    Base clase used by occlusion geometries in the environment.
    
    Objects are modeled with piece-wise linear boundaries.
    """
    def __init__(self):
        self._plot_lines = []
        self._plot_patches = []
        self._border = []

    def plot(self, axis):
        """
        Plot pathecs and lines for this geometry on the given matplotlib axis
        """
        for patch in self._plot_patches:
            axis.add_patch(patch)

        for line in self._plot_lines:
            axis.add_line(line)

    def border(self):
        """
        Return the list of 2D lines that form the geometry boundary
        """
        return self._border

    def _build_lines(self, points):
        """
        Helper func used to build an ordered list of 2D lines connenecting
        the 2D points in the numpy vector 'points'
        """
        assert isinstance(points, numpy.ndarray)
        assert points.shape[1] == 2

        i=0
        lines = []
        while i < points.shape[0]-1:
            p1 = points[i,:]
            p2 = points[i+1,:]
            line = Line2D(Point2D(p1[0], p1[1]), Point2D(p2[0], p2[1]))
            if not line.norm2():
                print "warning, skipping zero-length line on geometry border at " + p1
            else:
                lines.append(line)
            i += 1
        return lines


class Building(Geometry):
    """
    Geometry class for the building complex
    """
    def __init__(self):
        Geometry.__init__(self) 
        self._passage_open = False
        self._left_building_patch = None
        self._right_building_patch = None
        self._passage_closed_patch = None    
        self._passage_open_patch = None
        self._left_building_plot_lines = None
        self._right_building_plot_lines = None
        self._passage_closed_plot_lines = None
        self._passage_internal_plot_lines = None
        self._left_building_walls = None
        self._right_building_walls = None
        self._passage_closed_walls = None
        self._passage_open_walls = None
        self.transform(2, 1.5, 0)

    def _update_plot_lists(self):
        patches = [self._left_building_patch, self._right_building_patch]
        lines = self._left_building_plot_lines		\
                + self._right_building_plot_lines	\
                + self._passage_internal_plot_lines

        walls = self._left_building_walls + self._right_building_walls
       
        if self._passage_open:
            patches.append(self._passage_open_patch)
            walls += self._passage_open_walls
        else:
            patches.append(self._passage_closed_patch)
            lines.extend(self._passage_closed_plot_lines)
            walls += self._passage_closed_walls

        # defensinve checks
        assert None not in patches
        assert None not in lines
        assert None not in walls
        self._plot_patches = patches
        self._plot_lines = lines
        self._border = walls

    def set_passage_open(self, isopen):
        """
        update the state of the passage
        """
        self._passage_open = bool(isopen)
        self._update_plot_lists()

    def transform(self, posx, posy, rot):
        """
        Rotate and translate the building position in 2D space
        """
        self._posx = posx
        self._posy = posy
        self._rot = rot

        # shortcut into matplotlib's namespace used below
        Path = mpath.Path

        points_left_building = numpy.array([
             (0, 0),
             (1, 0),
             (1, 1),
             (3.5, 1),
             (3.5, 2),
             (2, 2),
             (2, 3),
             (0, 3),
             (0, 0),
        ])
        points_right_building = numpy.array([
             (4, 1),
             (4.5, 1),
             (4.5, 0),
             (5.5, 0),
             (5.5, 1),
             (5, 1),
             (5, 3.5),
             (4.5, 3.5),
             (4.5, 2),
             (4, 2),
             (4, 1),
        ]) 
        points_passage = numpy.array([
             (3.5, 1),
             (4, 1),
             (4, 2),
             (3.5, 2),
             (3.5, 1),             
        ])

        building_color = (0.5,0.5,0.5) 
        building_border_color = (0.25,0.25,0.25) 
        passage_color_closed = (0.3,0.3,0.3) 
        passage_color_open = (0.1,1.,0.1) 
    
        # rotate and tranlate the left building vertices
        verts = transform_points(points_left_building, self._rot, Point2D(posx, posy).array())
        codes = [Path.MOVETO] + [Path.LINETO]*(len(verts)-2) + [Path.CLOSEPOLY]
        self._left_building_patch = mpatches.PathPatch(mpath.Path(verts, codes), alpha=0.20, facecolor=building_color, linewidth=0)
    
        # compute the left building border excluding the wall in
        # common with the passage
        left_building_walls = self._build_lines(verts)
        left_building_walls.pop(3)
        self._left_building_walls = left_building_walls
        self._left_building_plot_lines = []
        for line in left_building_walls:
            x1, y1 = line.p1.tuple()
            x2, y2 = line.p2.tuple()
            mplib_line = matplotlib.lines.Line2D([x1, x2], [y1, y2], color=building_border_color)
            self._left_building_plot_lines.append(mplib_line)



        # rotate and tranlate the right building vertices
        verts = transform_points(points_right_building, self._rot, Point2D(posx, posy).array())
        codes = [Path.MOVETO] + [Path.LINETO]*(len(verts)-2) + [Path.CLOSEPOLY]
        self._right_building_patch = mpatches.PathPatch(mpath.Path(verts, codes), alpha=0.20, facecolor=building_color, linewidth=0)

        # compute the right building border excluding the wall in
        # common with the passage
        right_building_walls = self._build_lines(verts)
        right_building_walls.pop(9)
        self._right_building_walls = right_building_walls
        self._right_building_plot_lines = []
        for line in right_building_walls:
            x1, y1 = line.p1.tuple()
            x2, y2 = line.p2.tuple()
            mplib_line = matplotlib.lines.Line2D([x1, x2],[y1, y2], color=building_border_color)
            self._right_building_plot_lines.append(mplib_line)
    
    

        # rotate and tranlate the central passage vertices
        verts = transform_points(points_passage, self._rot, Point2D(posx, posy).array())
        codes = [Path.MOVETO] + [Path.LINETO]*(len(verts)-2) + [Path.CLOSEPOLY]
        self._passage_closed_patch = mpatches.PathPatch(mpath.Path(verts, codes), alpha=0.20, facecolor=passage_color_closed, linewidth=0)
        self._passage_open_patch = mpatches.PathPatch(mpath.Path(verts, codes), alpha=0.20, facecolor=passage_color_open, linewidth=0)
    
        # compute the passage walls
        passage_walls = self._build_lines(verts)
        passage_internal_walls = [passage_walls.pop(1)]
        passage_internal_walls.append(passage_walls.pop(2))
        self._passage_closed_walls = passage_walls
        self._passage_open_walls = passage_internal_walls

        self._passage_closed_plot_lines = []
        for line in passage_walls:
            x1, y1 = line.p1.tuple()
            x2, y2 = line.p2.tuple()
            mplib_line = matplotlib.lines.Line2D([x1, x2],[y1, y2], color=building_border_color)
            self._passage_closed_plot_lines.append(mplib_line)

        self._passage_internal_plot_lines = []
        for line in passage_internal_walls:
            x1, y1 = line.p1.tuple()
            x2, y2 = line.p2.tuple()
            mplib_line = matplotlib.lines.Line2D([x1, x2],[y1, y2], color=building_border_color, linestyle='--')
            self._passage_internal_plot_lines.append(mplib_line)

        self._update_plot_lists()

    
        


