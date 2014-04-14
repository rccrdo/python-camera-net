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
import matplotlib.path as mpath
import matplotlib.patches as mpatches

from math2D import *
from target import Target

class TargetData():
    def __init__(self, target_id, pos, vel, angle, area_id):
        assert isinstance(pos, Point2D)
        assert isinstance(vel, Point2D)
        assert numeric_real(angle)
        assert angle == normalize_angle(angle)
        assert area_id

        self.id = target_id
        self.pos = pos
        self.vel = vel
        self.angle = angle
        self.area_id = area_id        

class Camera():
    def __init__(self, orig, rot, radius, fullfov, camid="", bpoints=[], fovpoints=[]):
        assert isinstance(orig, Point2D)
        assert numeric_real(rot)
        assert numeric_real(radius)
        assert radius > 0.
        assert numeric_real(fullfov)
        assert fullfov > 0. and fullfov < math.pi

        self._targetdata = {}
        self._last_targetdata = {}
        self._dt = None
        self._orig = orig
        self._rot  = normalize_angle(rot)
        self._rotl = normalize_angle(rot - fullfov/2.)
        self._roth = normalize_angle(rot + fullfov/2.)
        self._radius = radius
        self._radius2 = radius**2
        self._fullfov = normalize_angle(fullfov)
        self._id = camid

        self._bangles = []
        self._bpoints = bpoints
        self._compute_bangles(bpoints)
        self._fovpoints = fovpoints
        
        # plot primitives
        self._fullcoverage_patch = self._coverage_patch(rot, fullfov, radius, (0.65,0.65,0.65))
        self._patches = []
        self._lines = []
        self._badge = self._camera_badge()

    #def _local_angle(self, p):
        # compute the angle in local coordinates with respect to the camera
        # depth axis
        #return normalize_angle(self._angle(p) - self._rot)

    def _angle(self, p):
        assert isinstance(p, Point2D)
        # compute the angle in global coordinates with respect to the camera
        # depth axis
        return Line2D(self._orig, p).angle()

    def _append_bangle(self, angle):
        assert angle == normalize_angle(angle)

        # check if the bisecting line actually bisects the coverage area
        if not normalized_angle_in_range(angle, self._rotl, self._roth):
            return

        bangles = self._bangles            
        bangles.append(angle)

        # reorder the list of bangles to make detection easier later on
        if self._rotl >=0 and self._roth <= 0:
            # in this case we must sort positive and negative angles separately
            posi = []
            nega = []
            for a in bangles:
                if a >= 0:
                    posi.append(a)
                else:
                    nega.append(a)
            bangles = sorted(posi) + sorted(nega)
        else:
            bangles = sorted(bangles)

        # purge bisecting lines with too similar angles
        if 0:
            done = False
            while not done:
                done = True
                for i in xrange(0,len(bangles)-1):
                    #print len(bangles), i
                    anglel = bangles[i]
                    angleh = bangles[i+1]
                    # pop bangles that are less than 2 degrees apart
                    union_rangeh = normalize_angle(anglel+math.pi/180.*2.)
                    #print self._id, anglel, union_rangeh, angleh, normalized_angle_in_range(angleh, anglel, union_rangeh)
                    if normalized_angle_in_range(angleh, anglel, union_rangeh):
                        #print "popping bangle at index", i+1
                        bangles.pop(i+1)
                        done = False
                        break

        self._bangles = bangles
        
    def _compute_bangles(self, bpoints):
        if not bpoints:
            return

        for p in bpoints:
            angle = self._angle(p)
            self._append_bangle(angle)
       
    def _camera_badge(self):
        # shortcut into matplotlib namespace
        Path = mpath.Path

        badge_size = 0.3
        points = numpy.array([
             (-badge_size/2., -badge_size/2.),
             ( badge_size/2., -badge_size/2.),
             ( badge_size/2.,  badge_size/2.),
             (-badge_size/2.,  badge_size/2.),
             (-badge_size/2., -badge_size/2.),
        ]) 

        # rotate the badge points and tranlate them to the camera origin
        verts = transform_points(points, self._rot, self._orig.array())
        codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
        #return mpatches.PathPatch(mpath.Path(verts, codes), alpha=1, facecolor=(0.5,0.5,0.75))
        #print self._orig.tuple()
        return matplotlib.patches.Circle(self._orig.tuple(), radius=.2, alpha=1, facecolor=(1,1,1))

    def _compute_blines(self):
        # create the bisceting lines
        self._lines = []
        clr_blines = (1,0,0)
        for angle in self._bangles:
            p1 = self._orig
            p2 = Point2D(self._radius, 0)
            p2 = rotate_point(p2, angle)
            p2.x += p1.x
            p2.y += p1.y
            if 0:
                self._lines.append(matplotlib.lines.Line2D([p1.x, p2.x], [p1.y, p2.y], color=clr_blines, linestyle=':'))

    def _update_patches(self):
        # namespace shortcut for the codes below
        Path = mpath.Path

        # build the camera badge patch
        badge_size = 0.3
        points = numpy.array([
             (-badge_size/2., -badge_size/2.),
             ( badge_size/2., -badge_size/2.),
             ( badge_size/2.,  badge_size/2.),
             (-badge_size/2.,  badge_size/2.),
             (-badge_size/2., -badge_size/2.),
        ]) 

        # rotate the badge points and tranlate them to the camera origin
        verts = transform_points(points, self._rot, self._orig.array())
        codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
        self._badge_patch = mpatches.PathPatch(mpath.Path(verts, codes), alpha=1, facecolor=(0.5,0.5,0.75))
        self._badge_patch = matplotlib.patches.Circle([self._posx, self._posy], radius=.2, alpha=1, facecolor=(1,0.75,0.75))

        # build the camera coverage patch
        r = self._radius
        w = self._fov
        nr_segments = 20
        points = [(0,0)]
        for i in numpy.arange(0, nr_segments+1):
            theta = -w/2. + i*w/float(nr_segments)
            points.append((math.cos(theta), math.sin(theta)))
        points.append((0,0))
        points = numpy.array(points)

        # rotate the cone patch points and tranlate them to the camera origin
        verts = transform_points(r*points, self._rot, self._orig.array())
        codes = [Path.MOVETO] + [Path.LINETO]*(nr_segments+1) + [Path.CLOSEPOLY]
        self._coverage_patch = mpatches.PathPatch(mpath.Path(verts, codes), alpha=0.10, facecolor=(0.25,0.25,0.25), linewidth=0.5)

        verts = transform_points(0.5*points, self._rot, self._orig.array())
        codes = [Path.MOVETO] + [Path.LINETO]*(nr_segments+1) + [Path.CLOSEPOLY]
        self._target_patch = mpatches.PathPatch(mpath.Path(verts, codes), alpha=0.45, facecolor=(1,0.2,0.2), linewidth=0.5)


        # create the bisceting lines
        self._plot_lines = []
        clr_blines = (1,0,0)
        for angle in self._bangles:
            p1 = self._orig
            p2 = Point2D(self._radius, 0)
            p2 = rotate_point(p2, angle)
            p2.x += p1.x
            p2.y += p1.y
            if 0:
                self._plot_lines.append(matplotlib.lines.Line2D([p1.x, p2.x], [p1.y, p2.y], color=clr_blines, linestyle=':'))

    def _coverage_patch(self, rot, arc, radius, color, alpha=0.1):
        assert numeric_real(alpha)
        assert alpha >0 and alpha <= 1
        # namespace shortcut for the codes below
        Path = mpath.Path

        # build the camera coverage patch
        nr_segments = 20
        points = [(0,0)]
        for i in numpy.arange(0, nr_segments+1):
            theta = -arc/2. + i*arc/float(nr_segments)
            points.append((math.cos(theta), math.sin(theta)))
        points.append((0,0))
        points = numpy.array(points)

        # rotate the cone patch points and tranlate them to the camera origin
        verts = transform_points(radius*points, rot, self._orig.array())
        codes = [Path.MOVETO] + [Path.LINETO]*(nr_segments+1) + [Path.CLOSEPOLY]
        return mpatches.PathPatch(mpath.Path(verts, codes), alpha=alpha, facecolor=color, linewidth=0.5)

    def origin(self):
        return self._orig

    def _angle_to_areaid(self, angle):
        assert numeric_real(angle)
        assert angle == normalize_angle(angle)

        # determine the exact coverage sub-area id
        area_id = self._id
        if len(self._bangles):
            anglel = self._rotl
            area_idx = 0
            i = 0
            for i in xrange(0,len(self._bangles)):
                angleh = self._bangles[i]
                if normalized_angle_in_range(angle, anglel, angleh):
                    area_idx = i + 1
                    break
                anglel = angleh
            if not area_idx:
                # if we haven't found it yet then it must be in the last area
                area_idx = len(self._bangles) + 1
            
            # use letters if possible, use numbers otherwise
            if area_idx < 26:
                area_id += chr(ord('a') + area_idx -1)
            else:
                area_id += '[%d]' % area_idx

        return area_id


    def detection_data(self):
        return self._targetdata.values()

    def plot(self, axis):
        assert self._badge
        assert self._fullcoverage_patch
        #assert self._target_patch
        
        if 1:#self._fovpoints:
        	for p in self._fovpoints:
	            axis.add_line(matplotlib.lines.Line2D([self._orig.x, p.x],[self._orig.y, p.y], color=(0.25,0.25,0.25,0.5), zorder=-100))
        	for p in self._bpoints:
        	    pass
	            #axis.add_line(matplotlib.lines.Line2D([self._orig.x, p.x],[self._orig.y, p.y], color=(1.,0.25,0.25, 1), zorder=-100, linestyle=':', linewidth=2))
        	
        else:
            pass
        if 1:
	        #axis.add_patch(self._fullcoverage_patch)

	        #for patch in self._patches:    
	        #    axis.add_patch(patch)

	        for patch in self._custom_patches():
        	    axis.add_patch(patch)
        
	        #if len(self._detected):    
	        #    axis.add_patch(self._target_patch)

	        #for line in self._lines:
        	#    axis.add_line(line)
        
	        #for data in self._detected:
	        #    p = data.pos
	            #print data.area_id
	        #    axis.add_line(matplotlib.lines.Line2D([self._orig.x, p.x],[self._orig.y, p.y], color=(1,0.,0.), alpha=1, zorder=-100))

        axis.add_patch(self._badge)

        if self._id:
            axis.text(self._orig.x, self._orig.y, self._id,
               verticalalignment='center', horizontalalignment='center', family='sans-serif',
               color='black', fontsize=15)
        return
        
    def _track_target(self, targetdata):
        assert isinstance(targetdata, TargetData)
            
        _id = targetdata.id
        try:
            oldpos = self._last_targetdata[_id].pos
            vx = (targetdata.x-oldpos.x)/self._dt
            vy = (targetdata.y-oldpos.y)/self._dt
            targetdata.vel.x = vx
            targetdata.vel.y = vy
        except:
            pass
            
        assert _id not in self._targetdata
        self._targetdata[_id] = targetdata


class CameraFixed(Camera):
    def __init__(self, orig, rot, radius, fov, camid="", bpoints=[], fovpoints=[]):
        Camera.__init__(self, orig, rot, radius, fov, camid, bpoints, fovpoints)

        self._compute_blines()
        self._dt = None

    def _custom_patches(self):
        return []

    def step(self, time, dt):
        self._dt = dt
        self._last_targetdata = self._targetdata
        self._targetdata = {}
        return
        
    def detect(self, target, walls):
        assert isinstance(target, Target)
        assert walls
        
        pos = target.pos()
        line = Line2D(self._orig, pos)
        
        # random detection artifacts
        # pre-refactoring <- this code must be checked for good
        if 0:#numpy.random.randint(50) == 0:
            area_id = self._id
            if len(self._bangles):
                area_idx = numpy.random.randint(1,len(self._bangles)+1)
    
                # use letters for up to 9 bisecting lines, use numbers otherwise
                if area_idx < 10:
                    area_id += chr(ord('a') + area_idx -1)
                else:
                    area_id += '.%d' % area_idx
            
            #print 'area_id', area_id
            self._detected.append(TargetData(target.id(), pos, Point2D(0,0), 0, area_id))
            return

        if line.norm2() > self._radius2:
            # the target is more distant than the radius of our coverage cone
            return

        # check if the target is within the coverage area without considering
        # occluding objects: simply check its angular position with respect to
        # our origin
        angle = self._angle(pos)
        if not normalized_angle_in_range(angle, self._rotl, self._roth):
            return
        
        # check if the line of sight intersects any occluding object
        for wall in walls:
            if line.intersects(wall):
                return
        
        area_id = self._angle_to_areaid(angle)
        self._track_target(TargetData(target.id(), pos, Point2D(0,0), angle, area_id))

_PTZ_STATE_PATROLLING = 0
_PTZ_STATE_TRACKING = 1
_PTZ_DIRECTION_RIGHT = 0
_PTZ_DIRECTION_LEFT = 1
class CameraPTZ(Camera):
    def __init__(self, orig, rot, radius, fullfov, fov, velocity, camid="", bpoints=[], fovpoints=[]):
        Camera.__init__(self, orig, rot, radius, fullfov, camid, bpoints, fovpoints)

        assert numeric_real(fov)
        assert fov > 0 and fov < math.pi
        assert fov <= self._fullfov
        assert numeric_real(velocity)
        assert velocity > 0.
        
        self._local_angle = fov/2
        self._fov = fov
        self._vel = velocity
        self._state = _PTZ_STATE_PATROLLING
        self._rotdir = _PTZ_DIRECTION_LEFT

	# split the extended fov in sectors
        #angle = normalize_angle(self._rotl + fov)
        #while True:
        #    self._append_bangle(angle)
        #    angle = normalize_angle(angle + fov)
        #    if not normalized_angle_in_range(angle, self._rotl, normalize_angle(self._roth - fov/2.)):
        #        break
                
        self._compute_blines()

    def _custom_patches(self):
        color = (1.,1.,0.75)
        if len(self._last_targetdata):
            color = (1.,0.1,0.1)
        return [self._coverage_patch(self._rotl + self._local_angle, self._fov, 0.9, color, 0.75)]

    def step(self, time, dt):
        self._dt = dt
        self._last_targetdata = self._targetdata
        self._targetdata = {}

        if len(self._last_targetdata):
            self._state = _PTZ_STATE_TRACKING
        else:
            self._state = _PTZ_STATE_PATROLLING

        assert self._state in (_PTZ_STATE_PATROLLING, _PTZ_STATE_TRACKING)
        if self._state == _PTZ_STATE_PATROLLING:
            #print "  patrolling"
            assert self._rotdir in (_PTZ_DIRECTION_LEFT, _PTZ_DIRECTION_RIGHT)
            if self._rotdir == _PTZ_DIRECTION_LEFT:
                #print "    patrolling left"
                new_angle = self._local_angle + dt*self._vel
                limit = self._fullfov -self._fov/2
                da = limit - new_angle
                #print limit, new_angle, da
                if da < 0:
                    #print "    patrolling changing direction"
                    new_angle = limit + da
                    self._rotdir = _PTZ_DIRECTION_RIGHT
                self._local_angle = new_angle
            else:
                #print "    patrolling right"
                assert self._rotdir == _PTZ_DIRECTION_RIGHT
                new_angle = self._local_angle - dt*self._vel
                limit = 0 + self._fov/2
                da = limit - new_angle
                #print limit, new_angle, da
                if da > 0:
                    #print "    patrolling changing direction"
                    new_angle = limit + da
                    self._rotdir = _PTZ_DIRECTION_LEFT
                self._local_angle = new_angle
        else:
            assert self._state == _PTZ_STATE_TRACKING
            #print "  tracking camid", self._id
            # currently this works only for 1 target at a time
            #print len(self._
            assert len(self._last_targetdata) == 1
            target = self._last_targetdata['target-hmm']
            pred = Point2D(target.pos.x + target.vel.x*dt, target.pos.y + target.vel.y*dt)
            #print "  pred", target.pos, pred
            angle = self._angle(pred)
            if not normalized_angle_in_range(angle, self._rotl, self._roth):
                # target exited our fullfov, go back to patrolling
                #print "  tracking -> back to patrolling"
                self._state = _PTZ_STATE_PATROLLING
                self.step(time, dt)

            cur_tracking_angle = normalize_angle(self._local_angle + self._rotl)
            #print "angle, cur_tracking_angle ", angle, cur_tracking_angle 

            if cur_tracking_angle < -math.pi/2. and angle > math.pi/2:
                #print "  case 1"
                da = -((math.pi-angle) + (cur_tracking_angle + math.pi))
            elif cur_tracking_angle > math.pi/2. and angle < -math.pi/2:
                #print "  case 2"
                da = (math.pi-cur_tracking_angle) + (angle + math.pi)
            else:
                da = angle - cur_tracking_angle
                #print "  case 3, da ", da

            # clamp to camera velocity
            avel = da / dt
            #print " avel", avel
            avel = max(-self._vel, min(self._vel, avel))
            #print " avel clamped", avel
            da = avel * dt
            #print " da clamped to speed", da
            loc_angle = self._local_angle + da
            #print " self._local_angle, new unclamped loc angle", self._local_angle, loc_angle
            
            limitl = self._fov/2.
            limith = self._fullfov - self._fov/2.
            if loc_angle < limitl:
                loc_angle = limitl
            elif loc_angle > limith:
                loc_angle = limith
                
            #print self._local_angle, loc_angle, limitl, limith, angle
            self._local_angle = loc_angle
            
    def detect(self, target, walls):
        assert isinstance(target, Target)
        assert walls
        
        pos = target.pos()
        line = Line2D(self._orig, pos)
        
        #if line.norm2() > self._radius2:
            # the target is more distant than the radius of our coverage cone
        #    return

        # check if the target is within the current coverage area without
        # considering occluding objects: simply check its angular position
        # with respect to our origin
        angle = self._angle(pos)
        fovl = normalize_angle(self._rotl + self._local_angle - self._fov/2)
        fovh = normalize_angle(self._rotl + self._local_angle + self._fov/2)
        #print "camid fovl, fovh, rotl, roth", self._id, fovl, fovh-0.0001, self._rotl, self._roth
        assert normalized_angle_in_range(fovl, self._rotl, self._roth)
        assert normalized_angle_in_range(fovh, self._rotl, self._roth, epsilon= 0.000001)
        if not normalized_angle_in_range(angle, fovl, fovh):
            return

        # check if the line of sight intersects any occluding object
        for wall in walls:
            if line.intersects(wall):
                return

        area_id = self._angle_to_areaid(angle)
        #print "detected in ", area_id
        self._track_target(TargetData(target.id(), pos, Point2D(0,0), angle, area_id))




