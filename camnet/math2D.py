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


def numeric_real(x):
    """
    Check if x is a scalar numeric type
    """
    return isinstance(x, (int, long, float))

def normalize_angle(phi):
    """
    Normalize the angle phi [rad] in the range (-pi,pi]
    """
    assert numeric_real(phi)
    while phi > math.pi or phi <= -math.pi:
        if phi > math.pi:
            phi -= 2*math.pi
        elif phi <= -math.pi:
            phi += 2*math.pi
    return phi

def normalized_angle_in_range(phi, low, high, epsilon = 0.):
    """
    Return True if the the normalized angle phi [rad]
    is in the range [low,high] \subseteq (-pi,pi].
    """
    assert phi == normalize_angle(phi)
    assert low == normalize_angle(low)
    assert high == normalize_angle(high)

    epsilon = 0.0000001
    if low < high:
        return phi>=(low - epsilon) and phi<= (high + epsilon)
    else:
        # the range contains \pi
        return (phi>=0 and phi >= low) or (phi<=0 and phi <= high)


def rotation_matrix(rot):
    """
    Return a numpy 2D rotation matrix with angle rot [rad].
    """
    assert numeric_real(rot)
    cosphi = math.cos(rot)
    sinphi = math.sin(rot)

    # anti-clockwise rotation for phi > 0
    return numpy.array([cosphi, -sinphi, sinphi, cosphi]).reshape((2,2))

def rotate_point(p, rot):
    """
    Rotate the Point2D p around the origin by an angle rot [rad]

    Returns a Point2D.
    """
    assert isinstance(p, Point2D)
    assert numeric_real(rot)
    ret = p.array().dot(rotation_matrix(rot).transpose())
    return Point2D(ret[0][0], ret[0][1])


def rotate_points(points, rot):
    """
    Rotate the numpy vector of 2D points in 'points' around the
    origin by an angle rot [rad].

    Returns a numpy vector of 2D points.
    """
    assert isinstance(points, numpy.ndarray)
    assert points.shape == (1,2)
    assert numeric_real(rot)
    return points.dot(rotation_matrix(rot).transpose())


def transform_points(points, rot, dxy):
    """
    Transform the numpy vector of 2D points in 'points' around the
    origin by an angle rot [rad] and then translate them by 'dxy'

    Returns a numpy vector of 2D points.
    """
    assert isinstance(points, numpy.ndarray)
    assert points.shape[1] == 2
    assert numeric_real(rot)
    assert isinstance(dxy, numpy.ndarray)
    assert dxy.shape == (1,2)
    return points.dot(rotation_matrix(rot).transpose()) + dxy


class Point2D(object):
    """
    A class representing a 2D point on a plane
    """
    def __init__(self, x, y):
        assert numeric_real(x)
        assert numeric_real(y)
        
        self.set_x(x)
        self.set_y(y)

    def get_x(self):
        return self._x

    def set_x(self, val):
        assert numeric_real(val)
        self._x = val

    x = property(get_x, set_x)

    def get_y(self):
        return self._y

    def set_y(self, val):
        assert numeric_real(val)
        self._y = val
       
    y = property(get_y, set_y)

    def tuple(self):
        return (self._x, self._y)
        
    def array(self):
        return numpy.array((self._x, self._y)).reshape(1,2)
        
    def __repr__(self):
        return "(%.2f,%.2f)" % (self._x, self._y)


class Line2D(object):
    """
    A class representing a 2D segment on a plane
    """
    def __init__(self, p1, p2):
        assert isinstance(p1, Point2D)
        assert isinstance(p2, Point2D)

        self.set_p1(p1)
        self.set_p2(p2)

    def get_x1(self):
        return self._x1

    def set_x1(self, val):
        assert numeric_real(val)
        self._x1 = val

    x1 = property(get_x1, set_x1)

    def get_y1(self):
        return self._y1

    def set_y1(self, val):
        assert numeric_real(val)
        self._y1 = val
       
    y1 = property(get_y1, set_y1)
    
    def get_x2(self):
        return self._x2

    def set_x2(self, val):
        assert numeric_real(val)
        self._x2 = val

    x2 = property(get_x2, set_x2)

    def get_y2(self):
        return self._y2

    def set_y2(self, val):
        assert numeric_real(val)
        self._y2 = val
       
    y2 = property(get_y2, set_y2)

    def get_p1(self):
        return Point2D(self._x1, self._y1)

    def set_p1(self, p):
        assert isinstance(p, Point2D)
        self.set_x1(p.x)
        self.set_y1(p.y)

    p1 = property(get_p1, set_p1)

    def get_p2(self):
        return Point2D(self._x2, self._y2)

    def set_p2(self, p):
        assert isinstance(p, Point2D)
        self.set_x2(p.x)
        self.set_y2(p.y)

    p2 = property(get_p2, set_p2)

    def points(self):
        return (self.p1, self.p2)

    def norm(self):
        return math.sqrt(self.norm2())

    def norm2(self):
        return (self._x2-self._x1)**2 + (self._y2-self._y1)**2
  
    def p1_tuple(self):
        return (self._x1, self._y1)

    def p2_tuple(self):
        return (self._x2, self._y2)
        
    def __repr__(self):
        return repr(self.p1) + " -> " + repr(self.p2)
        
    def intersects(self, oline):
        A = self.p1
        B = self.p2
        C = oline.p1
        D = oline.p2
        
        # non li becca tutti, vedi http://compgeom.cs.uiuc.edu/~jeffe/teaching/373/notes/x06-sweepline.pdf
        def ccw(A,B,C):
            return (C.y-A.y)*(B.x-A.x) > (B.y-A.y)*(C.x-A.x)

        return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

    def angle(self):
        dx = self._x2 - self._x1
        dy = self._y2 - self._y1
        return normalize_angle(math.atan2(dy, dx))


class Rectangle():
    """
    A class representing a 2D rectangle
    """
    def __init__(self, pos, rot, width, height):
        assert isinstance(pos, Point2D)
        assert isinstance(rot, int) or isinstance(rot, float)
        assert isinstance(width, int) or isinstance(width, float)
        assert width > 0.
        assert isinstance(height, int) or isinstance(height, float)
        assert height > 0.

        self._pos = pos
        self._rot = rot
        self._width = width
        self._height = height
        self._area = width*height
        self._lines = []
        self._points = [pos.tuple()]
        self._xlim = (pos.x, pos.x)
        self._ylim = (pos.y, pos.y)

        if self._area:
            rotmat = rotation_matrix(self._rot)
            points = numpy.array([
                 (   0.,     0.),
                 (width,     0.),
                 (width, height),
                 (   0., height),
            ]) 
        
            # rotate the rectangle points and tranlate them to the rectangle origin
            points = transform_points(points, self._rot, pos.array())

            #print points, points[:,0]
            x1, y1 = pos.x, pos.y
            x2, y2 = points[1,:]
            x3, y3 = points[2,:]
            x4, y4 = points[3,:]
        
            self._xlim = (numpy.min(points[:,0]), numpy.max(points[:,0]))
            self._ylim = (numpy.min(points[:,1]), numpy.max(points[:,1]))
        
            p1 = pos
            p2 = Point2D(x2, y2)
            p3 = Point2D(x3, y3)
            p4 = Point2D(x4, y4)
            self._points.extend([(x2, y2), (x3, y3), (x4, y4)])
            self._points.extend([p2.tuple(), p3.tuple(), p4.tuple()])
            self._lines.append(Line2D(p1, p2))
            self._lines.append(Line2D(p2, p3))
            self._lines.append(Line2D(p3, p4))
            self._lines.append(Line2D(p4, p1))
                
    def area(self):
        return self._area

    def lines(self):
        return self._lines
        
    def xlim(self):
        return self._xlim

    def ylim(self):
        return self._ylim
        
    def points(self):
        return self._points

    def contains_point(self, p):
        assert isinstance(p, Point2D)

        dxy = Point2D(p.x - self._pos.x, p.y - self._pos.y)
        loc = rotate_point(dxy, -self._rot)
        return loc.x >= 0. and loc.x <= self._width and loc.y >= 0. and loc.y <= self._height

    def contains_line(self, line):
        assert isinstance(line, Line2D)
        p1, p2 = line.points()
        return self.contains_point(p1) and self.contains_point(p2)

    def contains_rect(self, rect):
        assert isinstance(rect, Rectangle)
        for line in rect.lines():
            if not self.contains_line(line):
                return False

        return True
    
    # def plot(self, fig, axis):
    #     # namespace shortcut for the codes below
    #     Path = mpath.Path

    #     rotmat = rotation_matrix(self._rot)

    #     w, h = self._width, self._height
    #     points = numpy.array([
    #          (0, 0),
    #          (w, 0),
    #          (w, h),
    #          (0, h),
    #          (0, 0),
    #     ]) 
        
    #     # rotate the badge points and tranlate them to the camera origin
    #     verts = points.dot(rotmat) + numpy.array([self._posx, self._posy])
    #     codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY]
    #     patch = mpatches.PathPatch(mpath.Path(verts, codes), alpha=0.75, facecolor=(0.75,0.75,0.75))
    #     axis.add_patch(patch)


