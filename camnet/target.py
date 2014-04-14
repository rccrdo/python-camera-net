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
import matplotlib.pyplot as plt
import networkx as nx
from math2D import *

class Target():
    _target_id = 0
    def __init__(self, pos=Point2D(1,1), targetid = None):
        assert isinstance(pos, Point2D)
        self._pos = pos
        self._traj = []

        # set the target id
        Target._target_id += 1
        if targetid:
            self._id = targetid
        else:
            self._id = Target._target_id


    def step(self, time, dt, walls=[]):
        assert numeric_real(time)
        assert numeric_real(dt)
        assert dt > 0
        
        r = 3
        f = 0.01
        self._pos.x = 9.5/2. + 3.5*math.cos(time*f + math.pi)
        self._pos.y = 3 + 2.5*math.sin(time*f  + math.pi)

    def pos(self):
        return self._pos

    def id(self):
        return self._id

    def plot(self, axis):
        pos = self._pos.tuple()
        if 0:
            axis.add_patch(matplotlib.patches.Circle(pos, radius=0.1, alpha=0.5))
        
        for line in self._traj:
            assert isinstance(line, Line2D)
            axis.add_line(matplotlib.lines.Line2D([line.p1.x, line.p2.x],[line.p1.y, line.p2.y], color=(1,0.,0.), alpha=1))#, zorder=-100))

        
class RandomTarget(Target):
    def __init__(self, pos):
        Target.__init__(self, pos)

    def step(self, cur_time, step, walls=[]):
        assert numeric_real(time)
        assert numeric_real(step)
        assert step > 0

        loop = True
        while loop:
            old_point = Point2D(self._pos.x, self._pos.y) 
            dx, dy = numpy.random.normal(0, 0.25, 2)
            newx = numpy.clip(self._pos.x + dx, 0, 9)
            newy = numpy.clip(self._pos.y + dy, 0, 6)
            new_point = Point2D(newx, newy)
            line = Line2D(old_point, new_point)

            loop = False
            for wall in walls:
                if line.intersects(wall):
                    #print "RandomTarget intersected wall ", wall
                    # ops we bumped into a wall, retry :)
                    loop = True
                    break

        self._pos.x = newx
        self._pos.y = newy



class GraphTargetBase(Target):
    def __init__(self, targetid=None):
        Target.__init__(self, Point2D(0,0), targetid)

        self._graph = nx.Graph()
	self._graph.position = {}
	
        self._cur_node = None
        self._target_pos = None
        self._moving = False

    def step(self, cur_time, step, walls=[]):
        assert numeric_real(cur_time)
        assert numeric_real(step)
        assert step > 0

        if self._moving:
            VEL = 0.075
            STEP = VEL*step
            cur_pos = self._pos
            line = Line2D(cur_pos, self._target_pos)
            if line.norm() < STEP:
                newx = self._target_pos.x
                newy = self._target_pos.y
                self._moving = False
            else:
                dx = self._target_pos.x - cur_pos.x
                dy = self._target_pos.y - cur_pos.y
                dx = dx*(STEP/line.norm())
                dy = dy*(STEP/line.norm())
                newx = cur_pos.x + dx
                newy = cur_pos.y + dy
            
            self._pos.x = newx
            self._pos.y = newy        
        else: 
            self.plan(walls)
            self._moving = True
            self.step(cur_time, step, walls)


    def plot(self, axis):
        # Plot the target badge and trajectory first
        Target.plot(self, axis)

        # debug the transition graph
        if 0:
            #node_pos = []
            #for v in self._graph.nodes():
            #    p = self._graph.position[v]
            #    node_pos.append(p)
            node_pos = [self._graph.position[v] for v in self._graph.nodes()]
            nx.draw_networkx_edges(self._graph, self._graph.position, self._graph.edges(), edge_color='y', alpha=0.25, ax=axis)  
            nx.draw_networkx_nodes(self._graph, self._graph.position, self._graph.nodes(), 200, node_color='r', ax=axis)
            nx.draw_networkx_labels(self._graph, self._graph.position, ax=axis)




class MarkovTarget(GraphTargetBase):
    def __init__(self, targetid=None):
        GraphTargetBase.__init__(self, targetid)

        # build the transition graph
        self._graph.add_node(1)
        self._graph.position[1] = (7.75, 5.25)

        self._graph.add_node(2)
        self._graph.position[2] = (6.5, 5.25)

        self._graph.add_node(3)
        self._graph.position[3] = (5.75, 4)

        self._graph.add_node(4)
        self._graph.position[4] = (5, 4.75)

        self._graph.add_node(5)
        self._graph.position[5] = (3, 5.25)

        self._graph.add_node(6)
        self._graph.position[6] = (1.75, 5.5)

        self._graph.add_node(7)
        self._graph.position[7] = (1.5, 4.75)

        self._graph.add_node(8)
        self._graph.position[8] = (1.75, 3)

        self._graph.add_node(10)
        self._graph.position[10] = (1.5, 1.25)

        self._graph.add_node(11)
        self._graph.position[11] = (3, 1.)

        self._graph.add_node(12)
        self._graph.position[12] = (4, 2)

        self._graph.add_node(13)
        self._graph.position[13] = (4.5, 1)

        self._graph.add_node(14)
        self._graph.position[14] = (5.75, 2)

        self._graph.add_node(15)
        self._graph.position[15] = (7, 1.)

        self._graph.add_node(16)
        self._graph.position[16] = (8, 1.25)

        self._graph.add_node(17)
        self._graph.position[17] = (8.25, 2)

        self._graph.add_node(18)
        self._graph.position[18] = (7.5, 4.)

        self._graph.add_edge(1,2)
        self._graph.add_edge(2,3)
        self._graph.add_edge(2,4)
        self._graph.add_edge(3,4)
        self._graph.add_edge(4,5)
        self._graph.add_edge(4,7)
        self._graph.add_edge(5,6)
        self._graph.add_edge(5,7)
        self._graph.add_edge(6,7)
        self._graph.add_edge(7,8)
        self._graph.add_edge(7,10)
        self._graph.add_edge(8,10)
        self._graph.add_edge(10,11)
        self._graph.add_edge(11,12)
        self._graph.add_edge(11,13)
        self._graph.add_edge(12,13)
        self._graph.add_edge(12,14)
        self._graph.add_edge(13,14)
        self._graph.add_edge(13,15)
        self._graph.add_edge(14,3)
        self._graph.add_edge(14,13)
        self._graph.add_edge(14,15)
        self._graph.add_edge(15,16)
        self._graph.add_edge(15,17)
        self._graph.add_edge(16,17)
        self._graph.add_edge(17,18)
        self._graph.add_edge(17,1)
        self._graph.add_edge(18,1)
        
        self._cur_node = 10
        self._pos = Point2D(*self._graph.position[self._cur_node])

    def plan(self, walls):
        loop = True
        old_point = self._pos
        neighbors = self._graph[self._cur_node].keys()
        while loop:
            # select the next node
            next = neighbors[numpy.random.randint(len(neighbors))]
            xc, yc = self._graph.position[next]
                
            # 3 and 14 are the nodes at the entry/exit of the passage
            # we use a smaller variance to avoid bumping into passage lateral 
            # wall of
            sigma2 = 0.175
            if next in [3,14]:
                sigma2 = 0.1

            for i in xrange(0,10):
                dx, dy = numpy.random.normal(0, sigma2, 2)
                newx = numpy.clip(xc + dx, 0, 9)
                newy = numpy.clip(yc + dy, 0, 6)
                new_point = Point2D(newx, newy)
                line = Line2D(old_point, new_point)
                self._target_pos = new_point

                # check if the new segment in the target trajectory
                # intersects any walls
                hit = False
                for wall in walls:
                    if line.intersects(wall):
                        # ops we bumped into a wall, retry :)
                        hit = True
                        break
                if not hit:
                    self._cur_node = next
                    self._traj.append(Line2D(old_point, self._target_pos))
                    loop = False
                    break




class EightTarget(GraphTargetBase):
    def __init__(self, targetid=None):
        GraphTargetBase.__init__(self, targetid)

        # build the eight shaped transition graph
        self._graph.add_node(1)
        self._graph.position[1] = (5.75, 4.)

        self._graph.add_node(2)
        self._graph.position[2] = (4.5, 4.75)

        self._graph.add_node(3)
        self._graph.position[3] = (3., 4.95)

        self._graph.add_node(4)
        self._graph.position[4] = (1.5, 4.75)

        self._graph.add_node(5)
        self._graph.position[5] = (1.5, 3)

        self._graph.add_node(6)
        self._graph.position[6] = (1.5, 1.5)

        self._graph.add_node(7)
        self._graph.position[7] = (2.5, 1)

        self._graph.add_node(8)
        self._graph.position[8] = (3.5, 1.25)

        self._graph.add_node(9)
        self._graph.position[9] = (3.75, 2)

        self._graph.add_node(10)
        self._graph.position[10] = (5.75, 2)

        self._graph.add_node(11)
        self._graph.position[11] = self._graph.position[1]

        self._graph.add_node(12)
        self._graph.position[12] = (6, 5.25)

        self._graph.add_node(13)
        self._graph.position[13] = (7.5, 5.25)

        self._graph.add_node(14)
        self._graph.position[14] = (7.75, 3)

        self._graph.add_node(15)
        self._graph.position[15] = (8, 1.5)

        self._graph.add_node(16)
        self._graph.position[16] = (7, 1.25)

        self._graph.add_node(17)
        self._graph.position[17] = (6, 1.25)
        
        self._graph.add_node(18)
        self._graph.position[18] = self._graph.position[10]
       
        self._graph.add_edge(1,2)
        self._graph.add_edge(2,3)
        self._graph.add_edge(3,4)
        self._graph.add_edge(4,5)
        self._graph.add_edge(5,6)
        self._graph.add_edge(6,7)
        self._graph.add_edge(7,8)
        self._graph.add_edge(8,9)
        self._graph.add_edge(9,10)
        self._graph.add_edge(10,11)
        self._graph.add_edge(11,12)
        self._graph.add_edge(12,13)
        self._graph.add_edge(13,14)
        self._graph.add_edge(14,15)
        self._graph.add_edge(15,16)
        self._graph.add_edge(16,17)
        self._graph.add_edge(17,18)
        self._graph.add_edge(18,1)
        
        self._cur_node = 10
        self._pos = Point2D(*self._graph.position[self._cur_node])

    def plan(self, walls):
        loop = True
        old_point = self._pos
        neighbors = self._graph[self._cur_node].keys()
        while loop:
            # select the next node
            if self._cur_node + 1 in neighbors:
                next = self._cur_node + 1 
            else:
                assert 1 in neighbors
                next = 1
                
            xc, yc = self._graph.position[next]
                
            # 10 11 17 1 are the nodes at the entry/exit of the passage
            # we use a smaller variance to avoid bumping into the lateral walls
            # of the passage
            sigma2 = 0.175
            if next in [1,10,11,18]:
                sigma2 = 0.1
            for i in xrange(0,10):
                dx, dy = numpy.random.normal(0, sigma2, 2)
                newx = numpy.clip(xc + dx, 0, 9)
                newy = numpy.clip(yc + dy, 0, 6)
                new_point = Point2D(newx, newy)
                line = Line2D(old_point, new_point)
                self._target_pos = new_point

                # check if the new segment in the target trajectory
                # intersects any walls
                hit = False
                for wall in walls:
                    if line.intersects(wall):
                        # ops we bumped into a wall, retry :)
                        hit = True
                        break
                if not hit:
                    self._cur_node = next
                    self._traj.append(Line2D(old_point, self._target_pos))
                    loop = False
                    break

                print "Target traj. planning, discarding segment:", Line2D(old_point, self._target_pos)
                print " cur_node, next_node:", self._cur_node, next





