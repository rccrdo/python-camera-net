#! /usr/bin/python

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

import os
import sys
import time
import math
import shutil

from observation import ObservationAccu
from math2D import Point2D
from environment import Environment
from camera import CameraFixed, CameraPTZ
from target import Target, RandomTarget, MarkovTarget, EightTarget
from geometry import Building


# define the cameras
fixcam_radius = 3
ptzcam_radius = 2.35
ptzcam_fullfov = math.pi/180.*10.*12
ptzcam_fov = math.pi/180.*10.
ptzcam_vel = 0.10
cam1_pos = Point2D(0.546, 3.341)
cam2_pos = Point2D(7, 0.546)
cam3_pos = Point2D(3.461, 5.476)
cam4_pos = Point2D( 0.33,  0.33)
cam5_pos = Point2D(  2.9,  0.25)
cam6_pos = Point2D( 8.75,  5.75)

if 0:
    # all bisecting lines
    cam1_bpoints = []#[Point2D(2,1.5)]
    cam2_bpoints = [Point2D(5.5, 3.3)]#[Point2D(6.5,1.5), Point2D(7.5,1.5), Point2D(6,2.5)]
    cam3_bpoints = [Point2D(4.6,3.5), Point2D(6,3.)]#[Point2D(2,4.5), Point2D(4,4.5), Point2D(5.5,3.5)]
    cam4_bpoints = [Point2D(5.36,2.5), Point2D(6.,2.69)]
    cam5_bpoints = [Point2D(3.1,2.5)]
    cam6_bpoints = []
else:
    # bisecting-lines used in simulations
    cam1_bpoints = []
    cam2_bpoints = [Point2D(5.5, 3.3)]#[Point2D(6.5,1.5), Point2D(7.5,1.5), Point2D(6,2.5)]
    cam3_bpoints = [Point2D(4.6,3.5), Point2D(6,3.)]#[Point2D(2,4.5), Point2D(4,4.5), Point2D(5.5,3.5)]
    cam4_bpoints = [Point2D(6.,2.69)]
    cam5_bpoints = [Point2D(3.1,2.5)]
    cam6_bpoints = []

cam1_fovpoints = [Point2D(2,1.5), Point2D(2,4.5)]
cam2_fovpoints = [Point2D(7.5,1.5), Point2D(3,1.5)]
cam3_fovpoints = [Point2D(2,4.5), Point2D(6.5,5)]
cam4_fovpoints = [Point2D(2,4.5), Point2D(6.5,1.5)]
cam6_fovpoints = [Point2D(7.5,2.5), Point2D(7,5)]
cam5_fovpoints = [Point2D(2.35,1.5), Point2D(3.85,2.5)]

def test_get_ptz_cameras(use_blines=True):
    _cam1_bpoints = cam1_bpoints
    _cam2_bpoints = cam2_bpoints
    _cam3_bpoints = cam3_bpoints
    _cam4_bpoints = cam4_bpoints
    _cam5_bpoints = cam5_bpoints
    _cam6_bpoints = cam6_bpoints

    if not use_blines:
        _cam1_bpoints = []
        _cam2_bpoints = []
        _cam3_bpoints = []
        _cam4_bpoints = []
        _cam5_bpoints = []
        _cam6_bpoints = []
    
    return [
                CameraPTZ(  cam1_pos, -math.pi*0.0385, 2.5, math.pi*0.496, ptzcam_fov, ptzcam_vel, '1', \
                          _cam1_bpoints, cam1_fovpoints),
                CameraPTZ(  cam2_pos,  math.pi*0.635, 4.5, math.pi*0.58, ptzcam_fov, ptzcam_vel, '2',	\
                          _cam2_bpoints, cam2_fovpoints),
                CameraPTZ(  cam3_pos,  -math.pi*0.425, 3.5, math.pi*0.755, ptzcam_fov, ptzcam_vel, '3',	\
                	  _cam3_bpoints, cam3_fovpoints),
                CameraPTZ(  cam4_pos,   math.pi*0.218, 6.5, math.pi*0.317, ptzcam_fov, ptzcam_vel, '4',	\
                          _cam4_bpoints, cam4_fovpoints),
                CameraFixed(cam5_pos,      math.pi/2., 2.5,    math.pi/4., '5', _cam5_bpoints,		\
                	    cam5_fovpoints),
                CameraFixed(cam6_pos,   -2.975/4.*math.pi,   4,  math.pi/3.8, '6', _cam6_bpoints,	\
                           cam6_fovpoints)
          ]

def test_get_fixed_cameras(use_blines=True):
    _cam1_bpoints = cam1_bpoints
    _cam2_bpoints = cam2_bpoints
    _cam3_bpoints = cam3_bpoints
    _cam4_bpoints = cam4_bpoints
    _cam5_bpoints = cam5_bpoints
    _cam6_bpoints = cam6_bpoints

    if not use_blines:
        _cam1_bpoints = []
        _cam2_bpoints = []
        _cam3_bpoints = []
        _cam4_bpoints = []
        _cam5_bpoints = []
        _cam6_bpoints = []

    fixcam_fov = math.pi/3.
    return [
                CameraFixed(cam1_pos,    -math.pi/32., fixcam_radius, fixcam_fov, '1', _cam1_bpoints),
                CameraFixed(cam2_pos,  9.*math.pi/16., fixcam_radius, fixcam_fov, '2', _cam2_bpoints),
                CameraFixed(cam3_pos,     -math.pi/2., fixcam_radius, fixcam_fov, '3', _cam3_bpoints),
                CameraFixed(cam4_pos,      math.pi/4., fixcam_radius, fixcam_fov, '4', _cam4_bpoints),
                CameraFixed(cam5_pos,      math.pi/2., fixcam_radius, math.pi/4., '5', _cam5_bpoints),
                CameraFixed(cam6_pos, -11.*math.pi/16, fixcam_radius, math.pi/3., '6', _cam6_bpoints),
          ]

    #ptzcam_radius = 5
    #return [
    #            CameraFixed(cam1_pos,            0., ptzcam_radius, math.pi/4., '1', cam1_bpoints),
    #            CameraFixed(cam2_pos,     math.pi/2, ptzcam_radius, math.pi/4., '2', cam2_bpoints),
    #            CameraFixed(cam3_pos,   -math.pi/2., ptzcam_radius, math.pi/4., '3', cam3_bpoints),
    #            CameraFixed(cam4_pos,    math.pi/4., ptzcam_radius, math.pi/4., '4'),
    #            CameraFixed(cam5_pos,    math.pi/2., fixcam_radius, math.pi/4., '5', cam5_bpoints),
    #            CameraFixed(cam6_pos, -3./4*math.pi, fixcam_radius, math.pi/4., '6'),
    #      ]


if __name__ == "__main__":
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("-t", "--time", dest="sim_time",
                      help="simulation time in seconds. By default 10000 s",
                      default=10000, metavar="duration")
    parser.add_option("-s", "--skip", dest="skip_time",
                      help="do plot only every skip_time seconds", default=10,
                      metavar='skip_time')
    parser.add_option("--target", dest="target",
                      help="select the target implementation: markov (default) or 8", default="markov")
    parser.add_option("--no-blines", dest="use_blines", action="store_false",
                      help="do not use bisecting lines", default=True)
    parser.add_option("-f", "--fixcam", dest="fixcam", action="store_true",
                      help="Make all cameras fixed", default=False)
    parser.add_option("-v", "--video", dest="video", default=False, action="store_true", 	\
                      help="Make a video")

    (options, args) = parser.parse_args()

    if len(args):
        print "error: cannot parse these arguments %s\n" % args
        parser.print_help()
        sys.exit(0)

    try:        
        options.sim_time = float(options.sim_time)
        if options.sim_time <= 0:
            raise
    except:
        print "error: -t/--time arguments require a (positive) numeric value\n"
        parser.print_help()
        sys.exit(0)

    try:
        options.skip_time = math.fabs(float(options.skip_time))
    except:
        print "error: -s/--skip arguments require a (non-negative) numeric value\n"
        parser.print_help()
        sys.exit(0)

    if not options.target in ('markov', '8'):
        print "error: unknown value \"%s\" for argument --target\n" % options.target
        parser.print_help()
        sys.exit(0)

    # create the environment and add the building and cameras
    env = Environment(9.5, 6)
    building = Building()
    building.set_passage_open(True)
    env.add_geometry(building)
    if options.fixcam:
        cameras = test_get_fixed_cameras(options.use_blines)
    else:
        cameras = test_get_ptz_cameras(options.use_blines)
    env.add_cameras(cameras)

    # select target
    if options.target == 'markov':
        if 1:
            env.add_targets([MarkovTarget(targetid='target-hmm')])
    elif options.target == '8':
        env.add_targets([EightTarget(targetid='target-hmm')])     

    
    if options.video:
        # Create dirs in which frames and video will be saved
        simid_base = 'video-t=%d-s=%.2f-blines=%s-fixcam=%s-target=%s' 	\
        	% (options.sim_time, options.skip_time,			\
        	   repr(options.use_blines), repr(options.fixcam),	\
        	   options.target)
        	   
        sim_dir = os.path.join("./", simid_base)
        frames_dir = '%s' % sim_dir
        video_path = '%s/%s.avi' % (sim_dir, simid_base)
        video_cmd = "avconv -y -r 2 -i %s/%%.4d.jpg -b 2000k %s" % (frames_dir, video_path)

        if os.path.exists(sim_dir):
            shutil.rmtree(sim_dir)

        os.mkdir(sim_dir)
        #os.mkdir(frames_dir)

    # simulation step in [s]
    sim_step = 1.
    next_plot_time = options.skip_time
    accu = ObservationAccu()
    print "Simulation time: %.2fs" % options.sim_time
    i = 0
    while True:
        i += 1

        # step the simulation
        env.step(sim_step)

        cur_time = env.time()
        print " %*.2f %% " % (6, cur_time/options.sim_time *100)

        # gather observations for target 'target-hmm'
        new_observ = env.get_observation_for_target('target-hmm')
        accu.accumulate(new_observ)

        # do plots 
        if next_plot_time <= cur_time:
            next_plot_time += options.skip_time
            if options.video:
                env.plot(savepath=frames_dir + "/%.4d.png" % i)
            else:
                env.plot()

        # exit the sim. loop if we are past the choosen sim. lenght
        if cur_time >= options.sim_time:
            break

    # Save observation data to file
    filename = '-'.join(['simdata', options.target + ("-noblines", "")[options.use_blines] 	\
    			+ ("", "-fixcam")[options.fixcam], time.strftime("%Y.%m.%d-%H.%M")])
    savepath = os.path.join(os.getcwd(), filename) 
    accu.save(savepath)    

    if options.video:
        # Call avconv to make a video from the saved frames
        print video_cmd
        print video_cmd.replace("=", "\=")
        print video_cmd
        os.system(video_cmd.replace("=", "\="))

    sys.exit (0)

