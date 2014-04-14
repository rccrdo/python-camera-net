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
import numpy
from optparse import OptionParser
from hypo import HypoTest
from transition import TransitionalModel#HiddenMarkovModel, HiddenMarkovModelVNode, HiddenMarkovModelVNodeHypo


if __name__ == "__main__":

    if 1:
        print "warning, suppressing debug info"
        numpy.seterr(all='ignore')

    parser = OptionParser()

    parser.add_option("-v", "--vnodes", dest="vnodes", action="store_true", 		\
                      help="Use null nodes", default=False)

    parser.add_option("-y", "--hypo", dest="hypo", action="store_true", 		\
                      help="Use hypothesis test", default=False)

    parser.add_option("--max", dest="max_transitions", default=1000,			\
                      help="The max. number of transitions to use when building the models.")

    #parser.add_option("--coverage", dest="only_coverage", action="store_true", 		\
    #                  help="only plot the initial coverage model", default=False)

    parser.add_option("-s", "--save-plots", dest="saveplots", action="store_true", 		\
                      help="save all plots", default=False)

    (options, args) = parser.parse_args()

    # after parsing arguments we shold be left with the path to the
    # training and validation data sets
    datapaths = args
    if not len(datapaths):
        print "error: no paths given for either training or validation data\n"
        parser.print_help()
        sys.exit(0)
    else:
        trainpath = datapaths[0]
        validpaths = datapaths[1:]

    try:
        options.max_transitions = int(options.max_transitions)
        if not options.max_transitions >= 0:
            raise
    except:
        print "error: --max argument accepts only a non-negative integer\n"
        parser.print_help()
        sys.exit(0)

    hypotest = None
    if options.hypo:
        bargamma = 0.1
        alpha_0 = 0.05
        hypotest = HypoTest(bargamma, alpha_0)

    model = TransitionalModel(trainpath, validpaths, options.max_transitions,	\
                              vnodes=options.vnodes, hypotest=hypotest, saveplots=options.saveplots)
       
    #if not options.only_coverage:
    model.optimize()

    print "Press any key to exit ..."
    raw_input()
    os.sys.exit(0)

