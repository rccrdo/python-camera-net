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
import copy
import numpy

import matplotlib
from matplotlib import rc
rc('font', family='sans-serif')
import matplotlib.pyplot as plt
import networkx as nx

from util import *
from observation import ObservationSequence
from hmm import HiddenMarkovModel

class TransitionalModel():
    def __init__(self, trainpath, validpaths, max_trans=0, vnodes=False, hypotest=None, saveplots=False):
        assert isinstance(trainpath, str)
        assert len(trainpath)
        assert isinstance(validpaths, list)
        for path in validpaths:
            assert isinstance(path, str)
            assert len(path)
        assert isinstance(max_trans, int)
        assert max_trans >= 0
        assert isinstance(vnodes, bool)
        assert isinstance(saveplots, bool)

        self._modname = 'base'
        if vnodes:
            self._modname = 'vnode'
        if hypotest:
            self._modname += '+hypo'
        
        self._savepathbase = None
        if saveplots:
            self._savepathbase = "./%s-plots-%s-maxtrans%d/" % (trainpath, self._modname, max_trans)
            os.mkdir(self._savepathbase)

        # load datasets
        self._trainset = ObservationSequence(loadpath=trainpath, max_trans=max_trans, \
                                             vnodes=vnodes)

        if hypotest:
            #statespace, transitions = self._trainset.admissible_statespace_transitions()
            #print "orig:\n", statespace, transitions
            statespace, transitions = self._trainset.hypotest_admissible_statespace_transitions(hypotest)
            #print "hypotest:\n", statespace, transitions
            self._trainset = ObservationSequence(loadpath=trainpath, max_trans=max_trans, vnodes=vnodes, \
                                             admissible_statespace=statespace, \
                                             admissible_transitions=transitions)

        statespace, transitions = self._trainset.admissible_statespace_transitions()
        #self._trainset.dump()
        #print "trainset:\n", statespace, transitions

        if not self._trainset.outseq_len():
            print "warning, TransitionModel.__init__, empty trainset at path %s" % trainpath


        _FILTER_TRANS = 1
        self._validsets = []
        for validpath in validpaths:
            if _FILTER_TRANS:
                # filter trans
                dataset = ObservationSequence(loadpath=validpath, max_trans=max_trans, vnodes=vnodes, \
                                              admissible_statespace=statespace, \
                                              admissible_transitions=transitions)
            else:
                dataset = ObservationSequence(loadpath=validpath, max_trans=max_trans, vnodes=vnodes)

            if not dataset.outseq_len():
                print "TransitionModel, skipping empty validation dataset at path %s" % validpath
                continue
            self._validsets.append(dataset)

        # define the initial set of states and the set of output symbols
        self._init_states = sorted(self._trainset.states())
        self._N = len(self._init_states)
        self._symbols = sorted(self._trainset.symbols())
        self._M = len(self._symbols)

        self._symbols_ext = copy.copy(self._symbols)        
        if not _FILTER_TRANS:

              symbols = set()
              for dataset in self._validsets:
                  symbols.update(dataset.symbols())
                  
              print "symbols", symbols
              for sym in symbols:
                  if sym not in self._symbols_ext:
                      print "sym", sym
                      print "TransitionalModel: appending sym=\"%s\" in valset but not in trainset" % sym
                      self._symbols_ext.append(sym)
        self._M_ext = len(self._symbols_ext)
        
        print "self._symbols_ext", self._symbols_ext
        camsymbols = {}
        for sym in self._symbols_ext:
            print "sym = ", sym
            for _id in sym.split('-'):
                print "  _id = ", _id

                __id = int(_id.strip('abcdef'))
                try:
                    camsymbols[__id].add(sym)
                except:
                    print "new cam id", _id
                    camsymbols[__id] = set([sym])
        for camid in camsymbols.keys():
            symbols = camsymbols[camid]
            camsymbols[camid] = sorted(symbols)

        symbol2cam = {}
        for sym in self._symbols_ext:
            for _id in sym.split('-'):
                __id = int(_id.strip('abcdef'))
                try:
                    symbol2cam[sym].append(__id)
                except:
                    symbol2cam[sym] = [__id]
      
        symbol2cam_enc = {}
        for sym,camids in symbol2cam.iteritems():
            symenc = []
            symbol2cam_enc[self._symbol_enc(sym)] = camids

        print "camsymbols=",camsymbols
        print "symbol2cam=",symbol2cam        
        print "symbol2cam_enc=",symbol2cam_enc
        print "symbols_ext", self._symbols_ext
                
        # define the initial A, B matrices
        A, B = self._build_initial_model()

        # compute the node to be split
        self._compute_node_splits()
       
        # encode datasets (bijective mapping symobols <-> integers)
        self._trainset_enc = []
        for observ in self._trainset.outseq():
            self._trainset_enc.append(self._symbol_enc(observ))

        self._validsets_enc = []
        for dataset in self._validsets:
            encoded = []
            for observ in dataset.outseq():
                encoded.append(self._symbol_enc(observ))
            self._validsets_enc.append(encoded)

        # print some info
        self.print_info()

        # create the HMM model with node splitting
        self._hmm = HiddenMarkovModel(A, B, copy.copy(self._init_states),	\
              			      self._trainset_enc, validsets=self._validsets_enc, symnames=copy.copy(self._symbols_ext), symbol2cam=symbol2cam_enc)

        self._transition_graph_threshold = 0.01
        self._hmm.plot_transition_graph(0., savepathbase=self._savepathbase, savetag='coverage-model')

    def optimize(self, nu=0.01):
        # plot coverage
        self._hmm.plot_transition_graph(self._transition_graph_threshold,	\
                			savepathbase=self._savepathbase, savetag='coverage-model')

        # run BW before the first node-split
        while True:
            done = self._hmm.baumwelch(nu)
            self._hmm.plot_perf(savepathbase=self._savepathbase)
            self._hmm.plot_transition_graph(self._transition_graph_threshold,	\
                				savepathbase=self._savepathbase)
            if done:
                break

        while True:

            keys = sorted(self._pending_splits.keys())
            print "\n%d node-splits left" % len(keys)
            if not len(keys):
                break

            key = keys[0]
            splits = sorted(self._pending_splits[key], key=lambda split_data: split_data[0])
            #print "splits for node ", key, " (", len(splits), "):", splits
            assert len(splits)

            split = splits.pop()

            if not len(splits):
                self._pending_splits.pop(key)
            else:
                self._pending_splits[key] = splits

            # per ogni nodo corrente
            #   guarda alla backward star
            #   guarda alla forward star 
            #   costruisci matrice rettangolare transizioni
            #   identifica la riga piu' ortogonale a tutte le altre
            #
            #     se lo trovi splitta il nodo ribalanciando le probabilita'
           #(inner, bsym, bstar_syms, fstar_syms, row)
            #print split
            _inner = split[0]
            _bsym = split[1]
            _bstar_sym = split[2]
            _fstar_sym = split[3]
            rates = split[4]
            assert numpy.isclose(rates.sum(), 1)

            #print "   splitting node %s : inner=%.2f bsym=%s _bstar=%s _fstar=%s, rates=%s" % (key, _inner, _bsym, _bstar_sym, _fstar_sym, rates)

            state = self._state_enc(key)
            bstate = self._state_enc(_bsym)
            bstar = []
            for i in _bstar_sym:
                bstar.append(self._state_enc(i))
            fstar = []
            for i in _fstar_sym:
                fstar.append(self._state_enc(i))

            self._hmm.split_state(state, bstate, bstar, fstar, rates)
            self._hmm.plot_transition_graph(self._transition_graph_threshold,	\
            				    savepathbase=self._savepathbase, savetag='node-splitting')

            while True:
                done = self._hmm.baumwelch(nu)
                self._hmm.plot_perf(savepathbase=self._savepathbase)
                self._hmm.plot_transition_graph(self._transition_graph_threshold,	\
                				savepathbase=self._savepathbase)
                if done:
                    break

    def _symbol_enc(self, sym):
        #print sym
        #print self._symbols
        assert sym in self._symbols_ext
        return self._symbols_ext.index(sym)
        
    def _state_enc(self, state):
        assert state in self._init_states
        return self._init_states.index(state)

    def name(self):
        return self._modname

    def nr_training_transitions(self):
        return max(0, self._trainset.outseq_len() -1)

    def nr_validation_transitions(self, setid):
        assert isinstance(setid, int)
        assert setid >= 0
        assert setid < self._validsets.outseq_len()
        return max(0, self._validsets[setid].outseq_len() -1)

    def _build_initial_model(self):
        dataset = self._trainset
        A = numpy.zeros((self._N, self._N))

        for src in self._init_states:
            nr_from = dataset.count_transitions_from(src)
            #print "nr_from: %s = %d" % (src, nr_from)
            assert nr_from >= 0
            if nr_from == 0:
                continue

            src_idx = self._state_enc(src)

            for dest in self._init_states:
                nr_to = dataset.count_transitions_from_to(src, dest)
                #print "nr_from_to: %s->%s = %d" % (src, dest, nr_from)
                assert nr_to >= 0
                if nr_to == 0:
                    continue

                dest_idx = self._state_enc(dest)
                A[src_idx,dest_idx] = nr_to/nr_from

        for i in xrange(self._N):
            # potremmo non avere transizioni uscenti da un certo
            # nodo ed in quel caso le righe di A non sommano tutte all'unita'
            # usa una distribuzione uniforme
            if not numpy.isclose(A[i,:].sum(), 1.):
                A[i,:] = numpy.ones(self._N)/self._N
       
        B = numpy.zeros((self._M,self._N))
        for state in self._init_states:
            if state in self._symbols:
                B[self._symbol_enc(state), self._state_enc(state)] = 1.
            else:
                # this must be one of the vnull states
                assert '~0~' in state
                B[self._symbol_enc('0'), self._state_enc(state)] = 1.

        return A, B

    def _compute_node_splits(self):
        triplets = self._trainset.state_triplets_dict()
        fstar = self._trainset.state_fstar_dict()
        bstar = self._trainset.state_bstar_dict()

        #print "fstar: ", fstar
        #print "bstar: ", bstar

        splits = {}     # computed splits keyed by the id of the splitted node

        # compute splits for each state
        states = tuple(set(self._trainset.states()))
        for state in states:
            if (state not in fstar.keys()) or (state not in bstar.keys()):
                # this can happen for symbols in the first and last triplet
                continue

            fstar_states = sorted(fstar[state])
            bstar_states = sorted(bstar[state])
            nr_fstar_states = len(fstar_states)
            nr_bstar_states = len(bstar_states)
            if nr_fstar_states < 2 or nr_bstar_states < 2:
                # the 'transition rate' table must have be at least 2x2
                continue
            
            # build the table (using a matrix)
            mat = numpy.zeros((nr_bstar_states, nr_fstar_states))
            for i in xrange(nr_bstar_states):
                for j in xrange(nr_fstar_states):
                    bstate = bstar_states[i]
                    fstate = fstar_states[j]
                    try:
                        mat[i,j] = triplets[(bstate,state,fstate)]
                    except:
                        mat[i,j] = 0.
                
                # we have just filled the i-th row, normalize it to unit norm
                mat[i,:] /= numpy.linalg.norm(mat[i,:])
             
            # Evaluate the inner products
            matt = mat.transpose()
            for i in xrange(nr_bstar_states-1):
                bstate = bstar_states[i]
                row = mat[i,:]

                # We do all the inner products at once and then select just the smallest one
                # ! inners that have been evaluted already in past iterations are
                # ! skipped by setting them to 1.
                inners = row.dot(matt)
                inners[0:i] = 1.
                idx = numpy.argmin(inners)
                inner = inners[idx]
                if inner > 0.2:
                    #print "splitting: discarding split with to high inner for node %s : from=%s inner=%.2f" % (state, bstate, inner)
                    continue

                # add a new entry to the list of splits for the node with id 'sym'
                entry = (inner, bstate, bstar_states, fstar_states, row/row.sum())
                try:
                    #print "splitting: new split for node %s : %s" % (sym, entry)
                    splits[state].append(entry)
                except:
                    splits[state] = [entry]

        nr_splits = 0                        
        print "\nNode splitting:"
        for state, _splits in splits.iteritems():
            #print "  splits for node ", key
            nr_splits += len(_splits)
            print "  Nr. of splits for node %s: %d" % (state,len(_splits))
            #for entry in item:
                #print key, entry

        print "  Tot. nr. of splits %d" % nr_splits
        self._pending_splits = splits

    def print_info(self):
        print "transitional model dump: fill me"
        print "Symbols: "
        for i in xrange(len(self._symbols)):
            print "  ", i, self._symbols[i]
        print "States: "
        for i in xrange(len(self._init_states)):
            print "  ", i, self._init_states[i]



