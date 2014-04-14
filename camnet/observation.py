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
import collections
import copy
import math
import numpy

from math2D import numeric_real

from util import *
from hypo import HypoTest


class ObservationAccu():
    def __init__(self):
        self._sequence = []     # camera network output
        
    def __repr__(self):
        return __repr__(self._sequence)

    def accumulate(self, sym):
        assert isinstance(sym, str)
        assert len(sym)

        # defensive format checks
        assert '~' not in sym
        if '0' in sym:
            assert sym == '0'

        if len(self._sequence) == 0:
            self._sequence.append(sym)
            return

        last_sym = self._sequence[-1]
        if last_sym == sym:
            # discard successive coincident symbols
            return
        
        self._sequence.append(sym)

    def get_sequence(self):
        return self._sequence

    def save(self, savepath):
        # Save observation data to file
        print "Saving observations to ", savepath
        if os.path.exists(savepath):
            answer = raw_input("An experiment with the same name already exists. Do yout want to overwrite it ? [y/N] ")
            if answer.lower() == 'y':
                os.remove(savepath)
            else:
                return False

        file(savepath, 'w').write(repr(self._sequence))
        return True




class ObservationSequence():
    def __init__(self, loadpath=None, max_trans=0, name="NoName", vnodes=False, \
    		admissible_statespace=None, admissible_transitions=None):
        #assert sequence is None or isinstance(sequence, list)
        assert loadpath is None or isinstance(loadpath, str)
        assert isinstance(max_trans, int)
        assert max_trans >= 0
        assert isinstance(name, str)
        assert isinstance(vnodes, bool)
        assert not ((admissible_statespace is None) ^ (admissible_transitions is None))

        self._max_trans = max_trans	# max nr. of transition to accumulate (0 for no limit)
        self._name = name		# the sequence name
        self._vnodes = vnodes		# True if vnodes are used
        self._filterseq = None# filterseq     # this sequence must fit into the model of filterseq
        self._nr_skipped = 0		# nr of symbols skipped while processing input sequence
        self._matrix = None		# transition matrix
        self._idx = []			# state-id to row/col index in self._matrix
        self._orig_outseq = []		# un-processed output-sequence
        self._orig_stateseq = []	# un-processed state-sequence
        self._outseq = []		# processed output-sequence
        self._stateseq = []		# processed state-sequence
        self._accu_done = False         # true if max-trans transitions have been accumulated

        self._triplets = {}		# counts all triplets in self._sequence
        #self._last_triplet = None	# the last triplet inserted in self._triplets

        self._admissible_statespace = admissible_statespace
        self._admissible_transitions = admissible_transitions

        #if sequence and loadpath:
        #    print "warning, ObservationSequence.__init__(sequence=[...], loadpath=%s), " \
        #    	  "both sequence and loadpath set! discarding loadpath" % (loadpath)
        #    loadpath = None

        sequence = []
        if loadpath:
            if not os.path.isfile(loadpath):
                print "warning, ObservationSequence(loadpath=%s), no data file at path\n" % loadpath
            else:
                print "Loading ObservationSequence at \'%s\'\n" % loadpath

                if name == "NoName":
                    self._name = os.path.basename(loadpath)

                data = open(loadpath, 'r').read()
                sequence = eval(data)			# don't do this at home
                assert isinstance(sequence, list)

        if sequence:
            def _reduce(seq):
                _seq = []
                for entry in seq:
                    entry = entry.replace('.','-')
                    if not len(_seq):
                        _seq.append(entry)
                    else:
                        # discard successive coincident symbols
                        last_entry = _seq[-1]
                        if last_entry != entry:
                            _seq.append(entry)
                return _seq

            def _mod(seq):
                if not self._vnodes:
                    return seq
                _seq = []
                for i in xrange(len(seq)):
                    sym  = seq[i]
                    if sym == '0' and i >= 1 and i <= len(seq)-1:
                        prev = seq[i-1]
                        next = seq[i+1]
                        null = self._make_null_state(prev, next)
                        _seq.append(null)
                    else:
                        _seq.append(sym)
                return _seq

            self._orig_outseq = _reduce(sequence)
            if len(self._orig_outseq) and self._orig_outseq[0] == '0':
                self._orig_outseq.pop(0)
            if len(self._orig_outseq) and self._orig_outseq[-1] == '0':
                self._orig_outseq.pop(-1)

            self._orig_stateseq  = _mod(self._orig_outseq)
            assert len(self._orig_outseq) == len(self._orig_stateseq)
            for i in xrange(len(self._orig_outseq)):
                if self._accu_done:
                    break

                sym = self._orig_outseq[i]
                state = self._orig_stateseq[i]
                self._accumulate(sym, state)

            assert len(self._outseq) == len(self._stateseq)
            if len(self._outseq):
                assert len(self._outseq) >= 2
                assert self._outseq[0] != '0'
                if self._vnodes:
                    assert self._outseq[-1] != '0'
                #prev = self._sequence[-2]
                #last = '0'
                #self._untrack_last_triplet()
                #self._untrack_transition(prev, '0')
                #self._sequence.pop()
                #self._sequence_mod.pop()

        #if hypotest:
        #    self._do_hypothesis_test(hypotest)

        self._simplify()

    def orig_outseq_len(self):
        return len(self._orig_outseq)

    def orig_stateseq_len(self):
        return len(self._orig_stateseq)

    def outseq_len(self):
        return len(self._outseq)

    def stateseq_len(self):
        return len(self._stateseq)

    def orig_outseq(self):
        return self._orig_outseq

    def orig_stateseq(self):
        return self._orig_stateseq

    def outseq(self):
        return self._outseq

    def stateseq(self):
        return self._stateseq

    def orig_states(self):
        return tuple(set(self._orig_stateseq))

    def nr_orig_states(self):
        return len(self.orig_states())

    def orig_symbols(self):
        return tuple(set(self._orig_outseq))

    def nr_orig_symbols(self):
        return len(self.orig_symbols())

    def states(self):
        return tuple(set(self._stateseq))

    def nr_states(self):
        return len(self.states())

    def symbols(self):
        return tuple(set(self._outseq))

    def nr_symbols(self):
        return len(self.symbols())

    def dump(self):
        print "Transition matrix:"
        for i in xrange(self.nr_states()):
            print i, self._idx[i], numpy_matrix_row_str(self._matrix[i,:])

    def with_vnodes(self):
        return self._vnodes

    def state_triplets_dict(self):
        return self._triplets

    def state_fstar_dict(self):
        N, N1 = self._matrix.shape
        assert N == N1
        assert N == len(self._idx)
        fstar = {}
        for i in xrange(N):
            state_i = self._idx[i]
            for j in xrange(N):
                if self._matrix[i,j] > 0:
                    state_j = self._idx[j]
                    try:
                        fstar[state_i].append(state_j)
                    except:
                        fstar[state_i] = [state_j]
        return fstar
        
    def state_bstar_dict(self):
        N, N1 = self._matrix.shape
        assert N == N1
        assert N == len(self._idx)
        bstar = {}
        for j in xrange(N):
            state_j = self._idx[j]
            for i in xrange(N):
                if self._matrix[i,j] > 0:
                    state_i = self._idx[i]
                    try:
                        bstar[state_j].append(state_i)
                    except:
                        bstar[state_j] = [state_i]
        return bstar                

    def _make_null_state(self, src, dest):
        null = sorted([src, dest])
        null.insert(1, '0')
        return '~'.join(null)

    def _get_idx(self, state, grow=False):
        assert isinstance(state, str)
        assert len(state)
        assert isinstance(grow, bool)

        try:
            idx = self._idx.index(state)
        except:
            if not grow:
                return None
            idx = len(self._idx)
            self._idx.append(state)
            if self._matrix == None:
                self._matrix = numpy.zeros((1,1))
            else:
                # extend the transition matrix with one row and one column
                self._matrix = numpy_matrix_grow(self._matrix)
        return idx

    def _track_triplet(self, src, over, dest):
        assert isinstance(src, str)
        assert isinstance(over, str)
        assert isinstance(dest, str)
        assert src != over
        assert over != dest

        # increase count
        _tuple = (src, over, dest)
        try:
            self._triplets[_tuple] += 1
        except:
            self._triplets[_tuple] = 1

        #self._last_triplet = _tuple
        return self._triplets[_tuple]
        
    #def _untrack_triplet(self, src, over, dest):
    #    assert isinstance(src, str)
    #    assert isinstance(over, str)
    #    assert isinstance(dest, str)
    #    assert src != over
    #    assert over != dest
    #
    #    # decrease count
    #    _tuple = (src, over, dest)
    #    try:
    #        self._triplets[_tuple] -= 1
    #    except:
    #        # untracking an unknown triplet should never happen
    #        assert 0
    #
    #    val = self._triplets[_tuple]
    #    assert val >= 0
    #    return val
        
    #def _untrack_last_triplet(self):
    #    assert self._last_triplet is not None
    #    src, over, dest = self._last_triplet
    #    return self._untrack_triplet(src, over, dest)

    def _track_transition(self, src, dest):
        assert isinstance(src, str)
        assert isinstance(dest, str)
        src_idx = self._get_idx(src, grow=True)
        dest_idx = self._get_idx(dest, grow=True)
        self._matrix[src_idx][dest_idx] += 1
        return self._matrix[src_idx][dest_idx]

    def _untrack_transition(self, src, dest):
        assert isinstance(src, str)
        assert isinstance(dest, str)
        src_idx = self._get_idx(src)
        dest_idx = self._get_idx(dest)
        # untracking a transition from/to unknown states should never happen
        assert src_idx is not None
        assert dest_idx is not None
        self._matrix[src_idx][dest_idx] -= 1
        val = self._matrix[src_idx][dest_idx]
        assert val >= 0
        return val

    def count_transitions_from_to(self, src, dest):
        assert src is not None
        assert dest is not None
        src_idx = self._get_idx(src)
        dest_idx = self._get_idx(dest)
        if (src_idx is None) or (dest_idx is None):
            print "warning, ObservationSequence.count_transitions_from_to(%s,%s)," \
                  " src or dest node does not exist!" % (src, dest)
            return 0
        return self._matrix[src_idx, dest_idx]

    def count_transitions_from(self, nodeid):
        idx = self._get_idx(nodeid)
        if idx is None:
            print "warning, ObservationSequence.count_transitions_from(%s)," \
                  " node does not exist!" % (nodeid)
            return 0
        return self._matrix[idx, :].sum()

    def count_transitions_to(self, nodeid):
        idx = self._get_idx(nodeid)
        if idx is None:
            print "warning, ObservationSequence.count_transitions_to(%s)," \
                  " node does not exist!" % (nodeid)
            return 0
        return self._matrix[:, idx].sum()
        
    #def has_symbol(self, sym):
    #    idx = self._get_idx(sym)
    #    if idx == None:
    #        return False
    #    count = self.count_transitions_to(sym)
    #    assert count > 0 
    #    return True

    #def has_transition(self, src, dest):
    #    src_idx = self._get_idx(src)
    #    dest_idx = self._get_idx(dest)
    #    if (src_idx is None) or (dest_idx is None):
    #        return False
    #    count = self._matrix[src_idx, dest_idx]
    #    return count > 0

    def _accumulate_simple(self, sym):
        #print "accumulate_simple", self._name, sym, '5' in self._outseq

        if self._admissible_statespace:
            if sym not in self._admissible_statespace:
                print "warning, ObservationSequence \'%s\': filtering out symbol %s" % (self._name, sym)
                return

        if len(self._outseq) == 0:
            if sym == '0':
                print "warning, ObservationSequence \'%s\': filtering out first '0' symbol" % self._name
                return
            self._outseq.append(sym)
            self._stateseq.append(sym)
            return

        last_sym = self._outseq[-1]
        assert last_sym == self._stateseq[-1]
        if last_sym == sym:
            # discard successive coincident symbols
            return
        
        if self._admissible_transitions:
            if (last_sym, sym) not in self._admissible_transitions:
                print "warning, ObservationSequence \'%s\': filtering out transition %s -> %s" % (self._name, last_sym, sym)
                return

        self._outseq.append(sym)
        self._stateseq.append(sym)
        self._track_transition(last_sym, sym)
        
        if len(self._outseq) >= 3:
            triplet = self._outseq[-3:]
      	    src = triplet[0]
	    over = triplet[1]
	    dest = triplet[2]
	    self._track_triplet(src, over, dest)

    def _accumulate_vnodes(self, sym, state):
        #if self._filterseq:
        #    assert 0
        #    # implement me ...

        if self._max_trans:
            if (len(self._outseq) +1 > self._max_trans) and sym == '0':
                self._accu_done = True
                print "warning, ObservationSequence \'%s\': filtering out trailing '0' symbol" % self._name
                return
        else:
            if (len(self._outseq)  == len(self._orig_outseq) -1) and sym == '0':
                self._accu_done = True
                print "warning, ObservationSequence \'%s\': filtering out trailing '0' symbol" % self._name
                return
                
        if self._admissible_statespace:
            if state not in self._admissible_statespace:
                print "warning, ObservationSequence \'%s\': filtering out state %s" % (self._name, state)
                return

        if len(self._outseq) == 0:
            if sym == '0':
                print "warning, ObservationSequence \'%s\': filtering out first '0' symbol" % self._name
                return
            self._outseq.append(sym)
            self._stateseq.append(state)
            return

        last_sym = self._outseq[-1]
        last_state = self._stateseq[-1]
        if '~0~' in last_state:
            assert last_sym == '0'
        else:
            assert last_sym == last_state
        if last_sym == sym:
            # discard successive coincident symbols
            print sym, last_sym, last_state, state
            if not self._admissible_statespace:
                # can stil happen when not filtering
                assert '~0~' not in last_state 
            return

        if self._admissible_transitions:
            if (last_state, state) not in self._admissible_transitions:
                print "warning, ObservationSequence \'%s\': filtering out transition %s -> %s" % (self._name, last_state, state)
                return

        self._outseq.append(sym)
        self._stateseq.append(state)

        self._track_transition(last_state, state)

        if len(self._stateseq) >= 3:
            triplet = self._stateseq[-3:]
      	    src = triplet[0]
	    over = triplet[1]
	    dest = triplet[2]
	    self._track_triplet(src, over, dest)


    def _accumulate(self, sym, state):
        assert isinstance(sym, str)
        assert isinstance(state, str)
        assert len(sym)
        assert len(state)
        assert len(self._outseq) == len(self._stateseq)

        if self._max_trans and len(self._outseq) > self._max_trans:
                # accu_done should hove been set already
                assert 0

        # defensive format checks
        assert '~' not in sym
        if '0' in sym:
            assert sym == '0'

        if not self._vnodes:
            assert sym == state
            self._accumulate_simple(sym)
        else:
            if sym == '0':
                assert state == '0' or '~0~' in state
            self._accumulate_vnodes(sym, state)

        if self._max_trans and len(self._outseq) > self._max_trans:
            self._accu_done = True

    def _simplify(self):
        remove = []
        for i in xrange(len(self._idx)):
            state = self._idx[i]
            nr_from = self.count_transitions_from(state)
            nr_to = self.count_transitions_to(state)
            #print "simplify: state %s, nr_from %d, nr_to %d" % (state, nr_from, nr_to)
            if nr_from == 0 and nr_to ==0:
                remove.append(i)

        removed = 0
        for i in sorted(remove):
            idx = i-removed
            print "ObservationSequence %s: cleaning-up state %s" % (self._name, self._idx[idx])
            self._idx.pop(idx)

            self._matrix = numpy.delete(self._matrix, idx, 0)
            self._matrix = numpy.delete(self._matrix, idx, 1)

            removed += 1

    def __admissible_statespace_transitions(self, idx, A):
        assert isinstance(A, numpy.ndarray)
        N,N2 = A.shape
        assert N > 0
        assert N == N2
        assert len(idx) == N

        statespace = []
        transitions = []
        for i in xrange(N):
            for j in xrange(N):
                if A[i,j] > 0:
                    state_i = idx[i]
                    state_j = idx[j]
                    if 0:#'~0~' in state_i:
                        state_from, state_to = state_i.split('~0~')
                        _from = idx.index(state_from)
                        _to = idx.index(state_to)
                        assert _from is not None
                        assert _to is not None
                        transitions.append((_from,i))
                        transitions.append((i,_to))
                        statespace.append(_from)
                        statespace.append(i)
                        statespace.append(_to)
                    elif 0:#'~0~' in state_j:
                        state_from, state_to = state_j.split('~0~')
                        _from = idx.index(state_from)
                        _to = idx.index(state_to)
                        assert _from is not None
                        assert _to is not None
                        transitions.append((_from,j))
                        transitions.append((j,_to))
                        statespace.append(_from)
                        statespace.append(j)
                        statespace.append(_to)
                    else:
                        transitions.append((i,j))
                        statespace.append(i)
                        statespace.append(j)

        _statespace = sorted(tuple(set(statespace)))
        _transitions = sorted(tuple(set(transitions)))
        statespace = []
        transitions = []
        for state in _statespace:
            statespace.append(idx[state])
        for src,dest in _transitions:
            transitions.append((idx[src], idx[dest]))

        return statespace, transitions

    def admissible_statespace_transitions(self):
        return self.__admissible_statespace_transitions(self._idx, self._matrix)
        
    def hypotest_admissible_statespace_transitions(self, htest):
        assert isinstance(htest, HypoTest)

        matrix = copy.copy(self._matrix)
        N,N2 = matrix.shape
        assert N > 0
        assert N == N2
        assert len(self._idx) == N

        for src in xrange(N):
            state_src = self._idx[src]
            if not self._vnodes and state_src == '0':
                # when using the basic coverage model don't test transitions
                # starting form '0': there are many ramifications with relatively
                # few transitions each
                continue
            for dest in xrange(N):
                 state_dest = self._idx[dest]
                
                 nr_from = self.count_transitions_from(state_src)
                 nr_trans = self.count_transitions_from_to(state_src, state_dest)
                 if nr_trans > 0:
                     choose = htest.test(nr_from, nr_trans)
                     if choose == 0:
                         if '~0~' in state_src:
                             state_from, state_to = state_src.split('~0~')
                             _from = self._get_idx(state_from)
                             _to = self._get_idx(state_to)
                             assert _from is not None
                             assert _to is not None
                             matrix[_from,src] = 0
                             matrix[src,_to] = 0
                             matrix[src, _from] = 0
                             matrix[_to, src] = 0
                         elif '~0~' in state_dest:
                             state_from, state_to = state_dest.split('~0~')
                             _from = self._get_idx(state_from)
                             _to = self._get_idx(state_to)
                             assert _from is not None
                             assert _to is not None
                             matrix[_from,dest] = 0
                             matrix[dest,_to] = 0
                             matrix[_to, dest] = 0
                             matrix[dest,_from] = 0
                         else:                            
                             #print "_do_hypothesis_test: transition %s -> %s, nr_from %d, nr_trans %d to be removed" % (state_src, state_dest, nr_from, nr_trans)
                            matrix[src,dest] = 0

        return self.__admissible_statespace_transitions(self._idx, matrix)
    

