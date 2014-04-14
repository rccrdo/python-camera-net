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

import numpy

import matplotlib
from matplotlib import rc
rc('font', family='sans-serif')
import matplotlib.pyplot as plt
import networkx as nx

from util import *


def transition_matrix_to_graph(A, statenames, threshold=0.):
    assert isinstance(A, numpy.ndarray)
    assert isinstance(statenames, list)
    N,N1 = A.shape
    assert N == N1
    assert N > 0
    assert len(statenames) == N
    assert threshold >= 0.

    # color definitions used for pretty plotting the graph
    clr_node = (0.95,0.95,0.95,0.25)
    clr_node_null = (0.85,0.85,0.85,0.5)
    clr_edge = (0.,0.,0.,0.25)
    clr_split = (0.95,0.15,0.15,0.25)

    graph = nx.DiGraph()
    graph.node_color = {}
    graph.edge_color = {}

    for i in xrange(N):
        for j in xrange(N):
            if A[i,j] <= threshold:
                continue

            def _add_graph_node(i):
                if i not in graph.nodes():
                    graph.add_node(i)
                    # setup nodes color
                    clr = clr_node
                    if '\'' in i:
                        clr = clr_split
                    elif '0' in i:
                        clr = clr_node_null
                    graph.node_color[i] = clr

            state_i = statenames[i]
            state_j = statenames[j]

            _add_graph_node(state_i)
            _add_graph_node(state_j)

            edge = (state_i,state_j)
            if edge not in graph.edges():
                graph.add_edge(*edge)
                graph.edge_color[edge] = clr_edge

    pos = {}
    for node in graph.nodes():
        node = repr(node).strip('\'')
        #print "node is ", node

        r = 5
        try:
            _node = int(node)
            assert '-' not in node
            #print node, " -> ", _node
            
            if _node == 0:
                p = (r,r)
            else:
                assert _node - 1 >= 0
                assert _node - 1 <= 5

                remap = [3,6,2,4,5,1]
                px = numpy.cos(2*(numpy.pi/6)*remap[_node-1])*r + r
                py = numpy.sin(2*(numpy.pi/6)*remap[_node-1])*r + r
                p = (px,py)
            pos[node] = p
        except:
            #print node, " -|"
            continue
            
        
    #print pos
    # update the graph layout
    USE_POS = 1
    if not USE_POS:
        pos = None
    graph.node_pos = nx.spring_layout(graph, iterations=100, scale=10, pos=pos)

    return graph


class HiddenMarkovModel():
    def __init__(self, A, B, statenames, trainset, validsets=[], pi=None, symnames=None, symbol2cam=None):
        assert isinstance(A, numpy.ndarray)
        assert isinstance(B, numpy.ndarray)
        assert pi is None or isinstance(pi, numpy.ndarray)
        assert isinstance(statenames, list)
        assert symnames is None or isinstance(symnames, list)
        assert symbol2cam is None or isinstance(symbol2cam, dict)

        M, N = B.shape
        M_ext = len(symnames)
        assert M_ext >= M

        if pi is None:
            # uniform initial distribution
            pi = numpy.ones(N)/N
        
        self._check_hmm(A, B, pi)

        # HMM matrices
        self._A = A
        self._B = B
        self._pi = pi
        
        self._init_N = N
        
        self._statenames = statenames
        self._symnames = symnames
        self._symbol2cam = symbol2cam
        
        def _check_dataset(dataset, M):
            assert isinstance(dataset, list)
            assert isinstance(M, int)
            assert M > 0
            assert len(dataset) > 0
            for i in dataset:
                assert i < M

        # sanity checks on all datasets
        _check_dataset(trainset, M)
        for dataset in validsets:
            _check_dataset(dataset, M_ext)

        self._trainset = trainset
        self._validsets = validsets        

        # define the arrays that will contain the perf. indices
        self._perf_prediction_train = []
        self._perf_prediction_valid = []
        for dataset in validsets:
            self._perf_prediction_valid.append([])

        self._perf_ratio_train = []
        self._perf_ratio_valid = []
        for dataset in validsets:
            self._perf_ratio_valid.append([])


        self._split_children = {}

        # plotting ...
        self._figure_perf = None
        self._figure_perf_ratio = None
        self._figure_graph = None

        self._bw_iter = 0

    def _check_hmm(self, A, B, pi):
        assert isinstance(A, numpy.ndarray)
        assert isinstance(B, numpy.ndarray)
        assert isinstance(pi, numpy.ndarray)

        N,N2 = A.shape
        assert N == N2
        M, N2 = B.shape
        assert N == N2
        assert M > 0
        N2, = pi.shape
        assert N == N2

        rowsum_A = A.sum(1)
        colsum_B = B.sum(0)
        rowsum_pi = pi.sum()
        
        if not numpy.isclose(rowsum_A, 1).all():
            print rowsum_A
            assert 0
        
        if not numpy.isclose(colsum_B, 1).all():
            print colsum_B
            assert 0
            
        if not numpy.isclose(rowsum_pi, 1):
            print rowsum_pi
            assert 0

    def _list_state_with_split_children(self, state):
        assert state >= 0 and state < self._init_N
        l = [state]
        try:
            l.extend(self._split_children[state])
        except:
            # 'state' has no children
            pass
        return l

    def split_state(self, state, bstate, bstar, fstar, rates):
        assert isinstance(state, int)
        assert state >= 0 and state < self._init_N
        assert isinstance(bstate, int)
        assert bstate >= 0 and bstate < self._init_N
        assert isinstance(bstar, list)
        assert isinstance(fstar, list)
        assert isinstance(rates, numpy.ndarray)
        assert numpy.isclose(rates.sum(), 1)
	assert rates.size == len(fstar)

        for _state in bstar:
            assert isinstance(_state, int)
            assert _state >= 0 and _state < self._init_N
            
        for _state in fstar:
            assert isinstance(_state, int)
            assert _state >= 0 and _state < self._init_N

        print "\nHMM splitting state \'%s\'" % self._statenames[state]
        print "   transition from \'%s\'" % self._statenames[bstate]
        print "   backward star   %s" % repr([self._statenames[i] for i in bstar])
        print "   forward star    %s" % repr([self._statenames[i] for i in fstar])
        print "   rates           %s" % ', '.join(['%.2f' % i for i in rates])

        A = numpy_matrix_grow(self._A)
        B = numpy_matrix_grow(self._B, 0, 1)
        M, N = B.shape
        N1, N2 = A.shape

	# update the list of children
        try:
            self._split_children[state].append(N-1)
        except:
            self._split_children[state] = [N-1]

        self._statenames.append(self._statenames[state]+'\''*len(self._split_children[state]))

        # update A
        epsilon = 0.05
        _PATCHING_METHOD = 2
        if _PATCHING_METHOD == 1:
            # patch backward star
            for _bstate in bstar:
                indexes = self._list_state_with_split_children(_bstate)
                if _bstate == bstate:
                    for __bstate in indexes:
                        A[__bstate,N-1] = max(epsilon, A[__bstate, state])
                else:
                    for __bstate in indexes:
                        A[__bstate,N-1] = epsilon
            # patch forward star
            for i in xrange(len(fstar)):
                _fstate = fstar[i]
                indexes = self._list_state_with_split_children(_fstate)
                #print "indexes for %d are %s" % (_fstate, repr(indexes))
                val = max(epsilon, rates[i])
                for __fstate in indexes:
                    #print " fstar: setting A[%d,%d] = %.3f" % (N-1, idx, val)
                    A[N-1, __fstate] = val
        else:
            assert _PATCHING_METHOD == 2
            # patch backward star
            for _bstate in bstar:
                val = epsilon
                if _bstate == bstate:
                    val = max(epsilon, A[_bstate, state])
                A[_bstate,N-1] = val
            # patch backward star
            for i in xrange(len(fstar)):
                _fstate = fstar[i]
                A[N-1, _fstate] = max(epsilon, rates[i])
        
        # normalize A
        A = numpy.diag(1/A.sum(1)).dot(A)

        # update and normalize B
        B[:,N-1] = B[:,state]
        B = B/B.sum(0)

        # update pi
        pi = self._pi.tolist()
        pi.append(1./(len(pi)+1))
        pi = numpy.array(pi)
        pi /= pi.sum()

        self._check_hmm(A, B, pi)
        self._A = A
        self._B = B
        self._pi = pi

        #if 0:
        #    print "matrices after splitting"
        #    print A
        #    print "rowsum A: ", A.sum(1)
        #    print B
        #    print "colsum B: ", B.sum(0)
        #    print "pi=", pi

    def baumwelch_until_converged(self, nu=0.001):
        assert isinstance(nu, float)
        assert nu > 0.0001
        assert nu <= 0.01
        
        while True:
            done = self.baumwelch(nu)
            if done:
                break

    def baumwelch(self, nu=0.01):
        assert isinstance(nu, float)
        assert nu > 0.0001
        assert nu <= 0.01
        
        self._bw_iter += 1
        print  "\nBaum-Welch iter. %d ..." % self._bw_iter
        
        # iterate BW until parameters converge within epsilon
        A = self._A.copy()
        B = self._B.copy()
        pi = self._pi.copy()

        self._do_baumwelch()

        converged = numpy.allclose(self._pi, pi, atol=nu)     \
                    and numpy.allclose(self._B, B, atol=nu)   \
                    and numpy.allclose(self._A, A, atol=nu)

        self._compute_perf_indices()
        
        return converged

    def _do_baumwelch(self):
        A = self._A
        B = self._B
        pi = self._pi
        data = self._trainset
        M, N = B.shape
        T = len(data)
        assert T > 0

	_DEBUG = False

        # Compute the foward variables
        # here alpha[t,i] is alpha_{t+1}(i+1) in the paper
        alpha = numpy.zeros((T,N))

        # init alpha[0,i], 0 <= i < N
        O_0_idx = data[0]
        alpha[0,:] = pi*B[O_0_idx,:]
        
        # from t=0 to t=T-2
        for t in xrange(T-1):
            O_tp1_idx = data[t+1]
            alpha[t+1,:] = (alpha[t,:].dot(A))*B[O_tp1_idx,:]
        
        # backward variables
        # here beta[t,i] is alpha_{t+1}(i+1) in the paper
        beta = numpy.zeros((T,N))
        
        # init beta[T-1,i], 0 <= i < N
        beta[T-1,:] = numpy.ones(N)
        
        # from t=T-2 to t=0
        for t in xrange(T-2, -1, -1):
            O_tp1_idx = data[t+1]
            for i in xrange(N):
                beta[t,i] = (A[i,:]*B[O_tp1_idx,:]).dot(beta[t+1,:])

        P_OLambda = alpha[T-1,:].sum()
        t = int(N/2)
        P_OLambda2 = alpha[t,:].dot(beta[t,:])
        assert P_OLambda >= 0
        assert P_OLambda <= 1
        if not numpy.isclose(P_OLambda, P_OLambda2):
            print "warning different values for P_OLambda", P_OLambda, P_OLambda2
            assert 0
        if not P_OLambda > 0:
            print "warning P_OLambda is zero !", P_OLambda
            assert 0
            
        # gamma
        gamma = (alpha*beta)/P_OLambda
        if _DEBUG:
            for t in xrange(0,T):
                if not numpy.isclose(gamma[t,:].sum(),1.):
                    assert 0
            if 1:
                r,c = gamma.shape
                print "gamma"
                for row in xrange(r):
                    triplet = data[row:row+3]
                    if triplet == [6,0,6]:
                        print "**********************"
                
                    print row, "\t: ", data[row], "\t", numpy_matrix_row_str(gamma[row,:])
                print "A"
                r,c = A.shape
                for row in xrange(r):
                    print row, "\t: ", data[row], "\t", numpy_matrix_row_str(A[row,:])


        # zeta
        # ! this would be a rather big 3D matrix
        def _zeta(t, i, j):
            assert t >= 0
            assert t < len(data) -1
            assert i >= 0
            assert i < N
            assert j >= 0
            assert j < N
            
            # skip the following computations if the transition i->j
            # is not ammissible
            a_ij = A[i,j]
            if a_ij == 0.:
                return 0.
            
            O_tp1_idx = data[t+1]
            num = alpha[t,i]*a_ij*B[O_tp1_idx, j]*beta[t+1, j]
            return num/P_OLambda

        if _DEBUG:
            for t in xrange(T-1):
                for i in xrange(N):
                    x = numpy.zeros(N)
                    for j in xrange(N):
                        x[j] = _zeta(t,i,j)
                    
                    #print t, i, x, x.sum(), gamma[t,i]
                    if not numpy.isclose(x.sum(), gamma[t,i]):
                        print t, i, x, x.sum(), gamma[t,i]
                        assert 0

        sum_0_Tm1_gamma = gamma.sum(0)
        sum_0_Tm2_gamma = sum_0_Tm1_gamma - gamma[-1,:]

        if _DEBUG:
            for x in sum_0_Tm1_gamma:
                assert x > 0.
            for x in sum_0_Tm2_gamma:
                # in verita' qualche stato potrebbe essere raggiunto solo all'ultimo
                # istante in generale ...
                assert x > 0.
            print "sum_0_Tm1_gamma[]=", numpy_matrix_row_str(sum_0_Tm1_gamma)
            print "sum_0_Tm2_gamma[]=", numpy_matrix_row_str(sum_0_Tm2_gamma)
            print "sum_0_Tm1_gamma.sum()", sum_0_Tm1_gamma.sum(), T
            assert numpy.isclose(sum_0_Tm1_gamma.sum(), T)
            assert numpy.isclose(sum_0_Tm2_gamma.sum(), T -1)
        
        sum_0_Tm2_zeta = numpy.zeros((N,N))
        for i in xrange(N):
            for j in xrange(N):
                val = 0.
                # from t=0 to t=T-2
                for t in xrange(T-1):
                   O_tp1_idx = data[t+1]
                   val += alpha[t,i]*A[i,j]*B[O_tp1_idx, j]*beta[t+1,j]

                sum_0_Tm2_zeta[i,j] = val/P_OLambda

        if _DEBUG:
            if not numpy.isclose(sum_0_Tm2_zeta.sum(), T-1):
                print "sum_0_Tm2_zeta.sum", sum_0_Tm2_zeta.sum()
                assert 0

        #print "sums"
        #print sum_0_Tm2_zeta
        #print sum_0_Tm2_gamma
        tmp = sum_0_Tm2_gamma.copy()
        #tmp[numpy.isclose(sum_0_Tm2_gamma, 0)] = 1.
        #print tmp        

        barpi = gamma[0,:]
        barA = (sum_0_Tm2_zeta.transpose()/tmp).transpose()
        barB = numpy.zeros((M,N))
        for i in xrange(M):
            for j in xrange(N):
                den = sum_0_Tm1_gamma[j]
                assert den > 0
                
                #if not den:
                #    print "waring sum_0_Tm1_gamma[%d] is zero" % j
                #    print "sum_0_Tm1_gamma[]=", sum_0_Tm1_gamma
                #    print j, M, N
                #    print "A=", A
                #    print "B=", B 
                #    print 'probs to', A[:,j]
                #    print 'probs from', A[j,:]
                #    assert 0
                num = 0.
                for t in xrange(T):
                    O_t_idx = data[t]
                    if O_t_idx == i:
                        num += gamma[t,j]

                #val = 0.
                #if den:
                #    val = num/den
                val = num/den
                barB[i,j] = val
    
        self._check_hmm(barA, barB, barpi)
        
        self._A = barA
        self._B = barB
        self._pi = barpi
        
    def __compute_perf_index(self, dataset):
        A = self._A
        B = self._B
        M, N = B.shape
        pi = numpy.ones(N)/N
        self._check_hmm(A, B, pi)
        
        T = len(dataset)
        assert T > 0

        # forward variables
        # alpha[t,i] = alpha_{t+1}(i+1) nell'articolo
        alpha = numpy.zeros((T,N))

        # init alpha[0,i], 0 <= i < N
        O_0_idx = dataset[0]
        alpha[0,:] = pi*B[O_0_idx,:]

        # from t=0 to t=T-2
        for t in xrange(0,T-1):
            O_tp1_idx = dataset[t+1]
            if O_tp1_idx < M:
                alpha[t+1,:] = (alpha[t,:].dot(A))*B[O_tp1_idx,:]
            else:
                # The symbol O(t+1) is not in the vocabulary of \lambda
                alpha[t+1,:] = 0

        P_OLambda = alpha[T-1,:].sum()

        perf = numpy.power(P_OLambda, 1./T)
        return perf, 0
        
    def __compute_perf_index2(self, dataset):
        A = self._A
        B = self._B
        M, N = B.shape
        pi = numpy.ones(N)/N
        self._check_hmm(A, B, pi)
        
        T = len(dataset)
        assert T > 0

        t = 0
        alpha_t = None
        hits = 0
        ENOTRANS = False
        while t < T-1:
            #print "\n\n compute_perf_index2 t=%d, T=%d, M=%d, N=%d" % (t, T, M, N)
            alpha_tp1 = None
            if alpha_t is None:
                #print "  alpha_t is None"
                O_t = dataset[t]
                #print "    O_t = %d.id = \"%s\"" % (O_t, self._symnames[O_t])

                if O_t < M:
                    alpha_t = pi*B[O_t,:]
                else:
                    #print "    setting ENOTRANS!"
                    ENOTRANS = True
                    alpha_t = numpy.zeros(N)

            assert alpha_t is not None

            O_tp1 = dataset[t+1]
            #print "    O_tp1 = %d.id = \"%s\"" % (O_tp1, self._symnames[O_tp1])
            if O_tp1 < M:
                alpha_tp1 = (alpha_t.dot(A))*B[O_tp1,:]
            else:
                # The symbol O(t+1) is not in the vocabulary of \lambda
                alpha_tp1 = None
            
            
            t +=1

            if alpha_tp1 is None:
                #print "    alpha_tp1 is none, continue..."
                alpha_t = None
                continue

            P_O_1_t = alpha_t.sum()
            if P_O_1_t == 0.:
                #print "    P_O_1_t is 0"
                #print "    setting ENOTRANS!"
                ENOTRANS = True
                alpha_t = None
                continue

            symprob_tp1 = ((alpha_t.dot(A))*B).sum(1)/P_O_1_t
            #print "alpha_t ", alpha_t, alpha_t.sum()
            #print "alpha_tp1", alpha_tp1, alpha_tp1.sum()

            #print "symprob_tp1", symprob_tp1, symprob_tp1.sum()
            #for i in range(0,M):
            #    print "    P[hatO_tp1 = %s] = %.3f" % (self._symnames[i], symprob_tp1[i])

            #assert numpy.isclose(symprob_tp1.sum(), 1)
            _sum = symprob_tp1.sum()
            if not numpy.isclose(_sum, 1):
                # This can happen due to rounding errors in alpha_tp1
                print "! warning, symprob_tp1 does not sum to 1" , _sum
                symprob_tp1 /=_sum
                alpha_t = None
                ENOTRANS = True
                continue

            alpha_t = alpha_tp1

            camprob = {}
            for i in range(0, M):
                prob = symprob_tp1[i]
                cams = self._symbol2cam[i]
                for cam in cams:
                    try:
                         camprob[cam] += prob
                    except:
                         camprob[cam] = prob

            hatcam = None
            hatcam_prob = None
            for camid in camprob.keys():
                prob = camprob[camid]
                if hatcam is None:
                    hatcam = camid
                    hatcam_prob = prob
                elif prob > hatcam_prob:                    
                    hatcam = camid
                    hatcam_prob = prob

            assert hatcam is not None
            #print "hatcam %d with probability %.3f" % (hatcam, hatcam_prob)
            if not (hatcam_prob >= 0. and hatcam_prob <= 1.001):
                print "hatcam, hatcam_prob", hatcam, hatcam_prob
            if str(hatcam) in self._symnames[O_tp1]:
                #print "    Hit ! on symbol %s" % self._symnames[O_tp1]
                hits += 1

        if ENOTRANS:
            P_O_1_T = 0
        else:
            assert alpha_t is not None
            P_O_1_T = alpha_t.sum()

        perf = numpy.power(P_O_1_T, 1./T)
        perfratio = float(hits)/(T-1)
        #print perfratio
        return perf, perfratio



    def _compute_perf_indices(self):
        print "Computing perf indices:"
        if 0:
            compute_perf = self.__compute_perf_index
        else:
            compute_perf = self.__compute_perf_index2
        
        # training set
        perf, perfratio = compute_perf(self._trainset)
        self._perf_prediction_train.append(perf)
        self._perf_ratio_train.append(perfratio)
        print "  training (T=%d)\t: %.2f, %.2f" % (len(self._trainset)-1, perf, perfratio)

        # validation sets
        for i in xrange(0,len(self._validsets)):
            dataset = self._validsets[i]
            perf, perfratio = compute_perf(dataset)
            self._perf_prediction_valid[i].append(perf)
            self._perf_ratio_valid[i].append(perfratio)
            print "  validation-%d (T=%d)\t: %.2f, %.2f" % (i, len(dataset)-1, perf, perfratio)

    def plot_perf(self, savepathbase=None, savetag=""):
           
        _TEXT_SIZE = 18
        markers = 'os*'
        colors = 'brg'
        
        if 1:
            if not self._figure_perf:
                self._figure_perf = plt.figure(figsize=(6,6))
                plt.ion()
                plt.show()        
            fig = self._figure_perf
        
            # focus and setup the this figure
            plt.figure(fig.number)
            fig.clf()
            fig.patch.set_facecolor((1,1,1))

            axis = plt.subplot(1,1,1)
            plt.grid('on')
            plt.xlim(0., len(self._perf_prediction_train)-1)
            plt.ylim(0., 1)

            # plot trainset perf
            plt.plot(self._perf_prediction_train, 'k')
            plt.hold(True)

            # plot validsets perf
            for i in xrange(0, len(self._validsets)):
                clr = ""
                if i < len(colors):
                    clr = colors[i]
                plt.plot(self._perf_prediction_valid[i], clr+'--')

            plt.legend(['Training', 'Validation'], loc='upper left')

            plt.legend(['Training', 'Validation'], loc='upper left')
            plt.xlabel('iteration step')
            plt.ylabel('$\eta$', fontsize=_TEXT_SIZE)

            fig.canvas.draw()
        
            if savepathbase:
                savepathbase += "perf-%s-frame%04d.png" % (savetag, self._bw_iter)
                fig.savefig(savepathbase)
            
        if 0:
            if not self._figure_perf_ratio:
                self._figure_perf_ratio = plt.figure(figsize=(6,6))
                plt.ion()
                plt.show()

            fig = self._figure_perf_ratio
        
            # focus and setup the this figure
            plt.figure(fig.number)
            fig.clf()
            fig.patch.set_facecolor((1,1,1))

            axis = plt.subplot(1,1,1)
            plt.grid('on')
            plt.xlim(0., len(self._perf_ratio_train)-1)
            plt.ylim(0.5, 1)

            # plot trainset perf
            plt.plot(self._perf_ratio_train, 'k')
            plt.hold(True)

            # plot validsets perf
            for i in xrange(0, len(self._validsets)):
                clr = ""
                if i < len(colors):
                    clr = colors[i]
                plt.plot(self._perf_ratio_valid[i], clr+'--')
            
            plt.legend(['Training', 'Validation'], loc='upper left')
            plt.xlabel('iteration step')
            plt.ylabel(r'$\rho$', fontsize=_TEXT_SIZE)

            fig.canvas.draw()

            if savepathbase:
                savepathbase += "perf2-%s-frame%04d.png" % (savetag, self._bw_iter)
                fig.savefig(savepathbase)
        
            

    def plot_transition_graph(self, threshold=0., savepathbase=None, savetag=""):
        return

        #if self._bw_iter < 15:
        #    return

        if not self._figure_graph:
            self._figure_graph = plt.figure(figsize=(6,6))
            plt.ion()
            plt.show()

        fig = self._figure_graph
        
        # focus and setup the this figure
        plt.figure(fig.number)
        fig.clf()
        fig.patch.set_facecolor((1,1,1))
	gs = matplotlib.gridspec.GridSpec(1, 1, left=0.05, right=0.95, top=0.95, bottom=0.05, wspace=0, hspace=0)
        axis = plt.subplot(gs[0,0], aspect='equal')

        graph = transition_matrix_to_graph(self._A,  self._statenames, threshold)
        xmin = None
        xmax = None
        for x,y in graph.node_pos.values():
            if xmin == None:
                xmin = x
                xmax = x
                ymax = y
                ymin = y
                continue
            xmin = min(xmin, x)
            xmax = max(xmax, x)
            ymin = min(ymin, y)
            ymax = max(ymax, y)

        #plt.axis('off')
        #plt.grid('on')
        plt.xlim(xmin -0.5, xmax + 0.5)
        plt.ylim(ymin -0.5, ymax + 0.5)

        nx.draw(graph, pos=graph.node_pos, node_size=1000,			\
                node_color=[graph.node_color[v] for v in graph], 		\
                edge_color=[graph.edge_color[e] for e in graph.edges()],	\
                linewidths=0.5, ax=axis)

        fig.canvas.draw()

        if savepathbase:
            savepathbase += "graph-%s-frame%04d.png" % (savetag, self._bw_iter)
            fig.savefig(savepathbase)

