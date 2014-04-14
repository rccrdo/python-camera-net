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

# link filtering with an hypothesis test


import numpy
import scipy, scipy.stats


class HypoTest():
    def __init__(self, bargamma, alpha_0):
        assert isinstance(bargamma, float)
        assert isinstance(alpha_0, float)
        assert bargamma > 0.
        assert bargamma < 1.
        assert alpha_0 > 0.
        assert alpha_0 < 1.
        
        self._bargamma = bargamma
        self._alpha_0 = alpha_0
        self._thresholds = {}

    def _compute_thresholds(self, n):
	assert isinstance(n, int)
	assert n > 0

	bargamma = self._bargamma
        alpha_0 = self._alpha_0

        #print "_compute_thresholds(n=%d, bargamma=%.3f, alpha_0=%.3f)" % (n, bargamma, alpha_0)

        # compute the critical s and the remainder probability alpha
        for s in xrange(n+1):
            #print "s/n, %d,%d" % (s, n)
            prob1 = 1. - scipy.stats.binom.cdf(s-1, n, bargamma)
            prob2 = 1. - scipy.stats.binom.cdf(s, n, bargamma)
            #print "s=%d ?: prob1(>=%d) >= alpha0 >= prob2(>=%d)  %.2f >= %.2f >= %.2f" % (s, s, s+1, prob1, alpha_0, prob2)
            if prob1 >= alpha_0 and alpha_0 >= prob2:
                bars = s
                rem = alpha_0 - prob2
                
                assert rem >= 0.
                Ps = prob1-prob2
                assert Ps >= rem
                
                #print "rem=%.5f, Ps=%.5f" % (rem, Ps)
                #alpha = alpha_0 - prob2
                alpha = (prob1 - alpha_0)/Ps
                break
                


        assert bars >= 0
        assert bars <= n
        assert alpha >= 0
        assert alpha <= 1
        #print "hypo n, bars, alpha", n, bars, alpha
        self._thresholds[n] = (bars, alpha)
        return bars, alpha

    def test(self, nr_trials, nr_successes):
        _nr_trials = int(nr_trials)
        _nr_successes = int(nr_successes)
        assert _nr_trials > 0
        assert _nr_successes >= 0
        assert _nr_successes <= _nr_trials
        
        try:
            s, alpha = self._thresholds[_nr_trials]
        except:
            s, alpha = self._compute_thresholds(_nr_trials)

        if _nr_successes < s:
            # choose 0
            ret = 0
        elif _nr_successes > s:
            # choose 1
            ret = 1
        else:
            # if we get here then _nr_successes == s
            u = numpy.random.uniform()
            ret = int(u <= alpha)
            
        #print "test n=%d, s=%d choose=%d" % (nr_trials,nr_successes, ret)
        return ret
        
    def plot_power(self, trials, gammas):
        import matplotlib
        from matplotlib import rc
        rc('font', family='sans-serif')
        import matplotlib.pyplot as plt

        power = numpy.zeros((len(gammas), len(trials)))
        for i in range(len(trials)):
            n = trials[i]
            s, alpha = self._compute_thresholds(n)
            
            #print "test n=%d: s=%d, alpha=%.5f" % (n, s, alpha)
            
            for j in range(len(gammas)):
                gamma = gammas[j]
                prob = scipy.stats.binom.cdf(s-1, n, gamma)
                prob2 = scipy.stats.binom.cdf(s, n, gamma)
                
                power[j,i] = 1. - (prob + alpha*(prob2-prob))
                assert power[j,i] >= 0. 
                assert power[j,i] <= 1.
                #print "  power gamma=%.2f: %.5f" % (gamma, power[j,i])
                #print "  prob, prob2 = "
                
                prob = scipy.stats.binom.cdf(s-1, n, self._bargamma)
                prob2 = scipy.stats.binom.cdf(s, n, self._bargamma)
                size = 1. - (prob + alpha*(prob2-prob))
                #assert numpy.isclose(size, self._alpha_0)
                #print "  size  %.5f" % (size)
                
        X, Y = numpy.meshgrid(trials,gammas)
        fig = plt.figure(figsize=(6,6))
        plt.ion()
        plt.show()        
        _FONTSIZE = 22
        CS = plt.contour(X, Y, power, colors='k')
        plt.clabel(CS, inline=1, fontsize=12)
        plt.ylabel('$\gamma$', fontsize=_FONTSIZE)
        plt.xlabel('$n$', fontsize=_FONTSIZE)
        plt.title(r'$\bar \gamma=0.1$, $\alpha_0=0.05$', fontsize=_FONTSIZE)
        
        plt.draw()
        fig.savefig('./hypotest-power')



if __name__ == "__main__":
    nmax = 1000
    nstep = 20
    bargamma = 0.1
    alpha0 = 0.05
    trials = range(1,101)
    #trials = [1]
    
    gammas = numpy.arange(0,1.1,0.1)

    test = HypoTest(bargamma, alpha0)
    test.plot_power(trials, gammas)
    
    raw_input("Press a key...")
