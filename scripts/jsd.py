#!/usr/bin/env python
from scipy.stats import entropy
from scipy.spatial import distance
from numpy.linalg import norm
import numpy as np

"""
scipy.stats.entropy(pk, qk)
If only probabilities pk are given, the entropy is calculated as S = -sum(pk * log(pk), axis=axis).
If qk is not None, then compute the Kullback-Leibler divergence S = sum(pk * log(pk / qk), axis=axis).

JS-divergence and KL-divergence
JSD is a method of measuring the similarity between two probability distributions
JSD is a symmetrized and smoothed version of the Kullback–Leibler divergence
https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence

About the logarithm base
log base 2: 0 <= JSD(P||Q) <= 1
log base e:  0 <= JSD(P||Q) <= ln(2)
"""

def JSD(P, Q, base=2):
    # p, q are numpy arrays
    _P = P / norm(P, ord=1) # ord = 1 -> max(sum(abs(p)))
    _Q = Q / norm(Q, ord=1)
    _M = 0.5 * (_P + _Q)
    return 0.5 * (entropy(_P, _M, base=base) + entropy(_Q, _M, base=base))

def tissue_specificity_score(p, q):
    p = np.asarray(p)
    q = np.asarray(q)
    # normalize
    p /= p.sum()
    q /= q.sum()
    return 1- distance.jensenshannon(p, q, base=2)
def maximal_tss(tss):
    tss_idx = np.argmax(tss)
    score = tss[tss_idx]
    return tss_idx, score


