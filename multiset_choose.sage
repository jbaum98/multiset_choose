#!/usr/bin/env sage

import sys
from sage.all import *

Sym = SymmetricFunctions(QQ)
p = Sym.powersum()

def sub_cycle_index(Zout, Zin):
    """Substitutes Zin into Zout to produce Zout evaluated at Zin.

    This is accomplished by replacing p[i] in Zout with Zin, but with
    every p[j] in Zin replaced with p[i*j].
    """

    return p.sum(c * p.prod(Zin.frobenius(i) for i in mu) for mu, c in Zout)

def multiset_cycle_index(ms):
    """The cycle index of the given multiset, given by the formula

    .. math:: \prod_{\{k\}}\left( Z(S_{\sigma_k}; Z(S_k))\right)

    where :math:`\{k\}` are the elements of the multiset and
    :math:`\sigma_k` is the multiplicity of the element :math:`k`.
    """

    Z = lambda i: SymmetricGroup(i).cycle_index()
    return p.prod(sub_cycle_index(Z(sk), Z(k)) for k, sk in ms.items())

def list_to_multiset(els):
    """Converts a list of elements representing a multiset to a dictionary
    where the keys are the elements of the multiset and the values are
    the multiplicities.
    """

    ms = {}
    for x in els:
        ms[x] = ms.get(x,0) + 1
    return ms

def mset_choose(s, d):
    """Compute the "multiset coefficient" :math:`\binom{s}{d}`."""

    A = PolynomialRing(QQ, len(s), 'A').gens()
    mono = prod(a^i for a, i in zip(A, s))
    Z = multiset_cycle_index(list_to_multiset(d))
    return Z.expand(len(A), A).coefficient(mono)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: %s 's_1, s_2, ..' 'd_1, s_2, ..'" % sys.argv[0])
        print("Outputs the number of ways the multiset s can be partitioned into multisets of sizes d_i.")
        sys.exit(1)

    s = map(int, sys.argv[1].split(' '))
    d = map(int, sys.argv[2].split(' '))

    if sum(s) != sum(d):
        print("The sum of the elements of s must equal the sum of the elements of d")
        sys.exit(1)

    print(mset_choose(s, d))
