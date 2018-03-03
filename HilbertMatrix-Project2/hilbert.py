#!/usr/local/bin/python3

import sys
import numpy as np

def buildHilbert(n):
    """ 
        builds an n x n hilbert matrix
    """
    hilbert = np.zeros((n,n), dtype = np.float64)
    for i in range(n):
        for j in range(i, n):
            hilbert[j, i] = hilbert[i, j] = 1/ (i + 1 + j + 1 - 1)
    return hilbert

def binomialCoeff(n, r):
    """
        Computes n choose r
    """
    fac = {}
    nFact, rFact, nMinusRFact = factorial(n, fac), factorial(r, fac), factorial(n-r, fac)
    return nFact/(rFact * nMinusRFact)

def factorial(x, fac):
    """
        computes x factorial
    """
    if x in fac:
        return fac[x]
    if x <= 1:
        return 1
    fac[x] = x * factorial(x - 1, fac)
    return fac[x]

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: ./conditionNumber.py hilbertDimensions")
        sys.exit(1)
    n = int(sys.argv[1])
    hilbert = buildHilbert(n)
