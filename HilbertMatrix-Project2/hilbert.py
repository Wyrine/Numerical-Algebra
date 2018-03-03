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
            hilbert[j, i] = hilbert[i, j] = 1/ (i + 1 + j)
    return hilbert

def explicitInverse(hilbert, n):
    hilbInv = np.array(hilbert)

    for i in range(n):
        for j in range(i, n):
            hilbInv[j,i] = hilbInv[i, j] = (-1)**(i+j+2) * (i+j+1) \
                    * binomialCoeff(n+i, n-j-1) * binomialCoeff(n+j, n-i-1)\
                    * binomialCoeff(i+j, i)**2
    return hilbInv

def binomialCoeff(n, r, fac = {}):
    """
        Computes n choose r
    """
    nFact, rFact, nMinusRFact = factorial(n, fac), factorial(r, fac), factorial(n-r, fac)
    return nFact/(rFact * nMinusRFact)

def getCondition(hilb, invHilb, n):
    """
        get condition number using the infinity norm of a hilbert
    """
    return getInfinityNorm(hilb, n) * getInfinityNorm(invHilb, n)
            
def getInfinityNorm(mat, n):
    maxRow = np.float64("-inf")

    for i in range(n):
        curSum = 0
        for j in range(n):
            curSum += abs(mat[i, j])
        if maxRow < curSum:
            maxRow = curSum
    return maxRow

def factorial(x, fac):
    """
        computes x factorial
        Note: any non-positive values return 1
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
    invHilb = explicitInverse(hilbert, n)
    print("Condition number with infinity norm:", getCondition(hilbert, invHilb, n))
