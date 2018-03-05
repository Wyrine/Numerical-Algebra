#!/usr/local/bin/python3

import sys
import numpy as np
from math import sqrt

def buildHilbert(n):
    """ 
        builds an n x n hilbert matrix
    """
    hilbert = np.zeros((n,n), dtype = np.float64)
    for i in range(n):
        for j in range(i, n):
            hilbert[j, i] = hilbert[i, j] = 1/ (i + 1 + j)
    return hilbert

def explicitInverse(hilbert):
    """
        Calculate the inverse of a hilbert matrix using the explicit equation
    """
    n = len(hilbert)
    hilbInv = np.array(hilbert)

    for i in range(n):
        for j in range(i, n):
            #because i and j are 0 indexed, the values in each term are changed
            #to accomodate
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

def choleski(A):
    n = len(A)
    l = np.zeros((n, n), dtype=np.float64)
    
    l[0, 0] = np.sqrt(A[0, 0])
    for j in range(1, n):
        l[j, 0] = A[j, 0] / l[0,0]


    for i in range (1, n-1):
        sub = 0
        for k in range(i-1 ):
            sub += (l[i, k] ** 2)
        l[i, i] = np.sqrt(A[i, i] - sub)
        
        for j in range(i+1, n):
            sub = 0
            for k in range(i-1):
                sub += (l[j, k] / l[i, k])
            l[j, i] = (A[j, i] - sub) / l[i, i]

    sub = 0
    for k in range(n-1):
        sub += (l[n-1, k] ** 2)
    l[n-1, n-1] = np.sqrt(A[n-1, n-1] - sub)
    return l

def chol(a):
    n = len(a)

    for k in range(n):
        a[k, k] = np.sqrt(a[k,k] - np.dot(a[k, 0:k], a[k, 0:k]))
        
        for i in range(k+1, n):
            a[i,k] = (a[i,k] - np.dot(a[i, 0:k], a[k, 0:k])) / a[k,k]
    for k in range(1, n):
        a[0:k, k] = 0.0
    return a

def cholSolveSystem(L, b):
    n = len(b)
  # Solution of [L]{y} = {b}
    for k in range(n):
        b[k] = (b[k] - np.dot(L[k,0:k],b[0:k]))/L[k,k]
  # Solution of [L_transpose]{x} = {y}
    for k in range(n-1,-1,-1):
        b[k] = (b[k] - np.dot(L[k+1:n,k],b[k+1:n]))/L[k,k]
    return b

def qrSolveSystem(q, r, b):
    n = len(b)
    #forward subst
    for k in range(n):
        b[k] = (b[k] - np.dot(q[k, 0:k], b[0:k])) / q[k,k]

    for k in range(n-1, -1, -1):
        b[k] = (b[k] - np.dot(r[k+1:n,k], b[k+1:n])) / r[k,k]

    return b

def printMat(A):
    n = len(A)
    for i in range(n):
        for j in range(n):
            print(A[i, j], end=" ")
        print("\n")

def CompareMethods(hilbert):
    n = len(hilbert)
    invHilb = explicitInverse(hilbert)
    apprInv = np.zeros((n,n), dtype=np.float64)

    print("H inverse using equation 3:")
    printMat(invHilb)
    print("Condition number with infinity norm:", getCondition(hilbert, invHilb, n))

    #Using LL^T
    L = np.transpose(chol(np.array(hilbert)))
    for i in range(n):
        e = np.zeros(n)
        e[i] = 1
        apprInv[:, i] = cholSolveSystem(L, e)

    print("H inverse using LLT:")
    printMat(apprInv)

    print("Infinity norm using choleski:", getInfinityNorm(invHilb-apprInv, n))

    #Using QR
    q, r = np.linalg.qr(hilbert)
    for i in range(n):
        e = np.zeros(n)
        e[i] = 1
        apprInv[:, i] = qrSolveSystem(q, r, e)
    print("H inverse using QR:")
    printMat(apprInv)
    print("Infinity norm using QR:", getInfinityNorm(invHilb-apprInv, n))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: ./conditionNumber.py hilbertDimensions")
        sys.exit(1)
    CompareMethods(buildHilbert(int(sys.argv[1])))
