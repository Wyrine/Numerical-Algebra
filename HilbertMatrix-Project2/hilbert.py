#!/usr/local/bin/python3

#Kirolos Shahat
#Math 472
#Project 2 -- Hilbert Matrix
#Due Date: March 8th, 2018
#    -Computing the inverse hilbert matrix for an 8x8 hilbert matrix
#        using direct, Choleski, and QR factorization and solving systems
#        to do so.


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

def chol(a):
    """
        Choleski factorization of a = LL^T
    """
    n = len(a)

    for k in range(n):
        a[k, k] = np.sqrt(a[k,k] - np.dot(a[k, 0:k], a[k, 0:k]))
        
        for i in range(k+1, n):
            a[i,k] = (a[i,k] - np.dot(a[i, 0:k], a[k, 0:k])) / a[k,k]
    for k in range(1, n):
        a[0:k, k] = 0.0
    return a

def cholSolve(L, b):
    """
        Solves a system Ax = LL^Tx = b using forward and back substitution
    """

    n = len(b)
    b[0] = b[0]/ L[0,0]
    #forward substitution
    for i in range(1, n):
        s = 0.0
        for j in range(0, i):
            s += L[i, j] * b[j]
        b[i] = (b[i] - s) / L[i, i]

    #back substitution
    b[n-1] = b[n-1] / L[n-1, n-1]
    for i in range(n-2, -1, -1):
        s = 0.0
        for j in range(i+1, n):
            s += L[j, i] * b[j]
        b[i] = (b[i] - s) / L[i][i]
    return b

def qrSolveSystem(q, r, b):
    """
        Solves a system Ax = QRx = b
    """
    n = len(b)
    #Rx = Q^T * b
    b = np.matmul(np.transpose(q), b)
    #back substitution
    b[n-1] = b[n-1] / r[n-1, n-1]
    for i in range(n-2, -1, -1):
        s = 0.0
        for j in range(i+1, n):
            s += r[i, j] * b[j]
        b[i] = (b[i] - s) / r[i][i]
    return b

def printMat(A):
    """
        Prints A in readable format
        Note: Done mostly for testing purposes
    """
    n = len(A)
    for i in range(n):
        print("\t", end="")
        for j in range(n):
            print(A[i, j], end="  ")
        print()

def CompareMethods(hilbert):
    """
        Compares approximations of the inverse hilbert matrix using choleski and QR
        to see which is more accurate
    """
    n = len(hilbert)
    invHilb = explicitInverse(hilbert)
    apprInv = np.zeros((n,n), dtype=np.float64)

    print("H inverse using equation 3:")
    printMat(invHilb)
    print("Condition number with infinity norm:", getCondition(hilbert, invHilb, n))

    #Using LL^T
    L = chol(np.array(hilbert))
    for i in range(n):
        e = np.zeros(n)
        e[i] = 1
        apprInv[:, i] = cholSolve(L, e)

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
