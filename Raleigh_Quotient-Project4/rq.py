#!/usr/local/bin/python3

# Author: Kirolos Shahat
# Due Date: 4/19/17
# Course: Math 472
# Description: Comparison of Stabilized Power Method vs. Raleigh Quotient Iterations
#								To Approximate Eigenvalues and Vectors


# Output of running:
#Stabilized Power Method:
#	K = 28 , Lambda = 2.9999805022117965 , eigenVector = [-0.14286309  1.         -0.18366766]
#Raleigh Quotient Iteration:
#	K = 13 , Lambda = 3.0000000000011577 , eigenVector = [ 0.14285714 -1.          0.18367347]
# Therefore, Raleigh Quotient Iteration converges to the eigenpair faster than the Power Method 

from sys import exit
import numpy as np

def infNorm(A):
		maxRow = np.float64("-inf")
		for a in A:
				if abs(a) > maxRow:
						maxRow = abs(a)
		return maxRow

def raleighQuotient(A, q):
		qT = np.transpose(q)
		return np.matmul(qT, np.matmul(A, q))/ np.matmul(qT, q)

def stabilizedPowerMethod(A, q, eps= 1e-5):
		it = 0
		prevRQ = 0
		while True:
				it += 1
				sigNext = infNorm(np.matmul(A, q))
				q = np.array(np.matmul(A, q)) / sigNext
				rq = raleighQuotient(A, q) 
				if abs(rq - prevRQ) <= eps:
						return it, rq, q
				prevRQ = rq

def RQ_Iteration(A, q, eps=1e-5):
		it = 0
		while True:
				it += 1
				mu = raleighQuotient(A, q)
				sig = infNorm(np.matmul(A,q))
				B = np.array(mu * np.identity(3) - A)
				q = np.matmul(np.linalg.inv(B), np.array(q))/sig
				rq = raleighQuotient(A,q)
				if abs(rq - mu) <= eps:
						return it, rq, q/infNorm(q)
		
def main():
		A = np.array([[-149, -50, -154], [537, 180, 546], [-27, -9, -25]])
		x0 = np.array([1,1,1]).transpose()
		spm = stabilizedPowerMethod(np.array(A), np.array(x0))
		rqi = RQ_Iteration(A,x0)
		print("Stabilized Power Method:\n\tK =", spm[0], ", Lambda =", spm[1], ", eigenVector =", spm[2])
		print("Raleigh Quotient Iteration:\n\tK =", rqi[0], ", Lambda =", rqi[1], ", eigenVector =", rqi[2])
		return 0

if __name__ == "__main__":
		exit(main())
