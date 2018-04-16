#!/usr/local/bin/python3

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
				q = np.matmul(np.linalg.inv(B), np.array(q)/sig) 
				rq = raleighQuotient(A,q)
				if abs(rq - mu) <= eps:
						return it, rq, q
		
def main():
		A = np.array([[-149, -50, -154], [537, 180, 546], [-27, -9, -25]])
		x0 = np.array([1,1,1]).transpose()
		print(stabilizedPowerMethod(np.array(A), np.array(x0)))
		print(RQ_Iteration(np.array(A),np.array(x0)))
		return 0

if __name__ == "__main__":
		exit(main())
