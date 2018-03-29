#!/usr/local/bin/python3

import sys
import numpy as np

def jacobi(A, n, b, x, xReal, eps):
		k = 0
		x2 = np.ones(n, dtype=np.float64)
		while True:
				for i in range(n):
						pass	
				k += 1
				if infinityNorm(x2 - x) < eps:
						print("here")
						break
		return

def gaussSeidel(A, n, b, x, eps):
		return

def sor(A, n, b, x, w, eps):
		return

def infinityNorm(x):
		maxVal = -1
		for val in x:
				if abs(val) > maxVal:
						maxVal = abs(val)
		return maxVal

if __name__ == "__main__":
		A = np.array([-1, 2, -1], dtype=np.float64)
		nList = [25, 50, 100, 200]
		wList = [1.78486, 1.88402, 1.93968, 1.96922]
		eps = 10**(-6)
		for n, w in zip(nList, wList):
				h = 1/(n + 1)
				x0 = np.zeros(n, dtype=np.float64)
				xReal = np.zeros(n, dtype=np.float64)
				b = np.zeros(n, dtype=np.float64)
					
				#generating real solutions for this size matrix and the b vector
				for i in range(0, n):
						xReal[i] = (i+1)*h * (1. - (i+1) *h)
						b[i] = 2* (h**2)
				jacobi(A, n, b, x0, xReal, eps)
