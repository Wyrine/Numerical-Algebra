#!/usr/local/bin/python3

# Kirolos Shahat
# Computing Project 3
# Math 472
# March 29, 2018

import sys
import numpy as np
from math import exp, log

def jacobi(n, b, xReal, eps):
		x = np.zeros(n, dtype=np.float64)
		k, x2 = 0, np.zeros(n, dtype=np.float64)
		while True:
				x2[0] = 1/2*(x[1] + b[0])
				for i in range(1, n-1):
						x2[i] = 1/2*(x[i-1] + x[i+1] + b[i])
				x2[n-1] = 1/2*(x[n-2] + b[n-1])
				k += 1
				x = x2
				infNorm = infinityNorm(x-xReal)
				if infNorm < eps:
						print("\tJacobi:")
						print("\t\tK:", k)
						print("\t\tSpectral Radius:", exp(1/k * log(infNorm/ infinityNorm(xReal), 2)))
						return

def gaussSeidel(n, b, xReal, eps):
		k=0
		x = np.zeros(n, dtype=np.float64)
		while True:
				x[0] = 1/2*(x[1] + b[0])
				for i in range(1, n-1):
						x[i] = 1/2*(x[i-1] + x[i+1] + b[i])
				x[n-1] = 1/2*(x[n-2] + b[n-1])
				k += 1
				infNorm = infinityNorm(x-xReal)
				if infNorm < eps:
						print("\tGauss-Seidel:")
						print("\t\tK:", k)
						print("\t\tSpectral Radius:", exp(1/k * log(infNorm/ infinityNorm(xReal), 2)))
						return

def sor(n, b, xReal, w, eps):
		k=0
		x = np.zeros(n, dtype=np.float64)
		while True:
				x[0] = (1-w) * x[0] + w/2*(x[1] + b[0])
				for i in range(1, n-1):
						x[i] = (1-w) * x[i] + w/2*(x[i-1] + x[i+1] + b[i])
				x[n-1] = (1-w) * x[n-1] + w/2*(x[n-2] + b[n-1])
				k += 1
				infNorm = infinityNorm(x-xReal)
				if infNorm < eps:
						print("\tSOR:")
						print("\t\tK:", k)
						print("\t\tSpectral Radius:", exp(1/k * log(infNorm/ infinityNorm(xReal), 2)))
						return

def infinityNorm(x):
		maxVal = -1
		for val in x:
				if abs(val) > maxVal:
						maxVal = abs(val)
		return maxVal

if __name__ == "__main__":
		nList = [25, 50, 100, 200]
		wList = [1.78486, 1.88402, 1.93968, 1.96922]
		eps = 10e-6
		for n, w in zip(nList, wList):
				print("n =", n)
				h = 1/(n + 1)
				x0 = np.zeros(n, dtype=np.float64)
				xReal = np.zeros(n, dtype=np.float64)
				b = np.zeros(n, dtype=np.float64)
					
				#generating real solutions for this size matrix and the b vector
				for i in range(0, n):
						xReal[i] = (i+1)*h * (1. - (i+1) *h)
						b[i] = 2* (h**2)
				jacobi(n, b, xReal, eps)
				gaussSeidel(n,b,xReal, eps)
				sor(n,b,xReal, w, eps)
