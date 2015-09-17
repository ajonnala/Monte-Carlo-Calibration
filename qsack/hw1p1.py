# problem set 1
# iterative non-vectorized solution in python, for verification.
#
# Quang Sack
#
# Sept 13, 2015
#
#>>>>>>>>>>>>>>>>

import math
import random

# parameters
r_0 = 0.03
theta = 0.05
kappa = 0.4
beta = 0.05


def b(s):
  numerator = 1 - math.exp(-s * kappa)
  denominator = kappa
  return numerator/denominator


def a(s):
  value1 = theta
  value1 -= (beta**2) / (2 * (kappa**2))
  value2 = (s - b(s))
  value3 = ((beta**2) / (4*kappa)) 
  value3 *= b(s) ** 2
  return value1 * value2 + value3


def analyticPrice(t, T, X_t):
  return math.exp(-a(T-t) - b(T-t)*X_t)
  


def getBondPrice(t, T, iters, numSteps):
  total = 0
  for _ in range(iters):
    delta = float(t)/numSteps
    r = r_0
    for _ in range(numSteps):
      r = eulerStep(r, delta)
    total += analyticPrice(t, T, r)
  return total/iters


def eulerStep(X_prev, delta):
  epsilon = random.normalvariate(0, 1)
  X_next = kappa * (theta - X_prev) * delta
  X_next += beta * math.sqrt(delta) * epsilon
#  return X_next if X_next > 0 else 0
  return X_prev + X_next

# n: number of time steps
# t: initial time
# T: final time
def simulatePrice(n, t, T):
  delta = float(T-t)/n
  integral = r_0 * delta;
  X = r_0
  for i in range(n-1):
    # update X
    X = eulerStep(X, delta)
    # add path to integral
    integral +=  X * delta
  return math.exp((-1)*integral)


def monteCarlo(numSteps, iters, t, T):
  total = 0
  for i in range(iters):
    total += simulatePrice(numSteps, t, T)
  return total / iters



def solutionTable():
  M = 1000
  N = 1000
  print("\t\tAnalytic\t\tMonteCarlo")
  print("B(0, 1)\t\t" + str(analyticPrice(0, 1, r_0)) + "\t" + str(monteCarlo(N, M, 0, 1)))
  print("B(0, 2)\t\t" + str(analyticPrice(0, 2, r_0)) + "\t" + str(monteCarlo(N, M, 0, 2)))
  print("B(0, 3)\t\t" + str(analyticPrice(0, 3, r_0)) + "\t" + str(monteCarlo(N, M, 0, 3)))
  print("B(0, 4)\t\t" + str(analyticPrice(0, 4, r_0)) + "\t" + str(monteCarlo(N, M, 0, 4)))
  print("B(0, 5)\t\t" + str(analyticPrice(0, 5, r_0)) + "\t" + str(monteCarlo(N, M, 0, 5)))
  print("B(0, 6)\t\t" + str(analyticPrice(0, 6, r_0)) + "\t" + str(monteCarlo(N, M, 0, 6)))
  print("B(0, 7)\t\t" + str(analyticPrice(0, 7, r_0)) + "\t" + str(monteCarlo(N, M, 0, 7)))
  print("B(0, 8)\t\t" + str(analyticPrice(0, 8, r_0)) + "\t" + str(monteCarlo(N, M, 0, 8)))
  print("B(0, 9)\t\t" + str(analyticPrice(0, 9, r_0)) + "\t" + str(monteCarlo(N, M, 0, 9)))
  print("B(0, 10)\t" + str(analyticPrice(0, 10, r_0)) + "\t" + str(monteCarlo(N, M, 0, 10)))

solutionTable()





