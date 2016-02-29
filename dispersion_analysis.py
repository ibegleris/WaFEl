from __future__ import division
from dolfin import *
import numpy as np
from scipy.constants import c,pi
from scipy.linalg import eig
from scipy import sparse
from scipy.sparse import linalg
import matplotlib.pylab as plt

a = 1.0
b = 0.5

mu_r = 1.0
n = 1
k = 0.1
lamda = 1.55e-2

e_r = (n+k*1j)**2

k0 = 2*pi/lamda

class epsilon_r(Expression):
    def eval(self, values, x):
        if x[0] > 0.4 and x[0]<0.6 and x[1] < 0.3 and x[1] >0.2:
            values[0] = 4.0
        else:
            values[0] = 1.0


mesh =RectangleMesh(Point(0,0),Point(a,b),10,5)