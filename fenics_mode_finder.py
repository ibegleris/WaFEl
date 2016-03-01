from __future__ import division
import numpy as np
from scipy.constants import c,pi
from scipy.sparse.linalg import eigs, eigsh
from scipy.linalg import eig
from scipy import sparse
import matplotlib.pylab as plt
from dolfin import *

#The box domain
a = 1.0e-4
b = 1.0e-4

#constants
mu_r = 1.0


#Inputs
n = 1
k = 0.
lamda = 1.55e-6
neff_g = 1.445 # Guess of the modes
num= 100    #The number of modes
neff_g -= k*1j
bound =True

r = 1e-5
nclad = 1.444# + 0j
ncore = 1.445# + 0j
