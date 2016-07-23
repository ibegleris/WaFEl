from __future__ import division, print_function
import numpy as np
from scipy.constants import c,pi


class box_domain(object):
	def __init__(self,a,b):
		self.a = a
		self.b = b

class waveguide_inputs(object):
	def __init__(self,lam,refr,exti):
		self.lamda = lam
		self.ref = refr
		self.extinction = exti



class eigen_parameters(object):
	def __init__(self,num, neff_g)