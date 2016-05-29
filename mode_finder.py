from __future__ import division#, print_function
import numpy as np
from scipy.constants import c,pi
from scipy.sparse.linalg import eigs, eigsh
from scipy.linalg import eig
from scipy.sparse import csr_matrix, lil_matrix, csc_matrix
import matplotlib.pylab as plt
from scipy.integrate import simps
import os
from matplotlib.colors import from_levels_and_colors
from dolfin import *
import time
from functions_dispersion_analysis import *


a = 2e-4
b = 2e-4

mu_r = 1.0 
lamda = 1.55e-6
r_core = 1e-5 # radius of core
r_clad = 10e-5 #radius of the fibre
nclad = 1.444#- 0.1e-4j# ref index of cladding
ncore = 1.445 - 1e-4j # ref index of core
#neff_g = 1.4445 # Guess of the modes
num= 10   #The number of modes guess 
neff_g= ncore
mesh_refinement = 0 # number of times to uniformly refine the mesh (used for convergence plots and better results)

k = is_loss(ncore,nclad)
if k ==0:
    V = 2*pi/lamda*r_core*(ncore**2 - nclad**2)**0.5
    print(V)
k0 = 2*pi/lamda

def ref(x,values = np.zeros(1)):
    point = (x[0]**2+ x[1]**2)**0.5
    if  point<= r_core:
        values[0] = ncore.real**2 - ncore.imag**2
    elif point > r_core and point <= r_clad:
        values[0] = nclad.real**2 - nclad.imag**2
    else:
        values[0] = 1.
    return values

def extinction(x,values = np.zeros(1)):
    point = (x[0]**2+ x[1]**2)**0.5
    if  point<= r_core:
        values[0] = 2*ncore.imag*ncore.real
    elif point > r_core and point <= r_clad:
        values[0] = 2*nclad.imag*ncore.real
    else:
        values[0] = 0
    return values


x = np.linspace(-a,a,512)
y = np.linspace(-b,b,512)
n_prof,k_prof = goemetry_plot(x,y,a,b,ref,extinction,nclad,ncore,r_core,r_clad)


class epsilon_real(Expression):
    def eval(self, values, x):
       values = ref(x,values)

class epsilon_imag(Expression):
    def eval(self, values, x):
       values = extinction(x,values)


mesh = gmesh_mesh("original_geometry.geo",a,b,r_core,r_clad,mesh_refinement)
#plot(mesh,interactive=True)


plot(mesh, interactive=True)

### Define the orders of the fucntion spaces for vector and nodal basis functions
vector_order = 2
nodal_order = 3


combined_space, A,B, A_complex,B_complex = Matrix_creation(mesh,epsilon_real,epsilon_imag,mu_r,k,k0,vector_order,nodal_order)


A,B,A_complex,B_complex,electric_wall = Mirror_boundary(mesh,combined_space,A,B,A_complex,B_complex,k)
#free_dofs = boundary_marker_locator(A,electric_wall)
free_dofs = boundary_marker_locator(A,electric_wall)


eigen,ev = find_eigenvalues(A,B,A_complex,B_complex,neff_g,num,k0,free_dofs,k,sparse_=True)


beta =1j*(eigen)**0.5 
beta = np.abs(np.real(beta)) -1j*np.imag(beta)

sort_index = np.argsort(beta.real)[::-1]

propagating_modes = np.where(((beta[sort_index]/k0).real>nclad.real) & ((beta[sort_index]/k0).real<ncore))
propagating_modes = propagating_modes[0][:]

print("The effective index of the most dominant modes are:")
print(beta[sort_index][propagating_modes]/k0)



size1 = 512
size2 = 512




mode = 1
mode_idx = propagating_modes[mode]
x,y,E,E_axial = electric_field_interpolation(x,y,r_core,r_clad,mode,mode_idx,k,beta,k0,A,ev,sort_index,free_dofs,combined_space,ploting = True,sp=20)


mode=2
mode_idx = propagating_modes[mode]
x,y,E,E_axial = electric_field_interpolation(x,y,r_core,r_clad,mode,mode_idx,k,beta,k0,A,ev,sort_index,free_dofs,combined_space,ploting = True,sp=20)






"""
def electric_field_interpolation(x,y,r_core,r_clad,mode,mode_idx,k,beta,k0,A,ev,sort_index,free_dofs,combined_space,ploting = False,sp=10):
    E,E_axial =electric_field_full()
    if ploting:
        plot_electric_field(x,y,E,E_axial,mode,mode_idx,beta,sort_index,k0,sp)
    return E,E_axial
"""









class interp_plot(object):
	def __init__(self,mode,size1,size2,min_max,propagating_modes,beta,sort_index,k0):
		self.mode = mode
		self.mode_idx = propagating_modes[self.mode]
		self.neff = beta[sort_index][mode_idx]/k0
		self.xmin, self.xmax,self.ymin,self.ymax = min_max
		
		self.x = np.linspace(xmin,ymax,size1)
		self.y = np.linspace(ymin,ymax,size2)
		self.E = None
		self.E_axial = None


		
	def electric_field_full(self,k,A,ev,sort_index,free_dofs,combined_space):
	    """
	    Releases the electric field from the calculated eigenvalus and eigen vectors
	    
	    Returns::
	    E[size,size,2],E_axial(Ez)
	    """

	    #post-process the coefficients to map back to the full matrix
	    coefficiants_global = np.zeros(A.size(0),dtype=np.complex)
	    coefficiants_global[free_dofs] = ev[:,sort_index[self.mode_idx]]
	    #Create a Function on the combined space
	    mode_re = Function(combined_space)
	    mode_im = Function(combined_space)
	    #Assign the coefficients of the function to the calculated values
	    mode_re.vector().set_local(np.real(coefficiants_global))
	    mode_im.vector().set_local(np.imag(coefficiants_global))
	    #Split the function into the parts in each of the functions spaces in combined_space
	    #This is done using DOLFINs Function.split()
	    (TE_re,TM_re) = mode_re.split()
	    (TE_im,TM_im) = mode_im.split()

	    E = np.zeros([len(x),len(y),2],dtype = np.complex)
	    E_axial = np.zeros([len(x),len(y)], dtype= np.complex)
	    for i,xx in enumerate(x):
	        for j,yy in enumerate(y):
	            point = Point(xx,yy)
	            E[i,j,:]     =  TE_re(point) + 1j*TE_im(point)
	            E_axial[i,j] =  TM_re(point) + 1j*TM_im(point)
	    self.E = E
	    self.E_axial = E_axial
	    self.mode_field = np.transpose((np.abs(self.E[:,:,0])**2 + np.abs(self.E[:,:,1])**2+np.abs(self.E_axial[:,:])**2)**0.5)
	    maxi = np.max(self.mode_field1)
	    self.mode_field1 /=maxi
	    return E,E_axial,self.mode_field	

	def plot_electric_field(self,sp=10,**kwrds):
		if self.E == None:
			electric_field_full(self,k,A,ev,sort_index,free_dofs,combined_space)
		fig = plt.figure(figsize=(7.0, 7.0))
		X,Y = np.meshgrid(self.x,self.y)
		plt.contourf(X,Y,self.mode_field1,90)
		plt.quiver(X[::sp,::sp], Y[::sp,::sp], self.E[::sp,::sp,0], self.E[::sp,::sp,1],headlength=7,scale = 500000)
		plt.xlabel(r'$x(m)$')
		plt.ylabel(r'$y(m)$')
		plt.title(r'mode$=$'+str(mode)+', '+'  $n_{eff}=$'+str(self.neff.real)+str(self.neff.imag)+'j')
		plt.show()
		return None