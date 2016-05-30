# coding: utf-8
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
from scipy.io import savemat

# The box domain:
a = 2e-4
b = 2e-4
# Inputs of the problem
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

for mesh_refinement in range(4,5):
    k = is_loss(ncore,nclad)
    #if k ==0:
    V = 2*pi/lamda*r_core*(ncore.real**2 - nclad.real**2)**0.5
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
    #n_prof,k_prof = goemetry_plot(x,y,a,b,ref,extinction,nclad,ncore,r_core,r_clad)


    class epsilon_real(Expression):
        def eval(self, values, x):
           values = ref(x,values)

    class epsilon_imag(Expression):
        def eval(self, values, x):
           values = extinction(x,values)


    # ## Mesh

    # Load the gmsh file and if asked for refine the mesh.

    mesh = gmesh_mesh("original_geometry.geo",a,b,r_core,r_clad,mesh_refinement)
    #plot(mesh,interactive=True)





    num_cells = mesh.num_cells()



    vector_order = 2
    nodal_order = 3

    # Define the forms (matrix elements) for dispersion analysis into the basis functions

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
    size1,size2 = 512,512
    min_max = (-3*r_core,3*r_core,-3*r_core,3*r_core)
    Aeff = []
    for i in range(propagating_modes):
        mode0 = modes(i,size1,size2,min_max,propagating_modes,beta,sort_index,k0)
        mode0.electric_field_full(k,A,ev,sort_index,free_dofs,combined_space)
        mode0.effective_area(k,A,ev,sort_index,free_dofs,combined_space,r_clad)
        Aeff.append(mode0.Aeff)
    dicti = {}
    dicti['neff'] = beta[sort_index][propagating_modes]/k0
    dicti['cells'] = num_cells
    
    savemat('convergence'+str(mesh_refinement)+'.mat',dicti)
"""
# ### Plot the results

# In[75]:

size1,size2 = 512,512
min_max = (-3*r_core,3*r_core,-3*r_core,3*r_core)


# In[76]:

mode0 = modes(0,size1,size2,min_max,propagating_modes,beta,sort_index,k0)
mode0.electric_field_full(k,A,ev,sort_index,free_dofs,combined_space)


# In[77]:

mode0.plot_electric_field(scales = 100000,sp=40)


# In[78]:

mode1 = modes(1,size1,size2,min_max,propagating_modes,beta,sort_index,k0)
mode1.electric_field_full(k,A,ev,sort_index,free_dofs,combined_space)


# In[79]:

mode1.plot_electric_field(scales = 900000,sp=40)

"""
