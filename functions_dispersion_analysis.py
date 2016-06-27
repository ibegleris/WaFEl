from __future__ import division#, print_function
from dolfin import *
import numpy as np
from scipy.constants import c,pi
from scipy.sparse.linalg import eigs, eigsh
from scipy.linalg import eig
from scipy.integrate import simps,dblquad
from scipy.sparse import csr_matrix, lil_matrix, csc_matrix
import matplotlib.pyplot as plt
import os
from matplotlib.colors import from_levels_and_colors
from dolfin import *
import time

def gmesh_mesh(filename,a,b,r_core,r_clad,mesh_refinement,gmsh_ver = 'gmsh'):
    filename = os.path.join('fenics_mesh',filename) 
    with open(filename, 'r') as content_file:
        content = content_file.readlines()
    mesh_geom=os.popen(gmsh_ver + "-optimize_lloyd fenics_mesh/Output.geo -2 -o fenics_mesh/output_small.msh")
    print mesh_geom.read()
    new_content = []
    new_content.append('DefineConstant[ a = { '+str(a)+', Path "Gmsh/Parameters"}];\n')
    new_content.append('DefineConstant[ b = { '+str(b)+', Path "Gmsh/Parameters"}];\n')
    new_content.append('DefineConstant[ r_core = { '+str(r_core)+', Path "Gmsh/Parameters"}];\n')
    new_content.append('DefineConstant[ r_clad = { '+str(r_clad)+', Path "Gmsh/Parameters"}];\n')
    for i in range(4,len(content)):
        new_content.append(content[i])
    with open("fenics_mesh/Output.geo", "w") as text_file:
        for i in new_content:
            text_file.write(i)
    refine_list = ["output_small"]
    if mesh_refinement !=0:
        for i in range(mesh_refinement):
            refine_list.append("refine"+str(i+1))
        os.popen("rm fenics_mesh/refine*")
        for i in range(mesh_refinement):
            mesh_dolf = os.popen(gmsh_ver + " -refine fenics_mesh/"+str(refine_list[i])+".msh -o fenics_mesh/"+str(refine_list[i+1])+'.msh')
            print mesh_dolf.read()
            time.sleep(4)
    mesh_dolf = os.popen("dolfin-convert fenics_mesh/"+refine_list[-1]+".msh fenics_mesh/fibre_small.xml")
    print mesh_dolf.read()
    mesh = Mesh("fenics_mesh/fibre_small.xml")
    return mesh


def gmesh_mesh_new(filename,a,b,r_core,r_clad,mesh_refinement,lamda,num,gmsh_ver = 'gmsh'):
    filename = os.path.join('fenics_mesh',filename) 
    with open(filename, 'r') as content_file:
        content = content_file.readlines()
    
    
    new_content = []
    #new_content.append('DefineConstant[ a = { '+str(a)+', Path "Gmsh/Parameters"}];\n')
    #new_content.append('DefineConstant[ b = { '+str(b)+', Path "Gmsh/Parameters"}];\n')
    #new_content.append('DefineConstant[ r_core = { '+str(r_core)+', Path "Gmsh/Parameters"}];\n')
    #new_content.append('DefineConstant[ r_clad = { '+str(r_clad)+', Path "Gmsh/Parameters"}];\n')
    new_content.append('a = DefineNumber[ '+str(a)+', Name "Parameters/a" ];\n')
    new_content.append('b = DefineNumber[ '+str(b)+', Name "Parameters/b" ];\n')
    new_content.append('rcore = DefineNumber[ '+str(r_core)+', Name "Parameters/rcore" ];\n')
    new_content.append('rclad = DefineNumber[ '+str(r_clad)+', Name "Parameters/rclad" ];\n')
    new_content.append('lam = DefineNumber[ '+str(lamda)+', Name "Parameters/lam" ];\n')
    new_content.append('num = DefineNumber[ '+str(num)+', Name "Parameters/num" ];\n')
    
    for i in range(6,len(content)):
        new_content.append(content[i])
    with open("fenics_mesh/Output.geo", "w") as text_file:
        for i in new_content:
            text_file.write(i)
    refine_list = ["output_small"]
    mesh_geom = os.popen(gmsh_ver + " fenics_mesh/Output.geo -2 -o fenics_mesh/output_small.msh")
    print mesh_geom.read()
    if mesh_refinement !=0:
        for i in range(mesh_refinement):
            refine_list.append("refine"+str(i+1))
        os.popen("rm fenics_mesh/refine*")
        for i in range(mesh_refinement):
            mesh_dolf = os.popen(gmsh_ver + " -refine fenics_mesh/"+str(refine_list[i])+".msh -o fenics_mesh/"+str(refine_list[i+1])+'.msh')
            print mesh_dolf.read()
            time.sleep(4)
        


    mesh_dolf = os.popen("dolfin-convert fenics_mesh/"+refine_list[-1]+".msh fenics_mesh/fibre_small.xml")
    time.sleep(4)
    print mesh_dolf.read()
    mesh = Mesh("fenics_mesh/fibre_small.xml")
    return mesh




def geometry_plot(x,y,a,b,ref,extinction,nclad,ncore,r_core,r_clad):

    X,Y = np.meshgrid(x,y)
    n_plot = np.zeros(np.shape(X))
    k_plot = np.zeros(np.shape(X))
    for i,xx in enumerate(x):
        for j,yy in enumerate(y):
            n_plot[i,j] = ref([xx,yy])*10000000
            k_plot[i,j] = extinction([xx,yy])*10000000

    cmap1, norm1 = from_levels_and_colors([1,ref([r_core,0])*10000000,ref([r_clad,0])*10000000], ['blue', 'green','red'],extend='max')
    cmap2, norm2 = from_levels_and_colors([0,extinction([r_core,0])*10000000,extinction([r_clad,0])*10000000],  ['blue', 'green','red'],extend='max')

    fig = plt.figure(figsize=(20.0, 20.0))
    ax1 = fig.add_subplot(221)
    ax1.pcolormesh(X,Y,n_plot, cmap=cmap1, norm=norm1)
    ax1.set_title('real part profile')
    ax1.axis('equal')
    
    ax2 = fig.add_subplot(222)
    ax2.pcolormesh(X,Y,k_plot, cmap=cmap2, norm=norm2)
    ax2.set_title('Imaginary part  profile')
    ax2.axis('equal')
    return n_plot,k_plot


def scipy_sparse_eigensolver(A_np_sp,B_np_sp,neff_g,num,k0):
    """
    Uses the scipy eigs to calculate the eigenvalues and eigenvectors of the equation given an effective index guess and
    the number of modes needed.
    """
    eigen_g = -neff_g**2*k0**2
    eigen2, ev2 = eigs(A_np_sp,num,B_np_sp,sigma = eigen_g,which ='LM',v0=eigen_g*np.ones(np.shape(A_np_sp)[0]))
    
    return eigen2,ev2

def scipy_eigensolver(A_np,B_np):
    eigen, ev =  eig(A_np,B_np)
    return eigen, ev


def csr_creation(A,B,free_dofs):
    A_lil = lil_matrix((A.size(0), A.size(1)))
    B_lil = lil_matrix((B.size(0), B.size(1)))
    for i in range(A.size(1)):
        A_indices, A_values = A.getrow(i)
        B_indices, B_values = B.getrow(i)
        A_lil[A_indices, i] =  A_values[:, np.newaxis]
        B_lil[B_indices, i] =  B_values[:, np.newaxis]
    return A_lil.tocsr(), B_lil.tocsr()

from joblib import Parallel, delayed
def loop_hard(i,A,B,A_lil,B_lil):
    A_indices, A_values = A.getrow(i)
    B_indices, B_values = B.getrow(i)
    A_lil[A_indices, i] =  A_values[:, np.newaxis]
    B_lil[B_indices, i] =  B_values[:, np.newaxis]
    return 
#def get_rows(i,):
#    A_indices, A_values = A.getrow(i)
#    B_indices, B_values = B.getrow(i)
#    A_lil[A_indices, i] =  A_values[:, np.newaxis]
#    B_lil[B_indices, i] =  B_values[:, np.newaxis]
#    re






#def electric_field_interpolation(x,y,r_core,r_clad,mode,mode_idx,k,beta,k0,A,ev,sort_index,free_dofs,combined_space,ploting = False,sp=10):
#    E,E_axial =electric_field_full(mode_idx,x,y,k,A,ev,sort_index,free_dofs,combined_space)
#    if ploting:
#        plot_electric_field(x,y,E,E_axial,mode,mode_idx,beta,sort_index,k0,sp)
#    return x,y,E,E_axial



def conj_trans(A):
    return csr_matrix.conjugate(A).T


def is_loss(ncore,nclad):
    if (nclad.imag,ncore.imag) == (0,0):
        k = 0
    else:
        k = 1
    return k


def strip_boundary(free_dofs,A):
    A_np = A.array()[free_dofs,:][:,free_dofs]
    return A_np


def function_space(vector_order,nodal_order,mesh):
    "Define the function spaces"
    
    vector_space = FunctionSpace(mesh,"Nedelec 1st kind H(curl)",vector_order)
    nodal_space = FunctionSpace(mesh,"Lagrange",nodal_order)
    combined_space = vector_space*nodal_space
    return combined_space


def Matrix_creation(mesh,epsilon_real,epsilon_imag,mu_r,k,k0,vector_order = 2,nodal_order = 3):
    combined_space = function_space(vector_order,nodal_order,mesh)
    # Define the test and trial functions from the combined space here N_i and N_j are Nedelec 
    # basis functions and L_i and L_j are Lagrange basis functions
    (N_i,L_i) = TestFunctions(combined_space)
    (N_j,L_j) = TrialFunctions(combined_space)
    e_r_real = epsilon_real()
    
    s_tt_ij = 1.0/mu_r*inner(curl(N_i),curl(N_j))
    t_tt_ij = e_r_real*inner(N_i,N_j)
    s_zz_ij = (1.0/mu_r) * inner(grad(L_i),grad(L_j))
    t_zz_ij = e_r_real*inner(L_i,L_j)


    A_tt_ij = s_tt_ij - k0**2*t_tt_ij
    B_zz_ij = s_zz_ij - k0**2*t_zz_ij

    B_tt_ij = 1/mu_r*inner(N_i, N_j)
    B_tz_ij = 1/mu_r*inner(N_i, grad(L_j))

    B_zt_ij = 1/mu_r*inner(grad(L_i),N_j)
    #post-multiplication by dx will result in integration over the domain of the mesh at assembly time
    A_ij = A_tt_ij*dx
    B_ij = (B_tt_ij+B_tz_ij+B_zt_ij+B_zz_ij)*dx
    #assemble the system Matrices. If there is loss in the system then
    #we create a new set of matrixes and assemble them

    A = assemble(A_ij)
    B = assemble(B_ij)
    ####This is to try and introduce the complex part
    if k !=0:
        e_r_imag = epsilon_imag()
        A_ii_complex = e_r_imag*k0**2*inner(N_i,N_j)*dx
        B_ii_complex = e_r_imag*k0**2*inner(L_i,L_j)*dx
        A_complex = assemble(A_ii_complex)
        B_complex = assemble(B_ii_complex)
    else:
        A_complex, B_complex = None, None
    return combined_space, A,B, A_complex,B_complex


def Mirror_boundary(mesh,combined_space,A,B,A_complex,B_complex,k):
    boundary_markers = MeshFunction('size_t',mesh,1)
    boundary_markers.set_all(0)
    DomainBoundary().mark(boundary_markers,1)
    # Set zero electric field on the edges (electric wall) and mark the boundaries as 1
    electric_wall = DirichletBC(combined_space,Expression(("0.0","0.0","0.0"))
                            ,boundary_markers,1)
    # apply the boundary condition to the assembled matrices
    electric_wall.apply(A)
    electric_wall.apply(B)
    if k!=0:
        electric_wall.apply(A_complex)
        electric_wall.apply(B_complex)
    return A,B,A_complex,B_complex,electric_wall


def boundary_marker_locator(A,electric_wall):
    indicators = np.zeros(A.size(0))
    indicators[electric_wall.get_boundary_values().keys()]=1
    free_dofs = np.where(indicators == 0)[0]
    return free_dofs

def find_eigenvalues(A,B,A_complex,B_complex,neff_g,num,k0,free_dofs,k,sparse_,A_np=None,B_np=None):
    if not(sparse_):
        print('trying non sparse matrix')
        try:
            if k!=0:
                A_np = strip_boundary(free_dofs,A)
                B_np = strip_boundary(free_dofs,B)
                A_np_complex = strip_boundary(free_dofs,A_complex)
                B_np_complex = strip_boundary(free_dofs,B_complex)
                A_np = A_np+1j*A_np_complex
                B_np = B_np+1j*B_np_complex
            else:
                A_np = strip_boundary(free_dofs,A)
                B_np = strip_boundary(free_dofs,B)
            sparse_ = False
        except MemoryError:
            print "*****************The matrixes are way to large for this system.*****************"
            print "*********************The sparse Matrixes will now be tried**********************"
            sparse_ = True
            pass

    if sparse_:
        dot_sparse = csc_matrix.dot
        if A_np == None:       
            A_np, B_np = csr_creation(A,B,free_dofs)
            if k != 0:
                A_np_complex, B_np_complex = csr_creation(A_complex,B_complex,free_dofs)
                A_np += 1j*A_np_complex
                B_np += 1j*B_np_complex
                del A_np_complex,B_np_complex
            A_np = A_np[free_dofs,:][:,free_dofs]
            B_np = B_np[free_dofs,:][:,free_dofs]
        print "sparse eigenvalue time"
        eigen, ev = scipy_sparse_eigensolver(dot_sparse(conj_trans(B_np),A_np),dot_sparse(conj_trans(B_np),B_np),neff_g,num,k0)
    else:
        print("normal eigenvalue solver ")
        #eigen, ev = scipy_eigensolver(A_np,B_np)
        eigen, ev = scipy_eigensolver(np.dot(np.conjugate(B_np).T,A_np),np.dot(np.conjugate(B_np).T,B_np))
    return eigen,ev,A_np,B_np


def integration2d_simps(xx,yy,integrand):
    "integrated over two dimensions within the domain" 

    I = np.zeros(len(yy),dtype='complex')
    for i in range(len(yy)):
        I[i] = simps(integrand[i,:], yy)
    return simps(I,xx)

def effective_area_simps(E,x,y):
    integrand1 = (E[:,:,0].conjugate()*E[:,:,0] + E[:,:,1].conjugate()*E[:,:,1]).real   
    Over = integration2d_simps(x,y,integrand1)
        
    integrand2 = integrand1**2
    under = integration2d_simps(x,y,integrand2)
        
    return Over**2/under

def effective_area_simps(E,E_axial,x,y):
    integrand1 = np.conj(E[:,:,0])*E[:,:,0] + np.conj(E[:,:,1])*E[:,:,1] + np.conj(E_axial[:,:])*E_axial[:,:]    
    Over = integration2d_simps(x,y,integrand1)
        
    integrand2 = np.abs(np.abs(E[:,:,0])**2 + np.abs(E[:,:,1])**2 + np.abs(E_axial[:,:])**2)**2
    under = integration2d_simps(x,y,integrand2)
        
    return np.abs(Over)**2/under

def overlap_simps(En,E_axialn,Em,E_axialm,x,y):
    integrand1 = np.conjugate(En[:,:,0])*Em[:,:,0] + np.conjugate(En[:,:,1])*Em[:,:,1] + np.conjugate(E_axialn[:,:])*E_axialm[:,:]
    Over = integration2d_simps(x,y,integrand1)
    integrand2 = np.abs(En[:,:,0])**2 + np.abs(En[:,:,1])**2 + np.abs(E_axialn[:,:])**2
    under1 = integration2d_simps(x,y,integrand2)
       
    integrand3 = np.abs(Em[:,:,0])**2 + np.abs(Em[:,:,1])**2 + np.abs(E_axialm[:,:])**2
    under2 = integration2d_simps(x,y,integrand3)
    
    return np.abs(Over)**2/(under1*under2)

def Overlaps_simps(n,m,propagating_modes,x,y,r_core,r_clad,k,beta,k0,A,ev,sort_index,free_dofs,combined_space):
    """This function is set to calculate the overlaps between two modes using simpsons rule. If the modes are the same 
        then it calculates the effective area of the mode
        Inputs::
            E :
            E_axial :
            n,m (int,int) : The mode subscripts whose overlap integral is to be calculated. If 
                            n == m then the effective area of A_eff is calculated
        Local::
                
        Returns::
            Overlap of mode m and n
    """
    #Integrand = lambda x,y: E
    #propagating_modes[mode]
    #En,E_axialn = electric_field_interpolation(x,y,r_core,r_clad,propagating_modes[n],mode_idx,k,beta,k0,A,ev,sort_index,free_dofs,combined_space,False)
    En,E_axialn = electric_field_full(propagating_modes[n],x,y,k,A,ev,sort_index,free_dofs,combined_space)
    
    if n == m:
        res = effective_area_simps(En,E_axialn,x,y)
        return res
    else:
        Em,E_axialm = electric_field_full(propagating_modes[m],x,y,k,A,ev,sort_index,free_dofs,combined_space)
        res = overlap_simps(En,E_axialn,Em,E_axialm,x,y)
        return res
        #integrand1 = np.conjugate(En[:,:,0])*Em[:,:,0] + np.conjugate(En[:,:,1])*Em[:,:,1] + np.conjugate(E_axialn[:,:])*E_axialm[:,:]
        #Over = integration2d_simps(x,y,integrand1)
        #print(Over)
        ##integrand2 = np.abs(En[:,:,0])**2 + np.abs(En[:,:,1])**2 + np.abs(E_axialn[:,:])**2
        #under1 = integration2d_simps(x,y,integrand2)
        
        #integrand3 = np.abs(Em[:,:,0])**2 + np.abs(Em[:,:,1])**2 + np.abs(E_axialm[:,:])**2
        #under2 = integration2d_simps(x,y,integrand3)
    
        #print(integrand1,integrand2,integrand3)
        #print(under1,under2)
        #return np.abs(Over)**2/(under1*under2)


class modes(object):
    def __init__(self,mode,size1,size2,min_max,propagating_modes,beta,sort_index,k0):
        self.mode = mode
        self.mode_idx = propagating_modes[self.mode]
        self.neff = beta[sort_index][self.mode_idx]/k0
        
        self.xmin, self.xmax,self.ymin,self.ymax = min_max
        
        self.x = np.linspace(self.xmin,self.ymax,size1)
        self.y = np.linspace(self.ymin,self.ymax,size2)
        self.E = None
        self.E_axial = None


    def dolfin_functions(self,k,A,ev,sort_index,free_dofs,combined_space):
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
        self.TE_re = TE_re
        self.TE_im = TE_im
        self.TM_re = TM_re
        self.TM_im = TM_im
        return None#TE_re,TE_re,TM_re,TM_im

    def effective_area(self,k,A,ev,sort_index,free_dofs,combined_space,lim):
        try:
            temp = self.TE_re
        except AttributeError:
            self.dolfin_functions(k,A,ev,sort_index,free_dofs,combined_space)
            pass
        integrand1 = dblquad(self.Eabs2, -lim, lim, lambda x: -lim,lambda x: lim)
        integrand2 = dblquad(lambda y,x: self.Eabs2(y,x)**2, -lim, lim, lambda x: -lim,lambda x: lim)
        
        self.Aeff =  integrand1[0]**2/integrand2[0]


    def Eabs2(self,y,x):
        E_ = self.Efun(y,x)
        return (E_[0]*E_[0].conjugate() + E_[1]*E_[1].conjugate()).real
    def Efun(self,y,x):
        point = Point(x,y)
        E = self.TE_re(point)+1j*self.TE_im(point)
        return E[0],E[1]


    def effective_area_simps(self,k,A,ev,sort_index,free_dofs,combined_space):
        if self.E ==None:
            self.electric_field_full(k,A,ev,sort_index,free_dofs,combined_space)
        
        integrand1 = (self.E[:,:,0].conjugate()*self.E[:,:,0] + self.E[:,:,1].conjugate()*self.E[:,:,1]).real   
        Over = integration2d_simps(self.x,self.y,integrand1)
            
        integrand2 = integrand1**2
        under = integration2d_simps(self.x,self.y,integrand2)
        self.Aeff = Over**2/under    
        return Over**2/under
    def electric_field_full(self,k,A,ev,sort_index,free_dofs,combined_space):
        """
        Releases the electric field from the calculated eigenvalus and eigen vectors
        
        Returns::
        E[size,size,2],E_axial(Ez)
        """
        try:
            temp = self.TE_re
        except AttributeError:
            self.dolfin_functions(k,A,ev,sort_index,free_dofs,combined_space)
            pass
        
        E = np.zeros([len(self.x),len(self.y),2],dtype = np.complex)
        E_axial = np.zeros([len(self.x),len(self.y)], dtype= np.complex)
        for i,xx in enumerate(self.x):
            for j,yy in enumerate(self.y):
                point = Point(xx,yy)
                E[i,j,:]     =  self.TE_re(point) + 1j*self.TE_im(point)
                E_axial[i,j] =  self.TM_re(point) + 1j*self.TM_im(point)
        self.E = E
        self.E_axial = E_axial
        self.mode_field = np.transpose((np.abs(self.E[:,:,0])**2 + np.abs(self.E[:,:,1])**2+np.abs(self.E_axial[:,:])**2)**0.5)
        maxi = np.max(self.mode_field)
        self.mode_field /=maxi

        return None    

    def plot_electric_field(self,sp=10,scales = 500000,cont_scale=90,savefigs=False):

        fig = plt.figure(figsize=(7.0, 7.0))
        xplot = self.x*1e6
        yplot = self.y*1e6
        X,Y = np.meshgrid(xplot,yplot)
        try:
            plt.contourf(X,Y,self.mode_field,cont_scale)
        except AttributeError:
             raise NotImplementedError("interpolate before plotting")

        plt.quiver(X[::sp,::sp], Y[::sp,::sp], np.real(self.E[::sp,::sp,0]), np.real(self.E[::sp,::sp,1]),scale = scales,headlength=7)
        plt.xlabel(r'$x(\mu m)$')
        plt.ylabel(r'$y(\mu m)$')
        #plt.title(r'mode$=$'+str(self.mode)+', '+'  $n_{eff}=$'+str(self.neff.real)+str(self.neff.imag)+'j')
        if savefigs ==True:    
            plt.savefig('mode'+str(self.mode)+'.eps',bbox_inches ='tight')
        
            D = {}
            D['X'] = X
            D['Y'] = Y
            D['Z'] = self.mode_field
            D['u'] = np.real(self.E[::sp,::sp,0])
            D['v'] = np.real(self.E[::sp,::sp,1])
            D['scale'] = scales
            D['cont_scale'] = 90
            D['sp'] = sp
            savemat('mode'+str(self.mode)+'.mat',D)
        return None
