from functions_dispersion_analysis import *
from testing.Single_mode_fibre.Single_mode_theoretical import *
def test_main():
	mu_r = 1.0
	lamda = 1.55e-6
	r_core = 1e-5 # radius of core
	r_clad = 10e-5 #radius of the fibre
	nclad = 1.444#- 0.1e-4j# ref index of cladding
	ncore = 1.445# - 1e-4j # ref index of core
	neff_g = ncore # Guess of the modes
	num = 20   #The number of modes guess 

	lam_mult = 1 # Asks GMSH to create a mesh that has this number multiplied by the wavelength 
	mesh_refinement = 0 # number of times to uniformly refine the mesh (used for convergence plots and better results)
	vector_order = 3
	nodal_order = 3
	neff_th, Aeff_th = main_test(ncore,nclad,lamda,r_core,r_clad)


	waveguide = waveguide_inputs(lamda,ref,extinction)
	waveguide.fibre(r_core,r_clad,True,ncore,nclad)
	k = is_loss(ncore,nclad)
	box_domain = box_domains(a,b)
	modes_vec = main(box_domain,waveguide,vector_order,nodal_order,num,neff_g,lam_mult,k = 0,\
	                 size1 = 512,size2 = 512, mesh_plot = False,filename = 'geometry_test.geo')
	mesh = modes_vec[-1]
	modes_vec =modes_vec[:-1]
	for mode in modes_vec:
		abs_err_neff = np.abs(neff_th - mode.neff).real
		abs_err_aeff = 
	    assert np.abs(neff_th - mode.neff).real <= 1e-6
	    assert np.abs(neff_th - mode.neff).imag <= 1e-6