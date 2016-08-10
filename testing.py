from functions_dispersion_analysis import *
from testing.Single_mode_fibre.Single_mode_theoretical import *
def main_t(a,b,lamda,r_core,r_clad,ncore,nclad,neff_g,num,lam_mult,mesh_refinement,ref,extinction):
	mu_r = 1.0
	
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
		assert np.abs(neff_th - mode.neff).real <= 1e-6
		assert np.abs(neff_th - mode.neff).imag <= 1e-6

def ref_single_step(x,values = np.zeros(1)):
    point = (x[0]**2+ x[1]**2)**0.5
    if  point<= r_core:
        values[0] = ncore.real**2 - ncore.imag**2
    elif point > r_core and point <= r_clad:
        values[0] = nclad.real**2 - nclad.imag**2
    else:
        values[0] = 1.
    return values

def extinction_single_step(x,values = np.zeros(1)):
    point = (x[0]**2+ x[1]**2)**0.5
    if  point<= r_core:
        values[0] = -2*ncore.imag*ncore.real
    elif point > r_core and point <= r_clad:
        values[0] = -2*nclad.imag*ncore.real
    else:
        values[0] = 0
    return values

def test_single_mode_real():
	lamda = 1.55e-6
	r_core = 1e-5 
	r_clad = 10e-5 
	nclad = 1.444
	ncore = 1.445
	neff_g = ncore
	num = 20
	a = 2e-4
	b = 2e-4
	lam_mult = 1
	mesh_refinement = 0
	main_t(a,b,lamda,r_core,r_clad,ncore,nclad,neff_g,
			num,lam_mult,mesh_refinement,ref_single_step,extinction_single_step)
if __name__ == '__main__':
	test_single_mode_real()