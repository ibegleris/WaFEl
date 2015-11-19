# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%reset -f
from __future__ import division
import dolfin
from mshr import *
import numpy as np
sq_side1 = 15.
sq_side2 = 15.
core_r = 3.45
fibre_r = 5.0
n_core = 1.46
n_clad = 1.45
lamda = 1550e-9

clad_r = fibre_r - core_r
center1 = sq_side1/2
center2 = sq_side2/2
dolfin.set_log_level(dolfin.TRACE)

# Define 2D geometry
domain =   Rectangle(dolfin.Point(0., 0.), dolfin.Point(sq_side1, sq_side2)) \
         - Circle(dolfin.Point(center1, center2), core_r) \
         - Circle(dolfin.Point(center1, center2), fibre_r)
domain.set_subdomain(1, Rectangle(dolfin.Point(0., 0.), dolfin.Point(sq_side1, sq_side2)))
#domain.set_subdomain(1, Circle(dolfin.Point(center1, center2), core_r))
#omain.set_subdomain(2, Circle(dolfin.Point(center1, center2), fibre_r))
dolfin.info("\nVerbose output of 2D geometry:")
dolfin.info(domain, True)

# Generate and plot mesh
mesh2d = generate_mesh(domain, 45)
print mesh2d
dolfin.plot(mesh2d, "2D mesh")

# Convert subdomains to mesh function for plotting
mf = dolfin.MeshFunction("size_t", mesh2d, 2, mesh2d.domains())
dolfin.plot(mf, "Subdomains")

dolfin.interactive()

# <codecell>

# Generate and plot mesh
mesh2d = generate_mesh(domain, 45)
print mesh2d
dolfin.plot(mesh2d, "2D mesh")

# <codecell>

center

# <codecell>


