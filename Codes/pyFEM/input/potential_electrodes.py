"""
potential_electrodes.py

Input for a scalar Laplace/potential problem in an air rectangle with two
empty rectangular electrodes facing one another.

Problem:
    Laplacian(phi) = 0 in the air domain

Boundary conditions:
    phi = 1 on Electrode1
    phi = 0 on Electrode0
    zero flux on the outer boundary, imposed naturally by leaving tbc empty

Convention:
    use dir = 1 for the scalar potential phi.
"""

gmshfile = 'potential_electrodes.msh'

materials = [
    {'Name': 'Air', 'k': 1.0},
]

solid = {
    'Air': 'Air',
}

# Essential BC: prescribed potential phi on the electrode boundaries.
dbc = [
    {'physet': 'Electrode1', 'dir': 1, 'value': 1.0},
    {'physet': 'Electrode2', 'dir': 1, 'value': 0.0},
]

# No nodal essential BCs are required.
dbcn = []

# Natural BC: zero flux on all boundaries where no Dirichlet BC is prescribed.
# Therefore the outer boundary is zero flux without adding any tbc entry.
tbc = []
