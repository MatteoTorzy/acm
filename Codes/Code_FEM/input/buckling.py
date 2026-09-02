
import numpy as np

gmshfile = 'buckling.msh'

# Load/displacement sequence
delta = 1.0 / 16.0
LambdaT = np.arange(0.0, 1.0 + 0.5 * delta, delta)

# Newton controls.
toll = 1.0e-6
maxiter = 30


# materials: list of dicts, same style as example.py.
materials = [
    {'Name': 'Mat1', 'E': 1.0, 'nu': 0.0},
]

# solid: associates physical surface sets to materials
solid = {
    'Solid': 'Mat1',
}

# Dirichlet BCs on physical edge sets.
dbc = [
    {'physet': 'Lower', 'dir': 1, 'value': 0.0},
    {'physet': 'Lower', 'dir': 2, 'value': 0.0},
    {'physet': 'Upper', 'dir': 1, 'value': 0.0},
    {'physet': 'Upper', 'dir': 2, 'value': -0.3},
]

# Nodal Dirichlet BCs
dbcn = []

# Traction BCs.
tbc = [
    {'physet': 'Load', 'dir': 1, 'value': 1.0e-4},
]

# Buckling comparison plot parameters, matching the MATLAB post-processing block.
hist_mid_node = 6     # MATLAB numbering, node 6, x-displacement
hist_disp_node = 4    # MATLAB numbering, node 4, y-displacement
H = 1.0
L = 20.0
