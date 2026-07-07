gmshfile = 'dynamics_tyreincl.msh';

# materials: a list if dicts where each entry is one material
materials = [
    {'Name':'Mat', 'E': 500.0, 'nu': 0.2, 'rho': 1.0},
]

# solid: simple list that associate physical sets to materials
solid = {
    'Tyre': 'Mat',   # index into materials
}

# dbc along lines
dbc = []  # values

#  no dbcn
dbcn = []   # vals

# traction: here only used to define contact set
tbc = [
    {'physet': 'Border', 'dir': 0, 'value': 0.0},
]

#% initial velocity
v0 = -.2

# body forces val
bf = [0.0, -0.3]

# parameters for contact
Lambda = 100.       # penalty coefficient for contact
mu = 0.2             # friction coefficient

# parameters for time history
tf = 15.           # final time
Dt = .003       # time step stable
output_interval=.2    # output interval 

