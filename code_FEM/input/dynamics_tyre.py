gmshfile = 'dynamics_tyre.msh'

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

# parameters for time history
tf = 45.           # final time
#Dt = .0033285    # time step expl no contact
Dt = .00332       # time step stable
output_interval=.2    # output interval 

  
