gmshfile = 'elasticity_plate_compress.msh'


# materials: a list if dicts where each entry is one material
materials = [
    {'Name':'Mat', 'E': 100.0, 'nu': 0.3},
]

# solid: simple list that associate physical sets to materials
solid = {
    'Plate': 'Mat'
}

# Dirichlet BCs on edges. list of dicts
dbc = [
    {'physet': 'UpperEdge',  'dir': 1, 'value': 0.0},
    {'physet': 'UpperEdge',  'dir': 2, 'value': -0.1},
    {'physet': 'LowerEdge', 'dir': 1, 'value': 0.0},
    {'physet': 'LowerEdge', 'dir': 2, 'value': 0.0},
]

# nodal Dirichlet BCs
dbcn = []

# traction (or similar) BCs
tbc = []



   
  


