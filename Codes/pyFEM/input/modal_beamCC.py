"""
modal_beamCC.py
"""

gmshfile = 'modal_beamCC.msh'

# materials: a list if dicts where each entry is one material
materials = [
    {'Name':'Mat', 'E': 1000.0, 'nu': 0.2, 'rho': 1},
]

# solid: simple list that associate physical sets to materials
solid = {
    'Tyre': 'Mat',   # index into materials
}


# Dirichlet BCs on edges. list of dicts
dbc = [
    {'physet': 'LeftEdge',  'dir': 1, 'value': 0.0},
    {'physet': 'LeftEdge',  'dir': 2, 'value': 0.0},
    {'physet': 'RightEdge', 'dir': 1, 'value': 0.0},
    {'physet': 'RightEdge', 'dir': 2, 'value': 0.0},
]

#  no dbcn
dbcn = []   # vals


# no tbc 
tbc = []  # vals





