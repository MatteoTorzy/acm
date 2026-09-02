#
# file elasticity_plate_ellipse.py
#

gmshfile='elasticity_plate_ellipse.msh';


# materials: a list if dicts where each entry is one material
materials = [
    {'Name':'Mat', 'E': 1.0, 'nu': 0.3},
]

# solid: simple list that associate physical sets to materials
solid = {
    'Plate': 'Mat'
}

# Dirichlet BCs on edges. list of dicts
dbc = [
]

# nodal Dirichlet BCs
dbcn = [
    {'node': 1, 'dir': 1, 'value': 0.0},
    {'node': 1, 'dir': 2, 'value': 0.0},
    {'node': 2, 'dir': 2, 'value': 0.0},
]

# traction (or similar) BCs
tbc = [
    {'physet': 'LeftEdge', 'dir': 1, 'value': -1.0},
    {'physet': 'RightEdge', 'dir': 1, 'value': 1.0},
]




