"""
elasticity_example.py
"""

gmshfile = 'example.msh'

# materials: a list if dicts where each entry is one material
materials = [
    {'Name':'Mat1', 'E': 5.0, 'nu': 0.3},
    {'Name':'Mat2', 'E': 10.0,  'nu': 0.3},
]

# solid: simple dict that associate physical sets to materials
solid = {
    'LeftPart': 'Mat1',   # index into materials
    'RightPart': 'Mat2',
}

# Dirichlet BCs on edges. list of dicts
dbc = [
    {'physet': 'LeftEdge',  'dir': 1, 'value': 0.0},
    {'physet': 'RightEdge', 'dir': 1, 'value': 0.5},
]

# nodal Dirichlet BCs
dbcn = [
    {'node': 1, 'dir': 2, 'value': 0.0},
]

# traction (or similar) BCs
tbc = [
    {'physet': 'Hole', 'dir': 0, 'value': -2.0},
    {'physet': 'LowerEdge', 'dir': 2, 'value': -1.0},
    {'physet': 'UpperEdge', 'dir': 2, 'value': 1.0},
]

#The x constraints on two vertical sides prevent horizontal translation and in-plane rigid rotation. 
#Fixing one y displacement prevents vertical translation. The remaining system can have a unique solution.