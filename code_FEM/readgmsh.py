def readgmsh(gmshfile, materials, solid, dbc, dbcn, tbc):
    from dataclasses import dataclass
    import numpy as np

    @dataclass
    class class_analysis:
        NN: int
        NT6: int
        NB3: int

    @dataclass
    class class_nodes:
        coor: float
        dof: int
        U: float

    @dataclass
    class class_T6:
        nodes: int
        material: int

    @dataclass
    class class_B3:
        nodes: int
        loadxy: float
        press: float

    def read_until(f, target):
        for line in f:
            if line.strip() == target:
                return True
        return False

    def read_physical_names(f):
        physical_names = {}
        if not read_until(f, '$PhysicalNames'):
            return physical_names

        nphys = int(f.readline().strip())
        for _ in range(nphys):
            parts = f.readline().strip().split(maxsplit=2)
            dim = int(parts[0])
            tag = int(parts[1])
            name = parts[2].strip().strip('"')
            physical_names[(dim, tag)] = name

        read_until(f, '$EndPhysicalNames')
        return physical_names

    def read_nodes(f):
        if not read_until(f, '$Nodes'):
            raise ValueError('Section $Nodes not found in msh file.')

        Nn = int(f.readline().strip())
        Nid = np.zeros(Nn, int)
        Ncoor = np.zeros((Nn, 2), float)

        for i in range(Nn):
            parts = f.readline().strip().split()
            Nid[i] = int(parts[0]) - 1
            Ncoor[i, :] = [float(parts[1]), float(parts[2])]

        read_until(f, '$EndNodes')

        Nmax = max(Nid) + 1
        coor = np.zeros((Nmax, 2), float)
        coor[Nid, :] = Ncoor

        nodes = class_nodes(
            coor,
            np.full((Nmax, 2), -1),
            np.zeros((Nmax, 2), float),
        )
        return Nmax, nodes

    def read_elements(f, physical_names):
        if not read_until(f, '$Elements'):
            raise ValueError('Section $Elements not found in msh file.')

        En = int(f.readline().strip())
        elements = []

        for _ in range(En):
            parts = f.readline().strip().split()
            etype = int(parts[1])
            ntags = int(parts[2])
            tags = [int(x) for x in parts[3:3 + ntags]]
            conn = [int(x) - 1 for x in parts[3 + ntags:]]

            phys_tag = tags[0] if ntags > 0 else None
            if etype == 8:
                phys_name = physical_names.get((1, phys_tag))
            elif etype == 9:
                phys_name = physical_names.get((2, phys_tag))
            else:
                phys_name = None

            elements.append({
                'type': etype,
                'phys_name': phys_name,
                'conn': conn,
            })

        return elements

    with open(gmshfile, 'r') as f:
        physical_names = read_physical_names(f)
        Nmax, nodes = read_nodes(f)
        elements = read_elements(f, physical_names)

    T6n = 0
    B3n = 0
    for elem in elements:
        if elem['type'] == 9:
            T6n += 1
        elif elem['type'] == 8:
            for bc in tbc:
                if bc['physet'] == elem['phys_name']:
                    B3n += 1
                    break

    T6nodes = np.zeros((T6n, 6), int)
    T6mat = np.zeros(T6n, int)
    B3nodes = np.zeros((B3n, 3), int)
    B3loadxy = np.zeros((B3n, 2), float)
    B3press = np.zeros(B3n, float)

    mat_map = {mat['Name']: i for i, mat in enumerate(materials)}

    it6 = 0
    ib3 = 0
    for elem in elements:
        phys_name = elem['phys_name']
        conn = elem['conn']

        if elem['type'] == 8:
            has_tbc = False
            for bc in tbc:
                if bc['physet'] == phys_name:
                    has_tbc = True
                    break

            if has_tbc:
                B3nodes[ib3, :] = conn[:3]
                for bc in tbc:
                    if bc['physet'] == phys_name:
                        idir = bc['dir']
                        value = bc['value']
                        if idir > 0:
                            B3loadxy[ib3, idir - 1] += value
                        else:
                            B3press[ib3] = value
                ib3 += 1
            else:
                for bc in dbc:
                    if bc['physet'] == phys_name:
                        idir = bc['dir']
                        value = bc['value']
                        nodes.dof[conn[:3], idir - 1] = -2
                        nodes.U[conn[:3], idir - 1] = value

        elif elem['type'] == 9:
            T6nodes[it6, :] = conn[:6]
            T6mat[it6] = mat_map[solid[phys_name]]
            it6 += 1

    for bc in dbcn:
        nodes.dof[bc['node'] - 1, bc['dir'] - 1] = -2
        nodes.U[bc['node'] - 1, bc['dir'] - 1] = bc['value']

    analysis = class_analysis(Nmax, T6n, B3n)
    T6 = class_T6(T6nodes, T6mat)
    B3 = class_B3(B3nodes, B3loadxy, B3press)

    return analysis, T6, B3, nodes
