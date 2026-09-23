#Requires expected section order and MSH2 ASCII structure

def readgmsh(gmshfile, materials, solid, dbc, dbcn, tbc):
    from dataclasses import dataclass
    import numpy as np

    @dataclass
    class class_analysis:
        NN: int #number of nodes
        NT6: int #number of T6
        NB3: int #number of active B3

    @dataclass
    class class_nodes:
        coor: float #coordinates of nodes
        dof: int #numbering of dofs
        U: float #value of dofs (solution)

    @dataclass
    class class_T6:
        nodes: int #connectivity of elements
        material: int #material id of elements

    @dataclass
    class class_B3: #for active B3 only
        nodes: int #connectivity of elements
        loadxy: float #cartesian traction components of elements
        press: float #normal pressure on elements
        
    ######################### see functions during execution 

    def read_until(f, target): #reads a file row by row until finding the target string
        for line in f:
            if line.strip() == target: #search condition: .strip() removes blank spaces at the start and end of the row, including (\n)
                return True #Returns True if target is found - check variable line
        return False #Returns False if target is not found

    def read_physical_names(f): #extracts the physical groups defined in the mesh
        physical_names = {} #initialize empty dictionary for physical_names
        if not read_until(f, '$PhysicalNames'): #search for PhysicalNames line in the msh file
            return physical_names

        nphys = int(f.readline().strip()) #msh file reader is now located at the number of physical_names
        for _ in range(nphys): #for each physical name store three parts: dimension, physical_tag id, name
            parts = f.readline().strip().split(maxsplit=2)
            dim = int(parts[0]) #0 Point, 1 Curve/line, 2 Surface, 3 Volume
            tag = int(parts[1]) #physical entity id
            name = parts[2].strip().strip('"')
            physical_names[(dim, tag)] = name #dictionary entries with the key dim-tag tuple and value the name

        read_until(f, '$EndPhysicalNames') #advance until file reading utile finding EndPhysicalNames
        return physical_names

    def read_nodes(f): #reads nodes list and store coordinates
        if not read_until(f, '$Nodes'): #search for Nodes line in the msh file
            raise ValueError('Section $Nodes not found in msh file.')

        Nn = int(f.readline().strip()) #msh file reader is now located at the number of nodes
        Nid = np.zeros(Nn, int) # initialize array of nodes id
        Ncoor = np.zeros((Nn, 2), float) # initialize array of nodes coordinates

        for i in range(Nn): #read each node 
            parts = f.readline().strip().split() #for node store three parts: id, x_coor, y_coor
            Nid[i] = int(parts[0]) - 1 #numbering in Python starts from 0
            Ncoor[i, :] = [float(parts[1]), float(parts[2])] #extract x,y coordinates

        read_until(f, '$EndNodes') #advance until file reading utile finding EndNodes

        Nmax = max(Nid) + 1 # maximum node tag, ideally coincides with Nn
        coor = np.zeros((Nmax, 2), float) 
        coor[Nid, :] = Ncoor #coincides with Ncoor
        # ensures that the coordinate-array rows agree with the node indices
        # this is unnecessary if nodes are listed already in the order (1,2,...,Nid)
        
        nodes = class_nodes(
            coor, #fill nodal coordinates field
            np.full((Nmax, 2), -1), #initialize a matrix with dimensions (Nmax, 2) full of -1 for dof numbering
            # we will use the dof field to distinguish prescribed and unknown nodal displacements
            np.zeros((Nmax, 2), float), #initialize a matrix with dimensions (Nmax, 2) full of 0 for U field (dof values)
        ) #note coordinate storage is sized by Nmax, allowing gaps but potentially wasting rows
        return Nmax, nodes

    def read_elements(f, physical_names): #extracts the elements type, tags and connectivity
        if not read_until(f, '$Elements'): #search for Elements line in the msh file
            raise ValueError('Section $Elements not found in msh file.')

        En = int(f.readline().strip()) #msh file reader is now located at the number of elements
        elements = [] #initialize empty list for elements

        for _ in range(En): #loop over En with silent counter
            parts = f.readline().strip().split()
            etype = int(parts[1]) #element type (8 for 3-node boundary element (B3); 9 for 6-node triangle (T6))
            ntags = int(parts[2]) #number of tags for element
            tags = [int(x) for x in parts[3:3 + ntags]] #tags values for element 
            #(list comprehension, e.g., myList = [i for i in range(10)])
            
            # First tag: physical entity id (e.g., 5 = LowerEdge), 
            # Second tag: geometrical entity id (e.g., 1 = geometrical curve 1)
            
            # A physical group can contain multiple geometrical entities (not a one-to-one mapping)
            conn = [int(x) - 1 for x in parts[3 + ntags:]] #0-based connectivty from the remaining entries
            # For B3 the midside node is last. 
            # For T6, the first three are vertices; the next three belong to edges 1-2, 2-3, 3-1

            phys_tag = tags[0] if ntags > 0 else None #conditional assignment
            # interrogate physical_names dictionaries
            if etype == 8:
                phys_name = physical_names.get((1, phys_tag)) #get value for corresponding 1D physical entity (name)
            elif etype == 9:
                phys_name = physical_names.get((2, phys_tag)) #get value for corresponding 2D physical entity (name)
            else:
                phys_name = None

            elements.append({
                'type': etype,
                'phys_name': phys_name,
                'conn': conn,
            }) #append a dictionary for each element to form a list of dictionaries

        return elements

######### execution starts here #########

    with open(gmshfile, 'r') as f: #open gmshfile in reading mode, and refer to it as "f"
        physical_names = read_physical_names(f) #create a dictionary physical_names by reading the msh file
        Nmax, nodes = read_nodes(f) #return "number" of nodes and nodes object
        elements = read_elements(f, physical_names) #return TEMPORARY list of element dictionaries {'type', 'phys_name', 'conn'}

    # info from gmsh file has been now imported
    
    T6n = 0 #counter for elements of T6 type
    B3n = 0 #counter for elements of B3 type (only if physical group is loaded)
    
    #find number of elments required to initialize data structures
    for elem in elements:
        if elem['type'] == 9:
            T6n += 1
        elif elem['type'] == 8: #BCs are imposed only through B3 elements
            for bc in tbc: #loop over each declared TBC
                #check if element belongs to a physical entity with TBC
                if bc['physet'] == elem['phys_name']:
                    B3n += 1 #If used for TBC, B3 element is counted as "active" 
                    break #If already counted, exit first loop, without wasting computation

    #allocate (once) arrays for all T6 elements and (active) B3 elements
    T6nodes = np.zeros((T6n, 6), int) #table for T6 connectivities
    T6mat = np.zeros(T6n, int) #table for materials associated with T6 elements
    B3nodes = np.zeros((B3n, 3), int) #table for B3 connectivities
    B3loadxy = np.zeros((B3n, 2), float) #table for caterian traction components on B3 elements
    B3press = np.zeros(B3n, float) #table for signed coefficients multplying an oriented normal load on B3 elements

    mat_map = {mat['Name']: i for i, mat in enumerate(materials)} #dictionary comprehension
    #rapidly creates a dictionary that associates each material name to its id in the list
    #mat_map is a search map where I enter with names from the material list and get the id
    #enumerate returns index and content while loopoing over materials entries

    it6 = 0
    ib3 = 0
    for elem in elements:
        phys_name = elem['phys_name']
        conn = elem['conn']
        
        ############ B3 elements ##############
        
        if elem['type'] == 8:
            has_tbc = False #inizialization
            for bc in tbc: #loop over each declared TBC
                #check if element belongs to a physical entity with TBC
                if bc['physet'] == phys_name:
                    has_tbc = True
                    break #exit the first loop, without wasting computation

            if has_tbc: #if element is associated with a TBC - active
                B3nodes[ib3, :] = conn[:3] #store connectivity
                for bc in tbc: #loop over all TBCs, there can be more than one on single element (for type and direction)
                    #check if the element has the current TBC
                    if bc['physet'] == phys_name:
                        idir = bc['dir'] #get direction of relevant TBC
                        value = bc['value'] #get value of relevant TBC
                        if idir > 0: #cartesian traction
                            B3loadxy[ib3, idir - 1] += value #store relevant xy load at correct direction
                        else: #press load
                            B3press[ib3] = value #store relevant press load at correct direction
                ib3 += 1 #increase counter to move to next row at subsequent loaded B3 element
                
            else: #inactive B3 elements are not kept in the B3 structure.
                for bc in dbc: #loop over all DBCs, there can be more than one on single element
                    #check if element belongs to a physical entity with DBC
                    if bc['physet'] == phys_name: 
                        #if yes, the corresponding dofs are marked as prescribed
                        idir = bc['dir'] #get direction of relevant DBC
                        value = bc['value'] #get value of relevant DBC
                        nodes.dof[conn[:3], idir - 1] = -2 #mark the three edge nodes in the DBC direction
                        nodes.U[conn[:3], idir - 1] = value #store imposed value
                        #This applies DBCs only if that boundary has no TBC
                        #In current version, TBC suppresses DBC on the same group

        ############ T6 elements ##############

        elif elem['type'] == 9:
            T6nodes[it6, :] = conn[:6] #store connectivity
            T6mat[it6] = mat_map[solid[phys_name]] #solid associates physical name to material name
            #first physical tag is translated to a material name, then to a material-list index
            
            it6 += 1 #increase counter to move to next row at subsequent T6 element

    for bc in dbcn: #loop over all nodal DBCNs, there can be more than one on single node
        nodes.dof[bc['node'] - 1, bc['dir'] - 1] = -2 #mark constrained dof
        nodes.U[bc['node'] - 1, bc['dir'] - 1] = bc['value'] #store imposed value

    analysis = class_analysis(Nmax, T6n, B3n) #N_nodes, N_T6, N_B3
    T6 = class_T6(T6nodes, T6mat) #connectivity, material for each T6
    B3 = class_B3(B3nodes, B3loadxy, B3press) #connectivity, xyBC, press_BC for each B3

    return analysis, T6, B3, nodes
