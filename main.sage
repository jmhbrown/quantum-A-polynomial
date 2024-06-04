class QuantumAPolynomial:
    def __init__(self, knot, pachner_moves = []):
        """
        parameters:
        knot str - the knot we're working on. e.g. 4_1, 6_2
        pachner_moves list[dict] - the list of pachner moves to apply, each in the format {'tet_num':int, 'face_index':int}.
        """

        self.knot = knot
        self.pachner_moves = pachner_moves

        self.longitude = None
        self.meridian = None
        self.gens_dict = None
        self.omega_with_q = None
        self.quotient_omega = None
        self.weights_dict = None
        self.crossing_relations = []
        self.quotient_lattice = None

        # Define q and some of its roots.
        q = var('q')
        qrt2 = var('qrt2')
        qrt4 = var('qrt4')
        q = qrt4^4
        q = qrt2^2
        qrt2 = qrt4^2


        import snappy
        M = snappy.Triangulation(knot)
        ## Do Pachner Moves Here.
        for move in pachner_moves:
            M._two_to_three(move['tet'],move['face'])

        self.M = M


    # ## Helper methods

    import functools

    def std_basis(i,dim,power=1):
        # returns (0,...,0,power,0,...,0) in ZZ^dim
        vec = [0]*dim
        vec[i] = power
        return vector(vec)

    # list + list is concatenation, here's a method to add lists termwise:
    # needs the library functools for reduce.
    def add_termwise(*args):
        return list(reduce(lambda x,y: vector(x)+vector(y), args))

    # multiply_termwise(v,w) returns (v1*w1,...,vn*wn)
    def multiply_termwise(*args):
        return reduce(lambda vec1, vec2: list(map(lambda x,y: x*y, vec1,vec2)), args)

    def dict_monomial_to_list(monomial):
        coeff = LaurentPolynomialRing(ZZ,'qrt4')(list(monomial.values())[0])
        if not coeff.is_monomial():
            raise Exception('I can only deal with coefficients in the form q^a')
        else:
            lattice_coord = list(list(monomial.keys())[0])
            lattice_coord[0] += coeff.degree()/2 # the coeff is in terms of qrt4, qrt4^4 = q
            return lattice_coord


    def names_to_lattice_coordinate(names,gens_dict):
        """
        Gets the lattice coordinate for a monomial. Assumes the coordinates are
        Weyl normal order, e.g. ['x','y','z'] -> :xyz:
        
        Parameters:
            names (list of strings or dictionary specifying exponents) - generator names
            gens_dict (dict) - maps generator names to lattice coordinates.
        
        Returns:
            vector - the lattice coordinate of the Weyl normal order product of the given generator.
        """
        if isinstance(names,dict): # means powers are specified
            return sum([power*gens_dict[gen] for gen,power in names.items()])
        else: # its a list or set
            return sum([gens_dict[gen] for gen in names])
        
    def lattice_coord_to_dict(coord,gens_dict):
        """
        Converts a lattice coordinate to a dictionary with variable names as keys and powers as values.
        Might assume that each generator's coordinate has exactly one non-zero entry.
        
        Parameters:
            coord (list or vector) - the coordinate of an element in the quantum torus
            gens_dict (dict) - the dictionary of generator lattice coordinates.
        
        Returns:
            dict - a dictionary representation of the coordinate.
        """
        dict_form = {}
        for gen_name,gen_coord in gens_dict.items():
            # divide by that entry to deal with qrt2, which is (2,0,...,0)
            # does not handle roots of other generators - they'll be set to zero by int()
            power = int(sum(ii*jj for (ii,jj) in zip(gen_coord,coord))/max(gen_coord))
            if power != 0:
                dict_form[gen_name] = power
        
        return dict_form

    def lattice_coord_to_ring_element(coord,the_ring,relations):
        """
        Turns a lattice coordinate into a ring element.
        Assumes that the coordinate and the_ring.gens() use the same order.
        This recreates the functionality of the_ring.monomial for a quantum torus.
        
        Parameters:
            coord (List) - the lattice coordinate of a module element.
            the_ring (fgp_module) - the ambient module
            relations (matrix) - the relations matrix
        
        Returns:
            (the_ring.element_class) - the ring element corresponding to the lattice coordinate
        """
        

        dim = len(vector(ZZ,coord))
        tmp_scalar_power = 0
        for i in range(dim):
            if not (coord[i].is_zero() or vector(coord[i+1:]).is_zero()):
                leading_factor = vector([0]*dim)
                leading_factor[i] = coord[i]
                trailing_terms = vector([0]*dim)
                trailing_terms[i+1:] = coord[i+1:]
            
                tmp_scalar_power -= (matrix(leading_factor)*relations*matrix(trailing_terms).transpose())[0]

        q_power = vector(list(tmp_scalar_power/2) + [0]*(dim-1)) 
        #logger.debug("The q_power is {}".format(q_power))
        scaled_coord = vector(coord) + vector(q_power)
        return reduce(the_ring.product, [term[0]**term[1] for term in zip(the_ring.gens(),scaled_coord)])

    def product_from_lattice_coordinates(*coords, relations=matrix.identity(3), gens_dict={'qrt2':(2,0,0),'A':(0,1,0),'a':(0,0,1)}):
        #TODO - double check that this works like I think it should.
        """
        Gets the lattice coordinate for the product X^coord1 * X^coord2 * ... * X^coordn
        This isn't simple termwise addition because lattice coordinates are in normal order, so we have a factor of q introduced.
        
        If X^coord1 * X^coord2 = q^(2n) X^coord2 * X^coord1, then 
        X^coord1 * X^coord2 = q^(-n) X^(coord1 + coord2)
        
        Parameters:
            coords (vector or list) - an arbitrary number (>=1) of lattice coordinates, in the order they're being multiplied
            relations (matrix) - the relations matrix
            gens_dict (dict) - dictionary with generator names. used to specify q.
        
        Returns:
            vector - the lattice coordinate for the product
        """
        
        get_q_power = lambda c1, c2: -(matrix(c1)*relations*matrix(c2).transpose())[0,0]
        
        return reduce(lambda c1,c2: get_q_power(c1,c2)*vector(gens_dict['qrt2'])/2+vector(c1)+vector(c2),coords)

    def product_from_names(*monomials,relations=matrix.identity(3),gens_dict={'qrt2':(2,0,0),'A':(0,1,0),'a':(0,0,1)}):
        """
        Returns the product of the monomials (in normal order) in the form of a dictionary.
        
        Parameters:
            monomials (dict or list of strings) - an arbitrary number (2 or more?) of dictionarys, in the order they're being multiplied.
            relations (matrix) - gives the commutation relations
            gens_dict (dict) - maps variable names to lattice coordinates
        
        Returns:
            dict - the product in normal order.
        """
        
        coords = [QuantumAPolynomial.names_to_lattice_coordinate(m,gens_dict) for m in monomials]
        
        return QuantumAPolynomial.product_from_lattice_coordinates(*coords,relations=relations,gens_dict=gens_dict)

    def clear_denominator(element):
        """Clears denominator by multiplying through by a monomial."""
        ambient_ring = element.parent()
        
        return ambient_ring.monomial(*[
            -min(powers) for powers in matrix(element.exponents()).transpose().rows()
        ])*element


    def merge_or_add(dict1, dict2):
        """
        A version of dict1 | dict2 which adds the values of repeated keys. Does not overwrite the original dictionaries.
        
        Example: (compare 'b' value)
        {'a':1,'b':-1} | {'b':1,'c':2}
        >> {'a': 1, 'b': 1, 'c': 2}
        merge_or_add({'a':1,'b':-1}, {'b':1,'c':2})
        >> {'a':1, 'b':0, 'c':2}
        
        Parameters:
        dict1 (dict)
        dict2 (dict)
        
        Returns:
        dict
        """
        
        union_keys = set(dict1.keys()) | set(dict2.keys())
        
        return {k : dict1.get(k,0)+dict2.get(k,0) for k in union_keys}


    # Here's a commented version explaining some about the output of `snappy.Triangulation._to_string()`.
    # 
    # ```
    # 
    # 3                                       # <-- number of tetrahedron
    # 
    # ### tetrahedron 1 ###
    #    1    2    1    2                     # <-- which tetrahedron face 0 1 2 3 goes too. 
    #  1302 3012 0132 0132                    # <-- permutations on the vertices that defines
    #    0    0    0    0                     #     exactly where the face goes.
    #   0  0  0  0  0  0  0  0  0  1  0 -1  0  0  0  0
    #   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    #   0 -1  0  1 -1  0  1  0  0  1  0 -1  0  1 -1  0
    #   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    # 0.0 0.0
    # 
    # ### tetrahedron 2 ###
    #    2    0    2    0 
    #  1302 2031 0132 0132
    #    0    0    0    0 
    #   0  0  0  0  1  0 -1  0  0  0  0  0  0  0  0  0
    #   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    #   0  0  0  0  2  0 -1 -1 -1  0  0  1  0  1 -1  0
    #   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    # 0.0 0.0
    # 
    # ### tetrahedron 3 ###
    #    0    1    0    1 
    #  1230 2031 0132 0132
    #    0    0    0    0 
    #   0  0  0  0 -1  0  0  1  0  0  0  0  0 -1  1  0
    #   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    #   0  1 -1  0 -1  0  0  1 -1  0  0  1  1 -2  1  0
    #   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
    # 0.0 0.0
    # ```





    # ### Lattice Coordinates for Short and Long Edges


    def get_unglued_gens_dict(M):
        """"
        Builds a dictionary describing the short and long edges in our tetrahedra.
        Threads are done later.
        Parameters:
        M snappy.Triangulation

        Returns:
        dict {str: vector } - the non-thread generators of our internal skein algebra.
        """

        num_tet = M.num_tetrahedra()

        # First q. The first lattice is
        gens_dict = {'qrt2': vector([1] + [0]*(num_tet*18) + [0]*(num_tet*12))}

        # Generators for the skein algebra of the tetrahedra. Threads (from gluing) come later.
        for t in range(M.num_tetrahedra()):
                        # qrt2 + earlier tets
            leading_zeros = 1 + 18*t
                        # later tets       + threads
            trailing_zeros = 18*(num_tet-1-t) + num_tet*12
            # short edges
            gens_dict.update(
                {                        # q + before             + short                     +  long + after & threads
                    "a{0}01".format(t) : vector([0]*leading_zeros + [1,0,0,0,0,0,0,0,0,0,0,0] + [0]*(6+trailing_zeros)),
                    "a{0}02".format(t) : vector([0]*leading_zeros + [0,1,0,0,0,0,0,0,0,0,0,0] + [0]*(6+trailing_zeros)),
                    "a{0}03".format(t) : vector([0]*leading_zeros + [0,0,1,0,0,0,0,0,0,0,0,0] + [0]*(6+trailing_zeros)),
                    "a{0}10".format(t) : vector([0]*leading_zeros + [0,0,0,1,0,0,0,0,0,0,0,0] + [0]*(6+trailing_zeros)),
                    "a{0}13".format(t) : vector([0]*leading_zeros + [0,0,0,0,1,0,0,0,0,0,0,0] + [0]*(6+trailing_zeros)),
                    "a{0}12".format(t) : vector([0]*leading_zeros + [0,0,0,0,0,1,0,0,0,0,0,0] + [0]*(6+trailing_zeros)),
                    "a{0}20".format(t) : vector([0]*leading_zeros + [0,0,0,0,0,0,1,0,0,0,0,0] + [0]*(6+trailing_zeros)),
                    "a{0}21".format(t) : vector([0]*leading_zeros + [0,0,0,0,0,0,0,1,0,0,0,0] + [0]*(6+trailing_zeros)),
                    "a{0}23".format(t) : vector([0]*leading_zeros + [0,0,0,0,0,0,0,0,1,0,0,0] + [0]*(6+trailing_zeros)),
                    "a{0}30".format(t) : vector([0]*leading_zeros + [0,0,0,0,0,0,0,0,0,1,0,0] + [0]*(6+trailing_zeros)),
                    "a{0}32".format(t) : vector([0]*leading_zeros + [0,0,0,0,0,0,0,0,0,0,1,0] + [0]*(6+trailing_zeros)),
                    "a{0}31".format(t) : vector([0]*leading_zeros + [0,0,0,0,0,0,0,0,0,0,0,1] + [0]*(6+trailing_zeros))  
                })
            # long edges
            gens_dict.update(
                {
                "A{0}01".format(t) : vector([0]*(leading_zeros+12) + [1,0,0,0,0,0] + [0]*trailing_zeros),
                "A{0}02".format(t) : vector([0]*(leading_zeros+12) + [0,1,0,0,0,0] + [0]*trailing_zeros),
                "A{0}03".format(t) : vector([0]*(leading_zeros+12) + [0,0,1,0,0,0] + [0]*trailing_zeros),
                "A{0}12".format(t) : vector([0]*(leading_zeros+12) + [0,0,0,1,0,0] + [0]*trailing_zeros),
                "A{0}13".format(t) : vector([0]*(leading_zeros+12) + [0,0,0,0,1,0] + [0]*trailing_zeros),
                "A{0}23".format(t) : vector([0]*(leading_zeros+12) + [0,0,0,0,0,1] + [0]*trailing_zeros)
                })

        return gens_dict

    

    # ### Threads and Gluing Data

    def get_short_edges_in_face(tet,face):
        """
        Returns a list of the short edges contained in a specified face.

        Parameters:
        tet (Integer) - index of the tetrahedron containing the face.
        face (Integer) - index of the face.

        Returns:
        list of str - list of short edges contained in the face.
        """
        edge_indices = [0,1,2,3]
        edge_indices.remove(face)
        return ["a{0}{1}{2}".format(tet,i,face) for i in edge_indices]

    def get_long_edges_in_face(tet,face):
        """
        Returns a list of the long edges contained in a specified face.

        Parameters:
        tet (Integer) - index of the tetrahedron containing the face.
        face (Integer) - index of the face.

        Returns:
        list of str - list of long edges contained in the face.
        """
        indices = [0,1,2,3]
        indices.remove(face)
        return ["A{0}{1}{2}".format(tet,indices[i1],indices[i2]) for i1 in range(3) for i2 in range(i1+1,3)]   
            
            
        return ["A{0}{1}{2}".format(tet,i,face) for i in edge_indices]

    def get_long_edges_in_tet(tet):
        """
        Returns a list of the long edges contained in a specified face.

        Parameters:
        tet (Integer) - index of the tetrahedron.

        Returns:
        list of str - list of long edges contained in the tet.
        """
        return [name.format(tet) for name in ['A{0}01','A{0}02', 'A{0}03','A{0}12','A{0}13','A{0}23']]
        
        
    def get_local_thread_weights(tet,face,weights_dict):
        """
        Finds the weights for the threads that would glue this face to another.

        Parameters:
        tet (Integer) - index of the tetrahedron containing the face
        face (Integer) - index of the face inside the tetrehedron.

        Returns:
        dict {str : Integer} - gives the weights of the threads for each vertex in the face."""
        short_edges = QuantumAPolynomial.get_short_edges_in_face(tet,face)
        local_thread_weights = {}
        for short_edge in short_edges:
            local_thread_weights.update({
                vertex: -weight for vertex,weight in weights_dict[short_edge].items()
            })

        return local_thread_weights


    # turn the gluing data into a dictionary.

    def get_gluing_dict(M):
        """
        Reformats the gluing data from a snappy triangulation with n tetrahedra.

        Parameters:
        M (snappy.Triangulation)

        Returns
        dict - each item is a map:
            'r{i}' specifies which tetrahedrons the faces of the i-th tet ends up in.
            's{i}{j}' specifies the gluing map for the j-th face of the i-th tet.
        """
        gluing_data = M._get_tetrahedra_gluing_data()
        gluing_dict = {}
        for t in range(M.num_tetrahedra()):
            gluing_dict.update({"r{0}".format(t) : gluing_data[t][0] })
            for f in range(4):
                gluing_dict.update({"s{0}{1}".format(t,f) : gluing_data[t][1][f]})
                
        return gluing_dict

    # This could be more beautiful.
    def get_thread_weights_dict(M,weights_dict):
        """
        Finds the weights_dict for the threads, based on a snappy Triangulation.
        
        Parameters:
        M (snappy.Triangulation)
        
        Returns:
        dict of the form {str: vector} - gives the weights of the T-actions for the threads.
        """
        gluing_dict = QuantumAPolynomial.get_gluing_dict(M)
        thread_weights_dict = {}
        for t in range(M.num_tetrahedra()):
            for f in range(4):
                face_perm = gluing_dict["s{0}{1}".format(t,f)]
                local_thread_weights = QuantumAPolynomial.get_local_thread_weights(t,f,weights_dict)
                # add threads based on where they start. i.e. where they have weight 1.
                starting_here = list(filter(lambda k : local_thread_weights[k] == 1, local_thread_weights.keys()))
                for local_vertex in starting_here:
                    distant_vertex = "v{0}{1}{2}".format(
                        gluing_dict["r{0}".format(t)][f],
                        face_perm[Integer(local_vertex[-2])],
                        face_perm[Integer(local_vertex[-1])]
                    )
                    thread_weights_dict.update({
                        "x{0}_{1}".format(local_vertex[1:], distant_vertex[1:]) : {local_vertex: 1, distant_vertex: -1}   
                    })

        return thread_weights_dict


    # #### Lattice Coordinates for Threads

    def add_thread_lattice_coordinates(M,gens_dict,weights_dict):
        """ Adds the threads to the gens_dict.
        Assumes the the thread weights have already been added. Modifies gens_dict.
        
        Parameters:
        M (snappy.Triangulation()) - the triangulated knot complement
        gens_dict (dict {str: vector}) - the lattice coordinates of the generators
        weights_dict (dict {str: vector})- the weights of the generators 
        """
        num_tet = M.num_tetrahedra()
        thread_names = [k for k in weights_dict.keys() if k[0] == 'x']
        for i in range(len(thread_names)):
            gens_dict.update({
                thread_names[i] : vector([0]*(1+18*num_tet+i) + [1] + [0]*(12*num_tet-i-1))
            })





    def get_edges_adjacent_to_vertex(vertex,weights_dict):
        """
        Returns a list of generators with non-zero weight at a given vertex.
        
        Parameters:
        vertex (str) - the name of a vertex.
        weights_dict (dict) - the weights dictionary.
        
        Returns:
        list of strings - the names of all edges with non-zero weight at this vertex.
        """
        adjacent_edges = []
        for edge, weights in weights_dict.items():
            if vertex in weights and weights[vertex] != 0:
                adjacent_edges.append(edge)
                
        return adjacent_edges

    def get_thread_relations(gens_dict,weights_dict):
        """
        Builds the part of the relations matrix concerning threads.
        """
        thread_names = list(filter(lambda k : k[0] == 'x', gens_dict.keys()))
        thread_relations = []
        
        for thread in thread_names:
            this_threads_relations = []
            for vertex,weight in weights_dict[thread].items():
                tmp_adjacent_edges = QuantumAPolynomial.get_edges_adjacent_to_vertex(vertex,weights_dict)
                # don't add relations for this thread with itself. 
                tmp_adjacent_edges.remove(thread)
                for edge in tmp_adjacent_edges:
                    # commutation relation is the same for all non-thread generators at this vertex.
                    if edge[0] == 'x':
                        this_threads_relations.append(weight*gens_dict[edge])
                    else:
                        this_threads_relations.append(-weight*gens_dict[edge])
            thread_relations.append(sum(this_threads_relations))
        
        return 2*matrix(thread_relations)

    def get_relations_matrix(M,gens_dict,weights_dict):
        """Constructs the full relations matrix.
        
        Parameters:
        M (snappy.Triangulation)
        gens_dict (dict)
        weights_dict (dict)
        
        Returns:
        Matrix
        """
        num_tet = M.num_tetrahedra()
        thread_relations = QuantumAPolynomial.get_thread_relations(gens_dict,weights_dict)
        
        omega_thread_non_thread = thread_relations[:,1:18*num_tet+1]
        omega_thread_thread = thread_relations[:,18*num_tet+1:]
        
        # ### Relations for Short and Long Edges.
        # first do short edges around a puncture
        omega_short_one_puncture = 2*matrix([
            [ 0, 1,-1],
            [-1, 0, 1],
            [ 1,-1, 0]
        ])
        omega_short = matrix.block_diagonal([omega_short_one_puncture]*4)
        # columns are At01, At02, At03, At12, At13, At23
        omega_short_long = 2*matrix([
            [0,1,1,0,0,0],
            [1,0,1,0,0,0],
            [1,1,0,0,0,0],
            #
            [0,0,0,1,1,0],
            [1,0,0,1,0,0],
            [1,0,0,0,1,0],
            #
            [0,0,0,1,0,1],
            [0,1,0,0,0,1],
            [0,1,0,1,0,0],
            #
            [0,0,0,0,1,1],
            [0,0,1,0,1,0],
            [0,0,1,0,0,1]
        ])

        # columns are short then long
        omega_one_tet = block_matrix([
            [omega_short,omega_short_long],
            [-omega_short_long.transpose(),0]
        ])

        omega_with_q = block_diagonal_matrix(
            matrix([0]),
            block_matrix([
                [block_diagonal_matrix(*[omega_one_tet]*num_tet),-omega_thread_non_thread.transpose()],
                [omega_thread_non_thread,omega_thread_thread]
            ])
        )
        
        return omega_with_q

    def get_weights_matrix(vertices_dict,weights_dict):
        """
        Builds the weights matrix. This is useful for finding invariant sublattices.
        
        Parameters:
        vertices_dict - dictionary {vertex name: vertex coordinate}
        weights_dict - dictionary {generator name : weight}
        
        Returns
        Matrix - the weights matrix. it's (number generators) x (number vertices)
        """
        # start off with qrt2 - this is weight zero.
        weights_list = [[0]*(len(vertices_dict.keys()))]
        
        for weights in weights_dict.values():
            weights_list.append(sum([
                weight*vector(vertices_dict[vertex]) for vertex,weight in weights.items()
            ]))

        return Matrix(weights_list)

    def get_weights_from_names(names,gens_dict,weights_matrix,vertices_dict,weight_format='coordinate'):
        """
        Gets the weight of a monomial.
        
        Parameters:
        names Dict - specifies powers of generators. {'gen':exponent}
        gens_dict Dict - translates between generator names and coordinates
        weights_matrix Matrix - gives the T-weights of generators
        
        weight_format string - either 'coordinate' or 'dict'. Defaults to 'coordinate'
        """
        weight_coord = weights_matrix.transpose()*QuantumAPolynomial.names_to_lattice_coordinate(names,gens_dict)
        if weight_format == 'coordinate':
            return weight_coord
        else:
            return QuantumAPolynomial.lattice_coord_to_dict(weight_coord,vertices_dict)


    def get_end_point_from_names(monomial,gens_dict,weights_matrix,vertices_dict):
        """
        Finds all the negative weight vertices of a monomial 

        Parameters:
        monomial dict - specifics a monomial in the quantum torus.
        gens_dict dict - translates between generator names and cooridnates
        weights_matrix matrix - gives the T-weights of generators
        vertices_dict dict - gives the vertice coordinates.
        
        Returns:
        list of strings - the vertices where this monomial has negative weight.

        """

        weights = QuantumAPolynomial.get_weights_from_names(monomial,gens_dict,weights_matrix,vertices_dict,weight_format='dict')
        end_points = [v for v,k in weights.items() if k < 0]
        if len(end_points) != 1:
            logger.warning("{0} has more or less than one end-point!".format(monomial))
        return end_points 

    def weight_at_vertex(vertex,monomial,gens_dict,weights_matrix,vertices_dict):
        """
        gets the weight of a specified monomial at a specified vertex.
        """

        weights = QuantumAPolynomial.get_weights_from_names(monomial,gens_dict,weights_matrix,vertices_dict,weight_format='dict')
        #logger.debug("weights for {2} at {0}: {1}".format(vertex,weights.get(vertex,0),monomial))
        return weights.get(vertex,0)


    def find_gen_with_certain_endpoint(vertex,monomial,gens_dict,weights_matrix,vertices_dict):
        candidates_in_monomial = {
                k:v for k,v in monomial.items() 
                if QuantumAPolynomial.weight_at_vertex(vertex,{k:v}, gens_dict,weights_matrix,vertices_dict) >= 1
                }
        #logger.debug('candidates in monomial: {}'.format(candidates_in_monomial))
        
        return candidates_in_monomial

    def order_curve(curve,gens_dict,weights_matrix,vertices_dict,relations):
        """
        Takes the normal-order version of an oriented curve and returns the product of the generators in the order that they're encountered as we travel along the curve. Makes an arbitrary choice of starting point."""
        #TODO - this doesn't handle long edges, which have weight +1 everywhere.

        # make a copy since we do descructive actions:
        unordered_curve = curve.copy()

        tmp_first = unordered_curve.popitem()
        next_gen = {tmp_first[0]:tmp_first[1]}
        
        ordered_curve = []
        ordered_curve.append(next_gen)

        while unordered_curve != {}:
            leftover_gen = {}
            next_vertex = QuantumAPolynomial.get_end_point_from_names(next_gen,gens_dict,weights_matrix,vertices_dict)[0]
            next_gen = QuantumAPolynomial.find_gen_with_certain_endpoint(next_vertex,unordered_curve,gens_dict,weights_matrix,vertices_dict)
            
            if len(next_gen) > 1:
                # We have more than one option for the next generator! Just remove one and stick it back in unordered_curve.
                logger.debug("More than one option for the next generator!")
                tmp_leftover_gen = next_gen.popitem()
                leftover_gen.update({tmp_leftover_gen[0]:tmp_leftover_gen[1]})
                # unused - this version picks the other direction.
                #tmp_next_gen = next_gen.popitem()
                #leftover_gen.update(next_gen)
                #next_gen = {tmp_next_gen[0]:tmp_next_gen[1]}
            if next_gen == {}:
                # we didn't find a candidate! This (hopefully!) means that we've completed an ordered curve and need to restart
                tmp_next_gen = unordered_curve.popitem()
                next_gen = {tmp_next_gen[0]:tmp_next_gen[1]}
                # re-add our new start since we remove it at the very end.
                unordered_curve.update(next_gen)
            if list(next_gen.values())[0] > 1:
                # This generator shows up more than once!
                # We only want to move one copy of it from unordered_curve to ordered_curve.
                leftover_gen.update({k:v-sign(v) for k,v in next_gen.items() if v-sign(v) != 0})
                next_gen = {k:sign(v) for k,v in next_gen.items()}
            
            ordered_curve.append(next_gen)
            if unordered_curve == {}:
                break
            else:
                unordered_curve.pop(list(next_gen.keys())[0])
                unordered_curve.update(leftover_gen)

        return QuantumAPolynomial.product_from_names(*ordered_curve,relations=relations,gens_dict=gens_dict)





    # ## Peripheral Curves

    def get_peripheral_curve_intersection_dict(M,curve='meridian'):
        """
        Turns the peripheral curve data from SnapPy into a dictionary,
        The values are the intersection number of the peripheral curve with the given edge.
        
        Parameters:
        M snappy.Triangulation - the triangulation of a knot complement.
        curve string - either 'meridian' or 'longitude'
        
        Returns:
        dict - of the form {'short_edge': weight}.
        """
        all_periph_data = M._get_cusp_indices_and_peripheral_curve_data()[1]
        data_start = {'meridian':0, 'longitude':2}
        
        curve_data = [all_periph_data[i] for i in range(data_start[curve],len(all_periph_data), 4)]
        curve_dict = {'a{0}{1}{2}'.format(t,v,f):curve_data[t][4*v+f]
                for t in range(M.num_tetrahedra())
                for v in range(4)
                for f in range(4)
                if 0 != curve_data[t][4*v+f]
            }
        #logger.debug('intersection data for {0}: {1}'.format(curve,str(curve_dict)))

        return curve_dict 


    def get_thread_going_to_vertex(vertex_name,gens_dict):
        """
        Returns the thread with weight -1 at the given vertex.
        """
        vertex_index = vertex_name[1:]
        threads_list = [k for k in gens_dict.keys() if k[0] == 'x' and k.split("_")[1] == vertex_index]
        assert len(threads_list) == 1
        return threads_list[0]

    def get_peripheral_curve_monomial(M,gens_dict,vertices_dict,weights_dict,weights_matrix,curve='meridian'):
        """
        Gets an expression for the peripheral curve of the triangulation in terms of generators.
        There are multiple correct answers. We find this one based on SnapPy's suggestion.
        
        Parameters:
        M snappy.Triangulation - the triangulation of a knot complement.
        gens_dict Dict - the dictionary of generators.
        vertices_dict Dict - the dictionary of vertices
        weights_dict Dict - the dictionary of weights for generators.
        weights_matrix Matrix - the matrix of weights for generators
        
        curve string - either 'meridian' or 'longitude'

        returns:
        dict - the lattice representation of the curve.
        """
        
        
        intersection_dict = QuantumAPolynomial.get_peripheral_curve_intersection_dict(M,curve=curve)
        
        #### First find the threads from the intersection dictionary.
        peripheral_curve_threads = {}
        
        vertices_for_positive_crossings = { # the weight +1 vertex for each short edge we cross positively.
            [vert for vert,weight in weights_dict[short_edge].items() if weight==1][0] : se_weight
            for short_edge,se_weight in intersection_dict.items() if se_weight > 0
        }
        peripheral_curve_threads.update({
            QuantumAPolynomial.get_thread_going_to_vertex(v,gens_dict) : weight for v,weight in vertices_for_positive_crossings.items()
        })
        
        #logger.debug("Threads in the {curve}: {threads}".format(curve=curve,threads=peripheral_curve_threads))
        #### Now we need to fill in the gaps with short edges.
        # Make a dictionary of all the short edges we'll need to include.
        peripheral_curve_short_edges = {}
        
        # we work based on the weights of the threads:
        peripheral_curve_thread_weights = QuantumAPolynomial.get_weights_from_names(peripheral_curve_threads,gens_dict,weights_matrix,vertices_dict,weight_format='dict')
        
        # work one boundary triangle at a time:
        for tet in range(M.num_tetrahedra()):
            for vertex in range(4):
                this_triangle_weights = {
                    k:v for k,v in peripheral_curve_thread_weights.items() if k[:-1] == 'v{0}{1}'.format(tet,vertex)
                }
                if not this_triangle_weights == dict(): # if there's weight zero move on.
                    #logger.debug("Face {vertex} of Tetrahedra {tet} has non-zero weight: {weights}".format(tet=tet,vertex=vertex,weights=this_triangle_weights))
                    local_short_edge_weights = {short_edge : weights_dict[short_edge] for short_edge in gens_dict.keys() if short_edge[:-1] == 'a{0}{1}'.format(tet,vertex)}
                    for short_edge,se_weights in local_short_edge_weights.items():
                        if set(se_weights.keys()) == set(this_triangle_weights.keys()):
                            a_vertex = list(se_weights.keys())[0]
                            short_edge_power = -se_weights[a_vertex]*this_triangle_weights[a_vertex]
                            #logger.debug("Adding {short_edge}^{power} to the curve.".format(short_edge=short_edge,power=short_edge_power))
                            old_short_edge_weight = peripheral_curve_short_edges.get(short_edge,0)
                            peripheral_curve_short_edges.update({short_edge:old_short_edge_weight+short_edge_power})


        curve_dict = (peripheral_curve_threads | peripheral_curve_short_edges)
        curve_weights = QuantumAPolynomial.lattice_coord_to_dict(QuantumAPolynomial.get_weights_from_names(curve_dict,gens_dict,weights_matrix,vertices_dict),vertices_dict)
        logger.debug("{0}: {1}".format(curve,curve_dict))

        if curve_weights != dict():
            logger.error("Curve has non-zero T-weights! {weights}".format(weights=curve_weights))
        return curve_dict
    

    # # Relations and Quotients
    # 
    # Any quotient of the form $M^{x} = 1$ in the quantum torus can be expressed in terms of the lattice as $x = 0$. We want to quotient by
    # 1. Monodromy at each of the punctures
    # 1. The gluing constraints, 6 from each face so 12*num_tetrahedron total
    # 1. Monodromy around the long edge handles
    # 1. The great circle relations, 1 from each tetrahedron so 2 total.
    # 
    # The first three types of relations are monomial. I expect to get a lattice of rank num_tetrahedron+2 once those are imposed. The last relation will have to be done in terms of the quantum torus, not the lattice.
    # 

    # ## Relations Lists

    # ### T-Monodromy Around Punctures




    def get_T_monodromy_list(M,gens_dict):
        """
        Constructs the list of T-monodromy expressions. Assumes that gens_dict has the short edges in the right order.
        
        Parameters:
        M (snappy.Triangulation) - the triangulated knot complement
        gens_dict (dict {str:vector}) - the lattice coordinates for the generators of the skein algebra
        
        Returns:
        list of lists of strings - each element is a list of the short edges around a puncture.
        
        """
        num_tet = M.num_tetrahedra()
        # the short edges are listed in gens_dict in the right order.
        short_edge_names = list(filter(lambda name: name[0]=='a', gens_dict.keys()))
        T_monodromy_list = []
        for p in range(num_tet*4):
            T_monodromy_list.append({'qrt2':2} | {k:1 for k in short_edge_names[3*p:3*p+3]})
        
        return T_monodromy_list

    # ### Gluing Relations




    def get_long_edge_gluing_relations_list(M,gens_dict):
        """
        Parameters:
        M (snappy.Triangulation) - the triangulated knot complement.
        gens_dict (dict {str:vector}) - the lattice coordinates for the generators of the skein algebra

        Returns
        list of dicts - gluing relations. keys are generators and values are their powers.
        """

        long_edge_gluing_relations_list = []
        
        # threads have names starting with the letter 'x'
        threads_list = list(filter(lambda gen : gen[0] == 'x',gens_dict.keys()))
        th = threads_list[0]
        tmp_gluing_dict = {th:-1}
        
        index_reversal_map = lambda index : index[:-2] + ''.join(reversed(index[-2:]))
        index_sort_map = lambda index : index[:-2] + ''.join(sorted(index[-2:]))
        
        while threads_list:
            other_thread = "x{1}_{0}".format(*[index_reversal_map(index) for index in th[1:].split('_')])
            tmp_gluing_dict[other_thread] = 1
            
            long_edge_pos = "A" + index_sort_map(th[1:].split('_')[0])
            tmp_gluing_dict[long_edge_pos] = 1
            
            long_edge_neg = "A" + index_sort_map(th[1:].split('_')[1])
            if long_edge_pos == long_edge_neg:
                del tmp_gluing_dict[long_edge_pos]
            else:
                tmp_gluing_dict[long_edge_neg] = -1
            
            long_edge_gluing_relations_list.append(tmp_gluing_dict)
            
            # remove these two threads so we don't check them again
            threads_list = list(set(threads_list) - set(tmp_gluing_dict.keys()))
            
            if threads_list == list():
                return long_edge_gluing_relations_list
            else:
                # on to the next one
                th = threads_list[0]
                tmp_gluing_dict = {th:-1}
                

    def get_short_edge_gluing_relations_list(M,weights_dict):
        """
        Constructs a the list of constraints inclured by gluing together short edges.
        
        Parameters:
        M snappy.Triangulation - the triangulated knot complement.
        weights_dict - the weights of edges. Used to find where edges start/end.
        
        Returns:
        list of lists of strings - elements are lists of generator names. Should be in the right order. 
        """
        gluing_data = QuantumAPolynomial.get_gluing_dict(M)
        short_edge_gluing_relations_list = []
        
        for tet in range(M.num_tetrahedra()):
            for face in range(4):
                new_tet = gluing_data['r{0}'.format(tet)][face]
                if new_tet <= tet: # avoid double-lisitng relations. will still double-list self-foldings.
                    gluing_perm = gluing_data['s{0}{1}'.format(tet,face)]
                    new_face = gluing_perm[face]
                    for short_edge in QuantumAPolynomial.get_short_edges_in_face(tet,face):
                        distant_short_edge = 'a{0}{1}{2}'.format(new_tet,gluing_perm[Integer(short_edge[-2])],new_face)

                        tmp_local_weights = weights_dict[short_edge]
                        local_starting_vertex = [vertex[1:] for vertex,weight in tmp_local_weights.items() if weight == -1][0]  
                        local_ending_vertex = [vertex[1:] for vertex,weight in tmp_local_weights.items() if weight == 1][0]  

                        tmp_distant_weights = weights_dict[distant_short_edge]
                        distant_starting_vertex = [vertex[1:] for vertex,weight in tmp_distant_weights.items() if weight == -1][0]  
                        distant_ending_vertex = [vertex[1:] for vertex,weight in tmp_distant_weights.items() if weight == 1][0]  

                        short_edge_gluing_relations_list.append({
                            'qrt2' : 2,
                            short_edge : 1,
                            'x{0}_{1}'.format(local_starting_vertex,distant_ending_vertex) : 1,
                            distant_short_edge : 1,
                            'x{0}_{1}'.format(distant_starting_vertex,local_ending_vertex) : 1
                        })
        return short_edge_gluing_relations_list


    def get_internal_edge_monodromy(gens_dict):
        """Constructs expressions for the T-region monodromies around the 'extra handles'
        on the glued surface.
        
        Works based on the thread names. Assumes that the monodromies are each a cycle.
        if they cross themselves or something this probably won't work."""
        
        list_of_monodromies = []
        
        threads_list = [k for k in gens_dict.keys() if k[0] == 'x']
        th = threads_list[0]
        tmp_monodromy_dict = {'qrt2':2,th:1}

        while threads_list:
            th_end_vertex = th[1:].split('_')[1]
            next_thread = list(filter(lambda thread : thread[1:].split('_')[0] == th_end_vertex,threads_list))[0]

            if next_thread in tmp_monodromy_dict.keys():
                # we've closed the loop!
                list_of_monodromies.append(tmp_monodromy_dict)
                threads_list = list(set(threads_list) - set(tmp_monodromy_dict.keys()))
                if threads_list == list():
                    return list_of_monodromies
                th = threads_list[0]
                tmp_monodromy_dict = {'qrt2':2,th:1}
            else:
                tmp_monodromy_dict.update({next_thread:1})
                th = next_thread
        

    def fix_long_edge_indices(dictionary):
        """
        Fixes the indices on the long edge names in a dictionary.
        These sometimes end up being backwards, e.g. A020 instead of A002.
        NOTE: this won't handle the situation where a dictionary has the same edge under two different names.
        
        Parameters:
        dictionary (dict) - the dictionary that needs fixing.
        
        Returns:
        dict - the fixed dictionary.
        
        """
        index_sort_map = lambda index : index[:-2] + ''.join(sorted(index[-2:]))
        
        return {index_sort_map(k) : dictionary[k] for k in dictionary.keys() if k[0] == 'A'} | {k : v for k,v in dictionary.items() if k[0] != 'A'}


    def compute_skein_module(__self__):
        M = __self__.M
        num_tet = M.num_tetrahedra()

        gens_dict = QuantumAPolynomial.get_unglued_gens_dict(M)

        # Vertices

        # make a basis for the vertices. This will be useful for the weights matrix.
        vertices_dict = {}
        for t in range(M.num_tetrahedra()):
            vertices_dict.update(
                {                       
                    "v{0}01".format(t) : [0]*(12*t) + [1,0,0,0,0,0,0,0,0,0,0,0] + [0]*(12*(num_tet-1-t)),
                    "v{0}02".format(t) : [0]*(12*t) + [0,1,0,0,0,0,0,0,0,0,0,0] + [0]*(12*(num_tet-1-t)),
                    "v{0}03".format(t) : [0]*(12*t) + [0,0,1,0,0,0,0,0,0,0,0,0] + [0]*(12*(num_tet-1-t)),
                    "v{0}10".format(t) : [0]*(12*t) + [0,0,0,1,0,0,0,0,0,0,0,0] + [0]*(12*(num_tet-1-t)),
                    "v{0}13".format(t) : [0]*(12*t) + [0,0,0,0,1,0,0,0,0,0,0,0] + [0]*(12*(num_tet-1-t)),
                    "v{0}12".format(t) : [0]*(12*t) + [0,0,0,0,0,1,0,0,0,0,0,0] + [0]*(12*(num_tet-1-t)),
                    "v{0}20".format(t) : [0]*(12*t) + [0,0,0,0,0,0,1,0,0,0,0,0] + [0]*(12*(num_tet-1-t)),
                    "v{0}21".format(t) : [0]*(12*t) + [0,0,0,0,0,0,0,1,0,0,0,0] + [0]*(12*(num_tet-1-t)),
                    "v{0}23".format(t) : [0]*(12*t) + [0,0,0,0,0,0,0,0,1,0,0,0] + [0]*(12*(num_tet-1-t)),
                    "v{0}30".format(t) : [0]*(12*t) + [0,0,0,0,0,0,0,0,0,1,0,0] + [0]*(12*(num_tet-1-t)),
                    "v{0}32".format(t) : [0]*(12*t) + [0,0,0,0,0,0,0,0,0,0,1,0] + [0]*(12*(num_tet-1-t)),
                    "v{0}31".format(t) : [0]*(12*t) + [0,0,0,0,0,0,0,0,0,0,0,1] + [0]*(12*(num_tet-1-t))
                })



        # weights for long and short edges (threads are done later.)
        weights_dict = {}
        for t in range(M.num_tetrahedra()):
            weights_dict.update({
                'a{0}01'.format(t) : {'v{0}02'.format(t): 1, 'v{0}03'.format(t): -1},
                'a{0}02'.format(t) : {'v{0}03'.format(t): 1, 'v{0}01'.format(t): -1},
                'a{0}03'.format(t) : {'v{0}01'.format(t): 1, 'v{0}02'.format(t): -1},
                'a{0}10'.format(t) : {'v{0}13'.format(t): 1, 'v{0}12'.format(t): -1},
                'a{0}13'.format(t) : {'v{0}12'.format(t): 1, 'v{0}10'.format(t): -1},
                'a{0}12'.format(t) : {'v{0}10'.format(t): 1, 'v{0}13'.format(t): -1},
                'a{0}20'.format(t) : {'v{0}21'.format(t): 1, 'v{0}23'.format(t): -1},
                'a{0}21'.format(t) : {'v{0}23'.format(t): 1, 'v{0}20'.format(t): -1},
                'a{0}23'.format(t) : {'v{0}20'.format(t): 1, 'v{0}21'.format(t): -1},
                'a{0}30'.format(t) : {'v{0}32'.format(t): 1, 'v{0}31'.format(t): -1},
                'a{0}32'.format(t) : {'v{0}31'.format(t): 1, 'v{0}30'.format(t): -1},
                'a{0}31'.format(t) : {'v{0}30'.format(t): 1, 'v{0}32'.format(t): -1}
            })
            
            weights_dict.update({
                'A{0}01'.format(t) : {'v{0}01'.format(t): 1, 'v{0}10'.format(t): 1},
                'A{0}02'.format(t) : {'v{0}02'.format(t): 1, 'v{0}20'.format(t): 1},
                'A{0}03'.format(t) : {'v{0}03'.format(t): 1, 'v{0}30'.format(t): 1},
                'A{0}12'.format(t) : {'v{0}12'.format(t): 1, 'v{0}21'.format(t): 1},
                'A{0}13'.format(t) : {'v{0}13'.format(t): 1, 'v{0}31'.format(t): 1},
                'A{0}23'.format(t) : {'v{0}23'.format(t): 1, 'v{0}32'.format(t): 1},
            })

        # are there self gluings?
        if logger.level <= 10:
            tet_gluings = {k:v for k,v in QuantumAPolynomial.get_gluing_dict(M).items() if k[0] == 'r'}
            for k,v in tet_gluings.items():
                this_tet = int(k[1:])
                if v.count(this_tet) != 0:
                    logger.debug("Tet {0} is self-glued!".format(this_tet))

        # Add thread weights to the big weights dictionary.
        weights_dict.update(QuantumAPolynomial.get_thread_weights_dict(M,weights_dict))
        # expose this as an external variable
        __self__.weights_dict = weights_dict

        # add threads to the gens_dict
        QuantumAPolynomial.add_thread_lattice_coordinates(M,gens_dict,weights_dict)
        __self__.gens_dict = gens_dict 



        # check that the thread relations are skew symmetric
        if logger.level <= 30:
            thread_relations = QuantumAPolynomial.get_thread_relations(gens_dict,weights_dict)
            # split the thread relations into two parts, to help construct the full relations matrix.
            omega_thread_non_thread = thread_relations[:,1:18*num_tet+1]
            omega_thread_thread = thread_relations[:,18*num_tet+1:]

            if not thread_relations[:,18*num_tet+1:].is_skew_symmetric():
                logger.warning("The thread relations are not skew-symmetric!")


        omega_with_q = QuantumAPolynomial.get_relations_matrix(M,gens_dict,weights_dict)
        __self__.omega_with_q = omega_with_q

        # checking that the keys are in the same order.
        if list(weights_dict.keys()) != list(gens_dict.keys())[1:]:
            logger.warning("The weights_dict and gens_dict have different key orders! This will probably cause trouble.")
        if not omega_with_q.is_skew_symmetric():
            logger.warning("The relations matrix is not skew symmetric!")

        #logger.debug("Rank of the center: {0}".format(omega_with_q.right_kernel_matrix().rank() ) )



        lattice = FreeModule(ZZ, len(gens_dict['qrt2']))

        weights_matrix = QuantumAPolynomial.get_weights_matrix(vertices_dict,weights_dict)
        invariant_sublattice = weights_matrix.left_kernel()


        ### Get the meridian and longitude

        unordered_meridian = QuantumAPolynomial.get_peripheral_curve_monomial(M,gens_dict,vertices_dict,weights_dict,weights_matrix,curve='meridian')
        meridian =QuantumAPolynomial.order_curve(unordered_meridian,gens_dict,weights_matrix,vertices_dict,omega_with_q)
        __self__.meridian = meridian

        unordered_longitude = QuantumAPolynomial.get_peripheral_curve_monomial(M,gens_dict,vertices_dict,weights_dict,weights_matrix,curve='longitude')
        longitude = QuantumAPolynomial.order_curve(unordered_longitude,gens_dict,weights_matrix,vertices_dict,omega_with_q)
        __self__.longitude = longitude


        T_monodromy_variable_names_list = QuantumAPolynomial.get_T_monodromy_list(M,gens_dict)
        #logger.debug("T_monodromy_variable_names_list:{}".format(T_monodromy_variable_names_list))
        T_monodromy_lattice_coordinate_list = [
            QuantumAPolynomial.names_to_lattice_coordinate(v,gens_dict) for v in T_monodromy_variable_names_list
        ]

        for coord in T_monodromy_lattice_coordinate_list:
            if not (matrix(coord)*weights_matrix).is_zero():
                logger.error("T-monodromy is not invariant, and it should be! "+ str(coord))



        long_edge_gluing_relations_list = QuantumAPolynomial.get_long_edge_gluing_relations_list(M,gens_dict)
        #logger.debug("long_edge_gluing_relations_list:\n"+"\n".join([str(relation) for relation in long_edge_gluing_relations_list]))

        # Check that all relations are T-invariant, and incidentally that we've listed actual threads.
        for constraint in long_edge_gluing_relations_list:
            lat_coord = QuantumAPolynomial.names_to_lattice_coordinate(constraint,gens_dict)
            if not (lat_coord*weights_matrix).is_zero():
                logger.error("A gluing relation is not T-invariant or involves non-existant threads! " + str(lat_coord*weights_matrix))


        short_edge_gluing_relations_list = QuantumAPolynomial.get_short_edge_gluing_relations_list(M,weights_dict)
        #logger.debug("short_edge_gluing_relations_list:\n"+"\n".join([str(relation) for relation in short_edge_gluing_relations_list]))

        # Check that all relations are T-invariant, and incidentally that we've listed actual threads.
        for constraint in short_edge_gluing_relations_list:
            if not (sum([gens_dict[edge] for edge in constraint])*weights_matrix).is_zero():
                logger.error("A gluing relation is not T-invariant, or perhaps involves non-existant threads! "+ str(constraint))
        internal_edge_monodromy_list = QuantumAPolynomial.get_internal_edge_monodromy(gens_dict)
        #logger.debug("internal_edge_monodromy_list. These should come in pairs of equal length. \n"+"\n".join([str(relation) for relation in internal_edge_monodromy_list]))


        if len(internal_edge_monodromy_list) != 2*M.num_tetrahedra():
            # 2t because there are 2t pairs of glued faces, t-1 of which don't increase the genus.
            # so genus is g = 2t - (t-1) = t+1. We're interested in the T-region torus, so g-1 genus is killed using 2(g-1) = 2t constraints.
            logger.warn("There are {0} internal edge monodromy relations, we expect to have {1} if every face has been glued!".format(len(internal_edge_monodromy_list),2*M.num_tetrahedra()))
            
        for constraint in internal_edge_monodromy_list:
            if not (sum([gens_dict[edge] for edge in constraint])*weights_matrix).is_zero():
                logger.error("The follow internal edge monodromy constraint is not T-invariant! " +str(constraint))

        # ## Quotient Lattice

        # Take the quotient!
        quotient_lattice = invariant_sublattice.quotient(
            [QuantumAPolynomial.names_to_lattice_coordinate(v,gens_dict) for v in
            internal_edge_monodromy_list
            + short_edge_gluing_relations_list
            + long_edge_gluing_relations_list
            + T_monodromy_variable_names_list
            ]
        )

        __self__.quotient_lattice = quotient_lattice

        pi = quotient_lattice.coerce_map_from(quotient_lattice.V())


        # Next build the quotient basis!

        # this is the unordered verison.
        #T_region_basis = matrix([
            #pi(gens_dict['qrt2']), pi(QuantumAPolynomial.names_to_lattice_coordinate(meridian,gens_dict)),
            #pi(QuantumAPolynomial.names_to_lattice_coordinate(longitude,gens_dict))
            #])
        T_region_basis = matrix([
            pi(gens_dict['qrt2']), pi(meridian),
            pi(longitude)
            ])
        logger.debug("T_region_basis:\n{0}".format(T_region_basis))
        T_region_minor_indexes = [(i,j,k)
                for i in range(T_region_basis.ncols())
                for j in range(i+1,T_region_basis.ncols())
                for k in range(j+1,T_region_basis.ncols())
                ]

        minors = T_region_basis.minors(3)
        non_zero_minors = {}
        logger.debug("minors: {0}".format(minors))
        #TODO: What does it mean for us to need a negative minor?.
        for i in range(len(minors)):
            if minors[i] == 1 or minors[i] == -1: # we only need one minor
                non_zero_rows = set(range(T_region_basis.ncols()))-set(T_region_minor_indexes[i])
                logger.debug("non_zero_rows: {0}".format(non_zero_rows))
                num_new_rows = len(non_zero_rows)
                non_zero_entries = zip(range(num_new_rows),non_zero_rows)
                new_rows = zero_matrix(ZZ,nrows=num_new_rows, ncols=T_region_basis.ncols())
                for entry in non_zero_entries:
                    new_rows[entry] = 1
                break
            # store all non-zero minors in case we need to use extended Euclidean algorithm.
            if minors[i] != 0: 
                # TODO - implement this using xgcd.
                non_zero_minors.update({T_region_minor_indexes[i]:minors[i]})


        quotient_basis = block_matrix([[T_region_basis],[new_rows]]).transpose()
        logger.debug("quotient basis:\n{0}".format(quotient_basis))

        if logger.level <= 10:
            non_long_edge_generators = list(filter(lambda gen : gen[0] != 'A',gens_dict.keys()))
            non_long_edge_lattice = Matrix([gens_dict[gen] for gen in non_long_edge_generators]).row_module().intersection(invariant_sublattice)
            T_region_quotient_lattice = non_long_edge_lattice.quotient(
                [QuantumAPolynomial.names_to_lattice_coordinate(v,gens_dict) for v in
                internal_edge_monodromy_list
                + short_edge_gluing_relations_list
                + T_monodromy_variable_names_list
                ]
            )

            logger.debug("Full lattice rank: {0}".format(lattice.rank()))
            logger.debug("T-region (short-edges + threads) sublattice of that is rank: {0}".format(non_long_edge_lattice.intersection(lattice).rank()))
            logger.debug("T-invariant lattice is rank: {0}".format(invariant_sublattice.rank()))

            # Checks -
            logger.debug("Including q, the quotient lattice should be rank: {0}. It's rank {1}".format(M.num_tetrahedra()+2, quotient_lattice.ngens()))
            logger.debug("The T-region quotient lattice has {0} generators.".format(T_region_quotient_lattice.ngens()))

            logger.debug("Quotient basis matrix should have determinent 1. Det: {}".format(quotient_basis.det()))

            logger.debug("The generators:")
            # what do these generators look like?
            for g in quotient_lattice.gens():
                logger.debug(QuantumAPolynomial.lattice_coord_to_dict(g.lift(),gens_dict))
                
            logger.debug("\nT-region generators:")
            for g in T_region_quotient_lattice.gens():
                logger.debug(QuantumAPolynomial.lattice_coord_to_dict(g.lift(),gens_dict))
            

        relations_matrix = matrix([QuantumAPolynomial.names_to_lattice_coordinate(v,gens_dict) for v in
            internal_edge_monodromy_list
            + short_edge_gluing_relations_list
            + long_edge_gluing_relations_list
            + T_monodromy_variable_names_list
            ])

        logger.debug("The relations matrix has dimension {0} and rank {1}.".format(relations_matrix.dimensions(),relations_matrix.rank()))

        smith_form = relations_matrix.smith_form()[0]
        for i in range(smith_form.rank()):
            if smith_form[i,i] != 1:
                logger.warning("We have torsion in the quotient where we didn't expect it! {} != 1".format(smith[i,i]))

        quotient_ring = LaurentPolynomialRing(QQ,['qrt2','M','L'] + ['w{0}'.format(i-3) for i in range(3,quotient_lattice.ngens())]+ ['w{0}i'.format(i-3) for i in range(3,quotient_lattice.ngens())])
        polynomial_ring = PolynomialRing(QQ,['qrt2','M','L'] + ['w{0}'.format(i-3) for i in range(3,quotient_lattice.ngens())]+ ['w{0}i'.format(i-3) for i in range(3,quotient_lattice.ngens())])
        quotient_ring.inject_variables()
        __self__.polynomial_ring=polynomial_ring

        change_of_basis_matrix = quotient_basis.det()*quotient_basis.inverse()

        # make the relations matrix for the quotient
        quotient_omega = Matrix([[(Matrix(v.lift())*omega_with_q*Matrix(w.lift()).transpose())[0,0] for v in quotient_lattice.gens()] for w in quotient_lattice.gens()])
        __self__.quotient_omega = quotient_omega

        logger.debug("quotient_omega is:\n {}".format(quotient_omega))
        if not quotient_omega.is_skew_symmetric():
            logger.error("The quotient lattice does not have a skew symmetric bilinear form!")


        # ## Kernel for the Knot Complement
        # I'm assuming that the first generator in the quotient ring is just qrt2. I can check that here.
        if pi(gens_dict['qrt2'])[0] != 1:
            logger.error("The generator qrt2 isn't sent to (1,0,...,0) in the quotient! It's {0}".format(pi(gens_dict['qrt2'])))

        qrt2 = quotient_ring('qrt2')

        generic_crossing_relation = [ # this equals 1 and was found by hand
            # The +8 is to resolve crossings. If I got the power backward these should both be +0. 
            {'qrt2':1, 'A{t}13':-1, 'A{t}02':-1, 'A{t}03':1, 'a{t}32':1, 'a{t}01':1, 'A{t}12':1, 'a{t}23':1, 'a{t}10':1},
            {'qrt2':-1, 'A{t}13':-1, 'A{t}02':-1, 'A{t}01':1, 'a{t}03':-1,'a{t}12':-1, 'A{t}23':1, 'a{t}21':-1, 'a{t}30':-1}
        ]

        crossing_relations = [quotient_ring('w{0}'.format(i-3))*quotient_ring('w{0}i'.format(i-3))-1 for i in range(3,quotient_lattice.ngens())]
        for tt in range(M.num_tetrahedra()):
            specific_crossing_relation = sum([
                QuantumAPolynomial.lattice_coord_to_ring_element(change_of_basis_matrix*vector(pi(QuantumAPolynomial.names_to_lattice_coordinate({k.format(t=tt) : v 
                    for k,v in monomial.items()},gens_dict))),quotient_ring, quotient_omega)
                        for monomial in generic_crossing_relation
            ]).subs({quotient_ring('qrt2'):-1})+1

            logger.debug("Relation for tetrahedron {0}: {1}".format(tt,specific_crossing_relation))

            crossing_relations.append(polynomial_ring(QuantumAPolynomial.clear_denominator(specific_crossing_relation)))
            
        __self__.crossing_relations = crossing_relations

        logger.info("Relations:\n"+str(crossing_relations))
        A_poly_candidate = polynomial_ring.ideal(crossing_relations).elimination_ideal([polynomial_ring(str(g)) for g in quotient_ring.gens()[-2*(M.num_tetrahedra()-1):]]).gens()[0]
#        A_poly_candidate = polynomial_ring.ideal(crossing_relations).variety()
        # TODO debugging output to probe the elimination ideal
        if logger.level <= 10:
            logger.debug("Eliminating: {}".format(quotient_ring.gens()[-2*(M.num_tetrahedra()-1):]))
            #partial_elim_ideal  = polynomial_ring.ideal(crossing_relations).elimination_ideal([polynomial_ring(str(g)) for g in quotient_ring.gens()[-2*M.num_tetrahedra():]]).gens()[0]
            #logger.debug("I've eliminated everything but {0}: {1}".format(quotient_ring.gens()[-M.num_tetrahedra()+1:],partial_elim_ideal))
        if A_poly_candidate == 0:
            logger.warning("We have 0 for the A-polynomial!")
        else:
            logger.info("A-polynomial:\n"+str(A_poly_candidate.factor()))
            if __self__.knot == '4_1':
                if A_poly_candidate.reduce([polynomial_ring(L^2*polynomial_ring('M')^4+L*(-polynomial_ring('M')^8+polynomial_ring('M')^6+2*polynomial_ring('M')^4+polynomial_ring('M')^2-1)+polynomial_ring('M')^4)])==0:
                    logger.info('This candidate IS divisible by the A-polynomial known to KnotAtlas')
                else:
                    logger.warning('This candidate IS NOT divisible by the A-polynomial known to KnotAtlas')
                    if A_poly_candidate.reduce([polynomial_ring(L^2*polynomial_ring('M')^4-L*(-polynomial_ring('M')^8+polynomial_ring('M')^6+2*polynomial_ring('M')^4+polynomial_ring('M')^2-1)+polynomial_ring('M')^4)])==0:
                        logger.warning('However it is once we substitute L for (-L)')
            if __self__.knot == '6_2':
                if A_poly_candidate.reduce([polynomial_ring(L^5*polynomial_ring('M')^(26)+L^4*(-polynomial_ring('M')^(30)+2*polynomial_ring('M')^(28)-polynomial_ring('M')^(26)-2*polynomial_ring('M')^(24)+5*polynomial_ring('M')^(22)+5*polynomial_ring('M')^(20)-3*polynomial_ring('M')^(18))+L^3*(-polynomial_ring('M')^(28)+3*polynomial_ring('M')^(26)-polynomial_ring('M')^(24)-5*polynomial_ring('M')^(22)-3*polynomial_ring('M')^(20)+12*polynomial_ring('M')^(18)+13*polynomial_ring('M')^(16)-3*polynomial_ring('M')^(14)-8*polynomial_ring('M')^(12)+3*polynomial_ring('M')^(10))+L^2*(3*polynomial_ring('M')^(20)-8*polynomial_ring('M')^(18)-3*polynomial_ring('M')^(16)+13*polynomial_ring('M')^(14)+12*polynomial_ring('M')^(12)-3*polynomial_ring('M')^(10)-5*polynomial_ring('M')^(8)-polynomial_ring('M')^6+3*polynomial_ring('M')^4-polynomial_ring('M')^2)+L*(-3*polynomial_ring('M')^(12)+5*polynomial_ring('M')^(10)+5*polynomial_ring('M')^8-2*polynomial_ring('M')^6-polynomial_ring('M')^4+2*polynomial_ring('M')^2-1)+polynomial_ring('M')^4)])==0:
                    logger.info('This candidate IS divisible by the A-polynomial known to KnotAtlas')
                else:
                    logger.warning('This candidate IS NOT divisible by the A-polynomial known to KnotAtlas')
                    if A_poly_candidate.reduce([polynomial_ring(-L^5*polynomial_ring('M')^(26)+L^4*(-polynomial_ring('M')^(30)+2*polynomial_ring('M')^(28)-polynomial_ring('M')^(26)-2*polynomial_ring('M')^(24)+5*polynomial_ring('M')^(22)+5*polynomial_ring('M')^(20)-3*polynomial_ring('M')^(18))-L^3*(-polynomial_ring('M')^(28)+3*polynomial_ring('M')^(26)-polynomial_ring('M')^(24)-5*polynomial_ring('M')^(22)-3*polynomial_ring('M')^(20)+12*polynomial_ring('M')^(18)+13*polynomial_ring('M')^(16)-3*polynomial_ring('M')^(14)-8*polynomial_ring('M')^(12)+3*polynomial_ring('M')^(10))+L^2*(3*polynomial_ring('M')^(20)-8*polynomial_ring('M')^(18)-3*polynomial_ring('M')^(16)+13*polynomial_ring('M')^(14)+12*polynomial_ring('M')^(12)-3*polynomial_ring('M')^(10)-5*polynomial_ring('M')^(8)-polynomial_ring('M')^6+3*polynomial_ring('M')^4-polynomial_ring('M')^2)-L*(-3*polynomial_ring('M')^(12)+5*polynomial_ring('M')^(10)+5*polynomial_ring('M')^8-2*polynomial_ring('M')^6-polynomial_ring('M')^4+2*polynomial_ring('M')^2-1)+polynomial_ring('M')^4)])==0:
                        logger.warning('However it is once we substitute L for (-L)')

        return A_poly_candidate

    # ### Different bases for the quotient lattice (scratch!)



    def get_A_polynomial_using_basis(M,basis_matrix,gens_dict,print_steps=True):
        generic_crossing_relation = [ # this equals 1.
            {'qrt2':1, 'A{t}13':-1, 'A{t}02':-1, 'A{t}03':1, 'a{t}32':1, 'a{t}01':1, 'A{t}12':1, 'a{t}23':1, 'a{t}10':1},
            {'qrt2':-1, 'A{t}13':-1, 'A{t}02':-1, 'A{t}01':1, 'a{t}03':-1,'a{t}12':-1, 'A{t}23':1, 'a{t}21':-1, 'a{t}30':-1}
        ]
        tmp_change_of_basis_matrix = basis_matrix.det()*basis_matrix.inverse()
        crossing_relations = []
        for tt in range(M.num_tetrahedra()):
            specific_crossing_relation = sum([
                QuantumAPolynomial.lattice_coord_to_ring_element(tmp_change_of_basis_matrix*vector(pi(QuantumAPolynomial.names_to_lattice_coordinate({k.format(t=tt) : v 
                    for k,v in monomial.items()},gens_dict))),quotient_ring, quotient_omega) 
                        for monomial in generic_crossing_relation
            ]).subs({quotient_ring('qrt2'):-1})+1
            crossing_relations.append(polynomial_ring(QuantumAPolynomial.clear_denominator(specific_crossing_relation)))
        
        
        A_poly_candidate = polynomial_ring.ideal(crossing_relations).elimination_ideal([polynomial_ring(str(g)) for g in quotient_ring.gens()[-M.num_tetrahedra()+1:]]).gens()[0]
        
        if print_steps:
            print("Relations:\n",crossing_relations)
            print("A-polynomial:\n",A_poly_candidate.factor())
            
        return A_poly_candidate
        
        
    #A = QuantumAPolynomial.get_A_polynomial_using_basis(test_basis,gens_dict)

