class QuantumAPolynomial:
    def __init__(self, knot, pachner_moves = [], reverse_orientation=False, finish_variable_elimination=False, print_tex=False):
        """
        parameters:
        knot str - the knot we're working on. e.g. 4_1, 6_2
        pachner_moves list[dict] - the list of pachner moves to apply, each in the format {'tet_num':int, 'face_index':int}.
        """

        self.knot = knot
        self.pachner_moves = pachner_moves
        self.finish_variable_elimination = finish_variable_elimination

        self.longitude = None
        self.meridian = None
        self.omega_with_q = None
        self.quotient_omega = None
        self.classical_crossing_relations = []
        self.quotient_lattice = None
        self.ordered_quotient_lattice_gens = None
        self.pi = None
        self.polynomial_ring=None

        self.thread_monodromy_list = None 
        self.short_edge_gluing_relations_list = None 
        self.long_edge_gluing_relations_list = None 
        self.T_monodromy_variable_names_list = None 
        self.monomial_relations = None
        self.crossing_relations = None

        self.weights_matrix = None
        self.invariant_sublattice = None


        self.A_poly = None
        self.ref_A_poly = None

        with open('data/A_crossing/apolys/{}.txt'.format(self.knot), 'r') as file:
            self.ref_A_poly = file.read()

        qrt2 = var('qrt2')


        import snappy
        self.knot_comp = snappy.Triangulation(knot)
        ## Do Pachner Moves Here.
        for move in pachner_moves:
            self.knot_comp._two_to_three(move['tet'],move['face'])

        if reverse_orientation:
            self.knot_comp.reverse_orientation()

        self.num_tet = self.knot_comp.num_tetrahedra()

        self.gluing_dict = QuantumAPolynomial.get_gluing_dict(self.knot_comp)
        self.creased_edges = QuantumAPolynomial.get_creased_edges(self.gluing_dict)
        num_creased = len(self.creased_edges)
        logger.debug("Gluing dictionary: {}".format(self.gluing_dict))
        QuantumAPolynomial.tex_for_gluing_dict(self.gluing_dict)


        self.gens_dict = QuantumAPolynomial.get_unglued_gens_dict(self.knot_comp)

        # Vertices
        # make a basis for the vertices. This will be useful for the weights matrix.
        self.vertices_dict = {}
        for t in range(self.num_tet):
            self.vertices_dict.update(
                {                       
                    "v{0}01".format(t) : [0]*(12*t) + [1,0,0,0,0,0,0,0,0,0,0,0] + [0]*(12*(self.num_tet-1-t)),
                    "v{0}02".format(t) : [0]*(12*t) + [0,1,0,0,0,0,0,0,0,0,0,0] + [0]*(12*(self.num_tet-1-t)),
                    "v{0}03".format(t) : [0]*(12*t) + [0,0,1,0,0,0,0,0,0,0,0,0] + [0]*(12*(self.num_tet-1-t)),
                    "v{0}10".format(t) : [0]*(12*t) + [0,0,0,1,0,0,0,0,0,0,0,0] + [0]*(12*(self.num_tet-1-t)),
                    "v{0}13".format(t) : [0]*(12*t) + [0,0,0,0,1,0,0,0,0,0,0,0] + [0]*(12*(self.num_tet-1-t)),
                    "v{0}12".format(t) : [0]*(12*t) + [0,0,0,0,0,1,0,0,0,0,0,0] + [0]*(12*(self.num_tet-1-t)),
                    "v{0}20".format(t) : [0]*(12*t) + [0,0,0,0,0,0,1,0,0,0,0,0] + [0]*(12*(self.num_tet-1-t)),
                    "v{0}21".format(t) : [0]*(12*t) + [0,0,0,0,0,0,0,1,0,0,0,0] + [0]*(12*(self.num_tet-1-t)),
                    "v{0}23".format(t) : [0]*(12*t) + [0,0,0,0,0,0,0,0,1,0,0,0] + [0]*(12*(self.num_tet-1-t)),
                    "v{0}30".format(t) : [0]*(12*t) + [0,0,0,0,0,0,0,0,0,1,0,0] + [0]*(12*(self.num_tet-1-t)),
                    "v{0}32".format(t) : [0]*(12*t) + [0,0,0,0,0,0,0,0,0,0,1,0] + [0]*(12*(self.num_tet-1-t)),
                    "v{0}31".format(t) : [0]*(12*t) + [0,0,0,0,0,0,0,0,0,0,0,1] + [0]*(12*(self.num_tet-1-t))
                })



        # weights for long and short edges (threads are done later.)
        self.weights_dict = {}
        for t in range(self.num_tet):
            self.weights_dict.update({
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
            
            self.weights_dict.update({
                'A{0}01'.format(t) : {'v{0}01'.format(t): 1, 'v{0}10'.format(t): 1},
                'A{0}02'.format(t) : {'v{0}02'.format(t): 1, 'v{0}20'.format(t): 1},
                'A{0}03'.format(t) : {'v{0}03'.format(t): 1, 'v{0}30'.format(t): 1},
                'A{0}12'.format(t) : {'v{0}12'.format(t): 1, 'v{0}21'.format(t): 1},
                'A{0}13'.format(t) : {'v{0}13'.format(t): 1, 'v{0}31'.format(t): 1},
                'A{0}23'.format(t) : {'v{0}23'.format(t): 1, 'v{0}32'.format(t): 1},
            })



        # Add thread weights to the big weights dictionary.
        self.weights_dict.update(QuantumAPolynomial.get_thread_weights_dict(self.knot_comp,self.weights_dict))

        # add threads to the gens_dict
        QuantumAPolynomial.add_thread_lattice_coordinates(self.knot_comp,self.gens_dict,self.weights_dict)

        # check that the thread relations are skew symmetric
        if logger.level <= 30:
            thread_relations = QuantumAPolynomial.get_thread_relations(self.gens_dict,self.weights_dict)
            # split the thread relations into two parts, to help construct the full relations matrix.
            omega_thread_non_thread = thread_relations[:,1:18*self.num_tet+1]
            omega_thread_thread = thread_relations[:,18*self.num_tet+1:]

            if not thread_relations[:,18*self.num_tet+1:].is_skew_symmetric():
                logger.warning("The thread relations are not skew-symmetric!")


        logger.debug("gens: {}".format(self.gens_dict.keys()))
        # make the relations matrix!
        self.omega_with_q = 1/2*QuantumAPolynomial.get_relations_matrix(self.knot_comp,self.gens_dict,self.weights_dict)

        # checking that the keys are in the same order.
        if list(self.weights_dict.keys()) != list(self.gens_dict.keys())[1:]:
            logger.warning("The weights_dict and gens_dict have different key orders! This will probably cause trouble.")
        if not self.omega_with_q.is_skew_symmetric():
            logger.warning("The relations matrix is not skew symmetric!")

        #logger.debug("Rank of the center: {0}".format(self.omega_with_q.right_kernel_matrix().rank() ) )

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

    def tex_for_gluing_dict(gluing_dict):
        num_tet = len(gluing_dict)/5
        format_s = lambda s : " & s_{{{}}}".format(s[1:])
        format_r = lambda r : "\\\\ \\rho_{{{}}}".format(r[1:])
        list_to_perm_str = lambda l : " &= ({})".format(",".join([str(ii) for ii in l]))
        
        build_str = [(format_r(k) if k[0] == 'r' else format_s(k))+list_to_perm_str(v) for k,v in gluing_dict.items()]
        logger.debug("Gluing: {}".format(''.join(build_str)))

        

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
        
    def names_to_str(names,omega,gens_dict):
        """Converts the dictionary form of a variable to a string"""
        to_tex = []
        q_scalar = 2*QuantumAPolynomial.normal_to_lexi_order_scalar(names,omega,gens_dict)
        names['qrt2'] = names.get('qrt2',0) + q_scalar
        # make it print in lex order
        sorted_name_keys = list(names.keys())
        sorted_name_keys.sort(key=lambda val : gens_dict.get(val,0),reverse=True)
        sorted_names = {k:names[k] for k in sorted_name_keys}
        for name,power in sorted_names.items():
            if power != 0:
                if name[0] == 'q':
                    to_tex = ['q^{{{0}}}'.format(power/2)] + to_tex
                elif name[0] == 'a' or name[0] == 'A':
                    variable = '{0}_{{{1}}}'.format(name[0],name[1:])
                    if power == 1:
                        to_tex.append(variable)
                    else:
                        to_tex.append(variable+'^{{{0}}}'.format(power))
                elif name[0] == 'x':
                    variable = 'x^{{{0}}}_{{{1}}}'.format(*name[1:].split('_'))
                    if power == 1:
                        to_tex.append(variable)
                    else:
                        to_tex.append('(' + variable+')^{{{0}}}'.format(power))

        # TODO - this is converting between normal and lex ordering, so it needs a q-power!!
        #logger.debug("Tex expression: {}".format(' '.join(to_tex)))
        return ' '.join(to_tex)


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

    def normal_to_lexi_order_scalar(names,omega,gens_dict):
        coord = QuantumAPolynomial.names_to_lattice_coordinate(names,gens_dict)
        return (-1)*sum(flatten([
                [
                    coord[i]*coord[j]*omega[i,j] for i in range(j)
                ] for j in range(len(coord))
            ]))

    def lattice_coord_to_ring_element(coord,the_ring,omega):
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
        # match with the g algebra code
        second_version = (-1)*sum(flatten([
                [
                    coord[i]*coord[j]*omega[i,j] for i in range(j)
                ] for j in range(len(coord))
            ]))
        q_power = vector(list(second_version) + [0]*(dim-1)) 
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


    def get_unglued_gens_dict(knot_comp):
        """"
        Builds a dictionary describing the short and long edges in our tetrahedra.
        Threads are done later.
        Parameters:
        knot_comp snappy.Triangulation

        Returns:
        dict {str: vector } - the non-thread generators of our internal skein algebra.
        """

        num_tet = knot_comp.num_tetrahedra()

        # First q. The first lattice is
        gens_dict = {'qrt2': vector([1] + [0]*(num_tet*18) + [0]*(num_tet*12))}

        # Generators for the skein algebra of the tetrahedra. Threads (from gluing) come later.
        for t in range(knot_comp.num_tetrahedra()):
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

    def get_gluing_dict(knot_comp):
        """
        Reformats the gluing data from a snappy triangulation with n tetrahedra.

        Parameters:
        knot_comp (snappy.Triangulation)

        Returns
        dict - each item is a map:
            'r{i}' specifies which tetrahedrons the faces of the i-th tet ends up in.
            's{i}{j}' specifies the gluing map for the j-th face of the i-th tet.
        """
        gluing_data = knot_comp._get_tetrahedra_gluing_data()
        gluing_dict = {}
        for t in range(knot_comp.num_tetrahedra()):
            gluing_dict.update({"r{0}".format(t) : gluing_data[t][0] })
            for f in range(4):
                gluing_dict.update({"s{0}{1}".format(t,f) : gluing_data[t][1][f]})
                
        return gluing_dict

    def get_creased_edges(gluing_dict):
        """
        Looks for long edges that would be identified with themselves in the triangulation.

        Parameters:
        __self__ QuantumAPolynomial

        Returns:
        list[str] long edges that would be identified with themselves.
        """
        tet_gluings = {k:v for k,v in gluing_dict.items() if k[0] == 'r'}

        # first find self-glued tets.
        self_glued_faces = {}
        for k,v in tet_gluings.items():
            tet_index = int(k[1:])
            try:
                face_a = v.index(tet_index)
                face_b = v.index(tet_index,face_a + 1)
            except ValueError:
                # not self-glued
                pass
            else:
                self_glued_faces.update({tet_index: [face_a,face_b]})
                logger.warning("Tet {0} is self-glued at faces {1},{2}".format(tet_index, face_a,face_b))

        # next check if any of these have creased long edges
        creased_edges = []
        for tet,faces in self_glued_faces.items():
            if gluing_dict["s{0}{1}".format(tet,faces[0])] == gluing_dict["s{0}{1}".format(tet,faces[1])]:
                creased_edge_index = [str(x) for x in {0,1,2,3}.difference(set(faces))]
                creased_edges.append("A{0}{1}".format(tet,"".join(creased_edge_index)))
                logger.warning("Long edge {0} is self-identified!".format(creased_edges[-1]))

        return creased_edges


    # This could be more beautiful.
    def get_thread_weights_dict(knot_comp,weights_dict):
        """
        Finds the weights_dict for the threads, based on a snappy Triangulation.
        
        Parameters:
        knot_comp (snappy.Triangulation)
        
        Returns:
        dict of the form {str: vector} - gives the weights of the T-actions for the threads.
        """
        gluing_dict = QuantumAPolynomial.get_gluing_dict(knot_comp)
        thread_weights_dict = {}
        for t in range(knot_comp.num_tetrahedra()):
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
                    if local_vertex == distant_vertex:
                        pass
                        #weights = {local_vertex: 0, distant_vertex: 0}
                    else:
                        weights = {local_vertex: 1, distant_vertex: -1}

                    thread_weights_dict.update({
                        "x{0}_{1}".format(local_vertex[1:], distant_vertex[1:]) : weights   
                    })

        return thread_weights_dict


    # #### Lattice Coordinates for Threads

    def add_thread_lattice_coordinates(knot_comp,gens_dict,weights_dict):
        """ Adds the threads to the gens_dict.
        Assumes the the thread weights have already been added. Modifies gens_dict.
        
        Parameters:
        knot_comp (snappy.Triangulation()) - the triangulated knot complement
        gens_dict (dict {str: vector}) - the lattice coordinates of the generators
        weights_dict (dict {str: vector})- the weights of the generators 
        """
        num_tet = knot_comp.num_tetrahedra()
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
                try:
                    # don't add relations for this thread with itself. 
                    tmp_adjacent_edges.remove(thread)
                except ValueError:
                    # deals with threads that start/stop at a single vertex
                    pass

                #logger.debug("edges adjacent to {0} at {1}:  {2}".format(thread,vertex,tmp_adjacent_edges))
                for edge in tmp_adjacent_edges:
                    # commutation relations depend on both weights and the relative directions of approach.
                    if edge[0] == 'x':
                        #this_threads_relations.append(gens_dict[edge])
                        this_threads_relations.append(weight*gens_dict[edge])
                    elif edge[0] == 'A':
                        this_threads_relations.append(-1*gens_dict[edge])
                    else:
                        this_threads_relations.append(-weight*gens_dict[edge])
            thread_relations.append(sum(this_threads_relations))
        
        return 2*matrix(thread_relations)

    def get_relations_matrix(knot_comp,gens_dict,weights_dict):
        """Constructs the full relations matrix.
        
        Parameters:
        knot_comp (snappy.Triangulation)
        gens_dict (dict)
        weights_dict (dict)
        
        Returns:
        Matrix
        """
        num_tet = knot_comp.num_tetrahedra()
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
        
        return 1/2*omega_with_q

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
                #logger.debug("More than one option for the next generator!")
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

    def get_peripheral_curve_intersection_dict(knot_comp,curve='meridian'):
        """
        Turns the peripheral curve data from SnapPy into a dictionary,
        The values are the intersection number of the peripheral curve with the given edge.
        
        Parameters:
        M snappy.Triangulation - the triangulation of a knot complement.
        curve string - either 'meridian' or 'longitude'
        
        Returns:
        dict - of the form {'short_edge': weight}.
        """
        all_periph_data = knot_comp._get_cusp_indices_and_peripheral_curve_data()[1]
        data_start = {'meridian':0, 'longitude':2}
        
        curve_data = [all_periph_data[i] for i in range(data_start[curve],len(all_periph_data), 4)]
        curve_dict = {'a{0}{1}{2}'.format(t,v,f):curve_data[t][4*v+f]
                for t in range(knot_comp.num_tetrahedra())
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

    def get_peripheral_curve_monomial(knot_comp,gens_dict,vertices_dict,weights_dict,weights_matrix,curve='meridian'):
        """
        Gets an expression for the peripheral curve of the triangulation in terms of generators.
        There are multiple correct answers. We find this one based on SnapPy's suggestion.
        
        Parameters:
        knot_comp snappy.Triangulation - the triangulation of a knot complement.
        gens_dict Dict - the dictionary of generators.
        vertices_dict Dict - the dictionary of vertices
        weights_dict Dict - the dictionary of weights for generators.
        weights_matrix Matrix - the matrix of weights for generators
        
        curve string - either 'meridian' or 'longitude'

        returns:
        dict - the lattice representation of the curve.
        """
        
        
        intersection_dict = QuantumAPolynomial.get_peripheral_curve_intersection_dict(knot_comp,curve=curve)
        logger.debug('intersection dict for {0}: {1}'.format(curve,intersection_dict))
       ##TODO: this doesn't work when our peripheral curve interacts too much with a creased edge. 
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
        for tet in range(knot_comp.num_tetrahedra()):
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
            logger.error("The {curve} has non-zero T-weights! {weights}".format(curve=curve,weights=curve_weights))

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




    def get_T_monodromy_list(knot_comp,gens_dict):
        """
        Constructs the list of T-monodromy expressions. Assumes that gens_dict has the short edges in the right order.
        
        Parameters:
        knot_comp (snappy.Triangulation) - the triangulated knot complement
        gens_dict (dict {str:vector}) - the lattice coordinates for the generators of the skein algebra
        
        Returns:
        list of lists of strings - each element is a list of the short edges around a puncture.
        
        """
        num_tet = knot_comp.num_tetrahedra()
        # the short edges are listed in gens_dict in the right order.
        short_edge_names = list(filter(lambda name: name[0]=='a', gens_dict.keys()))
        T_monodromy_list = []
        for p in range(num_tet*4):
            T_monodromy_list.append({'qrt2':-3} | {k:1 for k in short_edge_names[3*p:3*p+3]})
        
        return T_monodromy_list

    # ### Gluing Relations




    def get_long_edge_gluing_relations_list(knot_comp,gens_dict):
        """
        Parameters:
        knot_comp (snappy.Triangulation) - the triangulated knot complement.
        gens_dict (dict {str:vector}) - the lattice coordinates for the generators of the skein algebra

        Returns
        list of dicts - gluing relations. keys are generators and values are their powers.
        """

        long_edge_gluing_relations_list = []
        # finds loops that start/end at a single vertex
        thread_finder_filter = lambda th : th[0] == 'x'
        loop_thread_filter = lambda th : th[1:].split('_')[0] != th[1:].split('_')[1]

        # threads have names starting with the letter 'x'
        threads_list = list(filter(loop_thread_filter,filter(thread_finder_filter,gens_dict.keys())))
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
                logger.debug("long_edge_gluing_relations_list:\n{}".format(long_edge_gluing_relations_list))
                return long_edge_gluing_relations_list
            else:
                # on to the next one
                th = threads_list[0]
                tmp_gluing_dict = {th:-1}
                

    def get_short_edge_gluing_relations_list(knot_comp,weights_dict):
        """
        Constructs a the list of constraints inclured by gluing together short edges.
        
        Parameters:
        knot_comp snappy.Triangulation - the triangulated knot complement.
        weights_dict - the weights of edges. Used to find where edges start/end.
        
        Returns:
        list of lists of strings - elements are lists of generator names. Should be in the right order. 
        """
        gluing_data = QuantumAPolynomial.get_gluing_dict(knot_comp)
        short_edge_gluing_relations_list = []
        qrt2_power = -4
        
        for tet in range(knot_comp.num_tetrahedra()):
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

                        # a self-identified long edge gives us two loops instead of one for the short edge relation.
                        if distant_starting_vertex == local_ending_vertex: 
                            short_edge_gluing_relations_list.append({
                                'qrt2' : qrt2_power,
                                short_edge : 1,
                                distant_short_edge : 1,
                                'x{0}_{1}'.format(local_starting_vertex,distant_ending_vertex) : 1,
                            })
                            #short_edge_gluing_relations_list.append({
                                #'qrt2' : 1,
                                #'x{0}_{1}'.format(distant_starting_vertex,local_ending_vertex) : 1
                            #})
                        elif local_starting_vertex == distant_ending_vertex:
                            short_edge_gluing_relations_list.append({
                                'qrt2' : qrt2_power,
                                short_edge : 1,
                                distant_short_edge : 1,
                                'x{0}_{1}'.format(distant_starting_vertex,local_ending_vertex) : 1
                            })
                            #short_edge_gluing_relations_list.append({
                                #'qrt2' : 1,
                                #'x{0}_{1}'.format(local_starting_vertex,distant_ending_vertex) : 1,
                            #})
                        else: # there's no self-folding to worry about here.
                            short_edge_gluing_relations_list.append({
                                'qrt2' : qrt2_power,
                                short_edge : 1,
                                'x{0}_{1}'.format(local_starting_vertex,distant_ending_vertex) : 1,
                                distant_short_edge : 1,
                                'x{0}_{1}'.format(distant_starting_vertex,local_ending_vertex) : 1
                            })
        logger.debug("short_edge_gluing_relations_list:\n{}".format(short_edge_gluing_relations_list))
        return short_edge_gluing_relations_list


    def get_thread_monodromy(self,num_factors=True,mirror_pairs=False):
        """Constructs expressions for the T-region monodromies around the 'extra handles'
        on the glued surface.
        
        Works based on the thread names. Assumes that the monodromies are each a cycle.
        if they cross themselves or something this probably won't work."""
        
        list_of_monodromies = []
        loop_thread_filter = lambda th : th[1:].split('_')[0] != th[1:].split('_')[1]
        loop_thread_finder = lambda th : th[1:].split('_')[0] == th[1:].split('_')[1]
        thread_finder_filter = lambda th : th[0] == 'x'
        
        
        #threads_list = [k for k in gens_dict.keys() if k[0] == 'x']
        threads_list = list(filter(loop_thread_filter,filter(thread_finder_filter,self.gens_dict.keys())))
        looped_threads_list = list(filter(loop_thread_finder,filter(thread_finder_filter,self.gens_dict.keys())))
        logger.debug("looped threads:\t{}".format(looped_threads_list))
        th = threads_list[0]
        tmp_monodromy_dict = {th:1}

        while threads_list:
            th_end_vertex = th[1:].split('_')[1]
            next_thread = list(filter(lambda thread : thread[1:].split('_')[0] == th_end_vertex,threads_list))[0]

            if next_thread in tmp_monodromy_dict.keys():
                # we've closed the loop!
                list_of_monodromies.append(tmp_monodromy_dict)
                threads_list = list(set(threads_list) - set(tmp_monodromy_dict.keys()))
                if threads_list == list():
                    break
                # we have another loop!
                th = threads_list[0]
                tmp_monodromy_dict = {th:1}
            else:
                tmp_monodromy_dict.update({next_thread:1})
                th = next_thread
        
        # set the q-powers. the thread monodromies come in pairs that should have opposite powers of q.
        if mirror_pairs:
            import itertools
            partnered_monodromies = []
            for m1, m2 in itertools.combinations(list_of_monodromies,2):
                # don't recheck a monodromy once we've found its partner
                if partnered_monodromies.count(m1) + partnered_monodromies.count(m2) > 0:
                    pass
                else:
                    # pick any thread in the monodromy, take the first index number
                    tester_thread = [k for k in m1.keys() if k[0] == 'x'][0].split('_')[0]
                    # swap the last two digits. the matching thread is at the other end of a long edge.
                    swapped_tester_thread = tester_thread[:-2]+tester_thread[-1]+tester_thread[-2]
                    matches = [k for k in m2.keys() if k.split('_')[0] == swapped_tester_thread]
                    if len(matches) > 0:
                        m1.update({'qrt2':-1})
                        m2.update({'qrt2': 1})
                        #logger.debug("Paired thread monodromies:{0},\t{1}".format(m1,m2))
                        #logger.debug("Their differences:\n {0}".format(QuantumAPolynomial.names_to_lattice_coordinate(m1,gens_dict)-QuantumAPolynomial.names_to_lattice_coordinate(m2,gens_dict)))
                        partnered_monodromies.append(m1)
                        partnered_monodromies.append(m2)

        if num_factors:
            for m in list_of_monodromies:
                tmp_threads = [k for k in m.keys() if k[0] == 'x']
                tmp_q_pow = m.get('qrt2',0)
                m.update({'qrt2':tmp_q_pow - len(tmp_threads)})

        # find the identified long edge groups and take the product:
        # this needs the list of long edge relations - not sure where the code should live
        # for th_m in thread_monodromy_list:
            #this_thread_le = []
            #for le in long_edge_gluing_relations_list:
                #if set(le.keys()).intersection(set(th_m.keys())) != set():
                    #this_thread_le.append(le)
            #logger.debug("The monodromy\n {0} \n with coordinate \n{1}\n has associated long edges\n {2} with total coordinate {3}".format(th_m,QuantumAPolynomial.names_to_lattice_coordinate(th_m,self.gens_dict),this_thread_le,sum([vector(QuantumAPolynomial.names_to_lattice_coordinate(le,self.gens_dict)) for le in this_thread_le])))
        for constraint in list_of_monodromies:
            if not (sum([self.gens_dict[edge] for edge in constraint])*self.weights_matrix).is_zero():
                logger.error("The follow thread monodromy constraint is not T-invariant! " +str(constraint))

        return list_of_monodromies

    def find_paired_thread_monodromies(list_of_monodromies):
        
        import itertools
        partnered_monodromies = []
        for m1, m2 in itertools.combinations(list_of_monodromies,2):
            # don't recheck a monodromy once we've found its partner
            if flatten(partnered_monodromies).count(m1) + flatten(partnered_monodromies).count(m2) > 0:
                pass
            else:
                # pick any thread in the monodromy, take the first index number
                tester_thread = [k for k in m1.keys() if k[0] == 'x'][0].split('_')[0]
                # swap the last two digits. the matching thread is at the other end of a long edge.
                swapped_tester_thread = tester_thread[:-2]+tester_thread[-1]+tester_thread[-2]
                matches = [k for k in m2.keys() if k.split('_')[0] == swapped_tester_thread]
                if len(matches) > 0:
                    partnered_monodromies.append((m1,m2))

        return partnered_monodromies


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



    def get_alg_element_from_coord(coord,alg,relations=None,double=True,q_override_value=None):
        """
        Converts from an exponent to a monomial. Assumes that 
        
        Parameters:
        coord vector -  This should be in terms of the generators of alg.
            If the (possibily doubled) length is more than the number of generators,
            then this assumes that the first coordinate is the q power.
        alg algebra - needs product and gens methods
        relations (default: None) - relations matrix for the algebra, if it's not commutative.
        double (default: True) - double the coordinate to deal with inverses
        qq_override_value (default: None) - if this is set then override the default value (alg.base().gen()) for qq

        """

        # use the base ring generator as q.
        if q_override_value == None:
            qq = alg.base().gen()
        else:
            qq = q_override_value

        # special handling for scalars
        if max(coord[1:]) == min(coord[1:]) == 0:
            return qq^coord[0]

        # find q-conversion
        if relations == None:
            conversion_q_factor = 0
        else:
            conversion_q_factor = (-1)*sum(flatten([
                [
                    coord[i]*coord[j]*relations[i,j] for i in range(j)
                ] for j in range(len(coord))
            ]))


        # deal with the case that inverses are handled as extra variables. 
        if double:
            tmp_coord = vector(flatten([ [max(coord[i],0), -min(coord[i],0)] for i in range(1,len(coord))]))
        # don't worry about inverses
        else:
            tmp_coord = vector(coord)

    
        #factors = [alg.gen(i-1)^tmp_coord[i] for i in range(1,len(tmp_coord)) if tmp_coord[i] != 0]

        monomial = reduce(alg.product,[alg.gen(i)^tmp_coord[i] for i in range(len(tmp_coord)) if tmp_coord[i] != 0])

        return qq^(coord[0]+conversion_q_factor)*monomial

    def setup_central_relations(self):
        T_monodromy_variable_names_list = QuantumAPolynomial.get_T_monodromy_list(self.knot_comp,self.gens_dict)
        #logger.debug("T_monodromy_variable_names_list:{}".format(T_monodromy_variable_names_list))
        T_monodromy_lattice_coordinate_list = [
            QuantumAPolynomial.names_to_lattice_coordinate(v,self.gens_dict) for v in T_monodromy_variable_names_list
        ]

        for coord in T_monodromy_lattice_coordinate_list:
            if not (matrix(coord)*self.weights_matrix).is_zero():
                logger.error("T-monodromy is not invariant, and it should be! "+ str(coord))



        long_edge_gluing_relations_list = QuantumAPolynomial.get_long_edge_gluing_relations_list(self.knot_comp,self.gens_dict)
        #logger.debug("long_edge_gluing_relations_list:\n"+"\n".join([str(relation) for relation in long_edge_gluing_relations_list]))

        # Check that all relations are T-invariant, and incidentally that we've listed actual threads.
        for constraint in long_edge_gluing_relations_list:
            lat_coord = QuantumAPolynomial.names_to_lattice_coordinate(constraint,self.gens_dict)
            if not (lat_coord*self.weights_matrix).is_zero():
                logger.error("A gluing relation is not T-invariant!" + str(constraint))


        short_edge_gluing_relations_list = QuantumAPolynomial.get_short_edge_gluing_relations_list(self.knot_comp,self.weights_dict)
        #logger.debug("short_edge_gluing_relations_list:\n"+"\n".join([str(relation) for relation in short_edge_gluing_relations_list]))

        # Check that all relations are T-invariant, and incidentally that we've listed actual threads.
        for constraint in short_edge_gluing_relations_list:
            if not (sum([self.gens_dict[edge] for edge in constraint])*self.weights_matrix).is_zero():
                logger.error("A gluing relation is not T-invariant, or perhaps involves non-existant threads! "+ str(constraint))
        #thread_monodromy_list = QuantumAPolynomial.get_thread_monodromy(self)
            

        self.relations_matrix = matrix([QuantumAPolynomial.names_to_lattice_coordinate(v,self.gens_dict) for v in
            short_edge_gluing_relations_list
            + long_edge_gluing_relations_list
            + T_monodromy_variable_names_list
            ])

        tmp_gluing_matrix = matrix([QuantumAPolynomial.names_to_lattice_coordinate(v,self.gens_dict) for v in
            short_edge_gluing_relations_list
            + long_edge_gluing_relations_list
            ])
        logger.debug("gluing relations matrix rank :{}, dimension: {}".format(tmp_gluing_matrix.rank(), tmp_gluing_matrix.dimensions()))

        # list of relations



        self.short_edge_gluing_relations_list = short_edge_gluing_relations_list
        self.long_edge_gluing_relations_list = long_edge_gluing_relations_list
        self.T_monodromy_variable_names_list = T_monodromy_variable_names_list


    def build_g_algebra(self,quotient_omega,q_override_value=None):
        S = FractionField(QQ['qrt2'])
        if q_override_value != None:
            qrt2 = q_override_value
        else:
            qrt2 = S.gen()

        ngens = self.quotient_lattice.ngens()

        doubled_quotient_omega = (quotient_omega[1:,1:]).tensor_product(Matrix([[1,-1],[-1,1]]))

        #A = FreeAlgebra(S,2*ngens-1,['M','L']+['x{}'.format(i) for i in range(ngens-3)]+['Mi','Li']+['xi{}'.format(i) for i in range(ngens-3)])
        num_skel = self.num_tet -1
        num_thread_mon = self.num_tet - 1
        logger.debug("Quotient ngens: {0}, Skeletal variables: {1}, Thread monodromies: {2}".format(ngens,num_skel,num_thread_mon))
        if num_skel + num_thread_mon + 3 != ngens:
            logger.error("The generators are not all accounted for in the g algebra! ngens: {0}, we're using: {1}.".format(ngens,num_skel+num_thread_mon+3))

        A = FreeAlgebra(S,2*ngens-2,
                        ['M','Mi','L','Li']
                        +flatten([['r{}'.format(i+1),'ri{}'.format(i+1)] for i in range(num_thread_mon)])
                        +flatten([['Y{}'.format(i+1),'Yi{}'.format(i+1)] for i in range(num_skel)])
            )
        A_gens_dict  = {A.gen(i) : [0]*(i) + [1] + [0]*(A.ngens()-i-1) for i in range(A.ngens())}

        import itertools
        relations = {}
        for g2,g1 in itertools.combinations(A_gens_dict.keys(),2):
            power = 2*(matrix(A_gens_dict[g2])*doubled_quotient_omega*(matrix(A_gens_dict[g1]).transpose()))[0][0]
            if power != 0:
                relations[g1*g2] = qrt2**int(-power)*g2*g1

        self.qrt2 = qrt2
        return A.g_algebra(relations)



    def get_bulk_relations_in_g_algebra(self,q_override_value=None):

        if q_override_value != None:
            qrt2 = q_override_value
        else:
            qrt2 = self.H.base().gen()

        generic_bulk_relation = [ # this was found by hand. its morally 0 = q^(3/2) A13 A02 + q^(3/2) A03 A12 + q^(5/2) A01 A23, the last +- 2 is from multiplying through by A13 A02
            {'qrt2':4+1-2, 'A{t}13':-1, 'A{t}02':-1, 'A{t}03':1, 'a{t}32':1, 'a{t}01':1, 'A{t}12':1, 'a{t}23':1, 'a{t}10':1},
            {'qrt2':4-1+2, 'A{t}13':-1, 'A{t}02':-1, 'A{t}01':1, 'a{t}03':-1,'a{t}12':-1, 'A{t}23':1, 'a{t}21':-1, 'a{t}30':-1},
            {'qrt2':3}
        ]

        self.bulk_relations = [
                [ QuantumAPolynomial.names_to_lattice_coordinate({k.format(t=tt) : v for k,v in monomial.items()},self.gens_dict) for monomial in generic_bulk_relation ]
                for tt in range(self.knot_comp.num_tetrahedra())
        ]
         

        # Checks -
        if logger.level <= 10: # check commutation relations
            bulk_lhs = [ # this equals A13*A02 and was found by hand
                {'A{t}03':1, 'a{t}32':1, 'a{t}01':1, 'A{t}12':1, 'a{t}23':1, 'a{t}10':1},
                {'A{t}01':1, 'a{t}03':-1,'a{t}12':-1, 'A{t}23':1, 'a{t}21':-1, 'a{t}30':-1}
            ]
            bulk_rhs = {'A{t}13':-1, 'A{t}02':-1}

            for tt in range(self.knot_comp.num_tetrahedra()):
                this_bulk_lhs = [QuantumAPolynomial.names_to_lattice_coordinate({k.format(t=tt) : v for k,v in monomial.items()},self.gens_dict) for monomial in bulk_lhs]
                this_bulk_rhs = QuantumAPolynomial.names_to_lattice_coordinate({k.format(t=tt) : v for k,v in bulk_rhs.items()},self.gens_dict)
                logger.debug("q-scaling factors for tetrahedra {0}: {1}, {2}".format(tt, -this_bulk_rhs*self.omega_with_q*matrix(this_bulk_lhs[0]).transpose(), -this_bulk_rhs*self.omega_with_q*matrix(this_bulk_lhs[1]).transpose()))


        g_algebra_crossing_relations = [
                reduce(self.H.summation,[
                           QuantumAPolynomial.get_alg_element_from_coord(self.quotient_basis.solve_left(vector(self.pi(mon))),self.H,relations=self.quotient_omega,double=True)
                           for mon in relation ]
                       ) for relation in self.bulk_relations
                ]

        return g_algebra_crossing_relations

    def restrict_to_subalgebra(elem,alg):
        coeff = lambda m : alg.base()(str(elem.monomial_coefficient(m)))
        restricted_monomials = [
                coeff(m) * alg(alg.free_algebra()(str(m))) for m in elem.monomials()
                ]

        return reduce(alg.summation,restricted_monomials)


    def compute_skein_module(self):

        # find the weights matrix and T-invariant lattice
        self.weights_matrix = QuantumAPolynomial.get_weights_matrix(self.vertices_dict,self.weights_dict)



        ### Build the central quotient.
        
        # Get the central relations
        QuantumAPolynomial.setup_central_relations(self)
        to_vec = lambda name : QuantumAPolynomial.names_to_lattice_coordinate(name,self.gens_dict)
        T_monodromy_lattice_coordinate_list = [to_vec(name) for name in self.T_monodromy_variable_names_list]
        gluing_relation_coordinate_list = [to_vec(name) for name in self.short_edge_gluing_relations_list + self.long_edge_gluing_relations_list]

        self.thread_monodromy_list = QuantumAPolynomial.get_thread_monodromy(self)
        self.paired_thread_monodromies = QuantumAPolynomial.find_paired_thread_monodromies(self.thread_monodromy_list)

        logger.info('TEX thread monodromies: {}'.format([
            'r_{0} = {1}'.format(ii+1,QuantumAPolynomial.names_to_str(flatten(self.paired_thread_monodromies)[ii],self.omega_with_q,self.gens_dict))
            for ii in range(len(self.thread_monodromy_list))
            ]))

        for m1,m2 in self.paired_thread_monodromies:
            if self.omega_with_q*(vector(to_vec(m1))-vector(to_vec(m2))) != 0:
                logger.warning("Thread monodromy difference isn't central: {}")

        self.thread_monomials = matrix([
            to_vec(name) for name in flatten(self.paired_thread_monodromies)
            ])

        logger.debug("thread_monodromy_list. These should come in pairs of equal length. \n"+"\n".join([str(relation) for relation in self.thread_monodromy_list]))


        self.monomial_relations = block_matrix([
                [matrix([to_vec(name) for name in self.thread_monodromy_list])],
                [matrix([to_vec(name) for name in self.short_edge_gluing_relations_list])],
                [matrix([to_vec(name) for name in self.long_edge_gluing_relations_list])],
                [matrix([to_vec(name) for name in self.T_monodromy_variable_names_list])]
                ])

        if logger.level <= 10:
            logger.debug("monomial relation q-values:\n{}".format(self.monomial_relations.columns()[0]))

            logger.debug("Full lattice rank: {0}".format(len(self.gens_dict)))

        self.central_relations = block_matrix([
                [matrix([to_vec(name) for name in self.short_edge_gluing_relations_list])],
                [matrix([to_vec(name) for name in self.long_edge_gluing_relations_list])],
                [matrix([to_vec(name) for name in self.T_monodromy_variable_names_list])]
                #[matrix([to_vec(m1)-to_vec(m2) for m1,m2 in self.paired_thread_monodromies])]
                ])

        
        # Get the meridian and longitude

        unordered_meridian = QuantumAPolynomial.get_peripheral_curve_monomial(self.knot_comp,self.gens_dict,self.vertices_dict,self.weights_dict,self.weights_matrix,curve='meridian')
        q_meridian_scaling = sum([-abs(v) for v in unordered_meridian.values()]) + unordered_meridian.get('qrt2',0)
        unordered_meridian['qrt2'] = q_meridian_scaling
        self.meridian =QuantumAPolynomial.names_to_lattice_coordinate(unordered_meridian,self.gens_dict)
        logger.debug("Ordered Meridian: {}".format(self.meridian))
        logger.info("TEX Meridian : {}".format(QuantumAPolynomial.names_to_str(
            QuantumAPolynomial.lattice_coord_to_dict(self.meridian,self.gens_dict),self.omega_with_q,self.gens_dict
            )))

        unordered_longitude = QuantumAPolynomial.get_peripheral_curve_monomial(self.knot_comp,self.gens_dict,self.vertices_dict,self.weights_dict,self.weights_matrix,curve='longitude')
        q_longitude_scaling = sum([-abs(v) for v in unordered_longitude.values()]) + unordered_longitude.get('qrt2',0)
        unordered_longitude['qrt2'] = q_longitude_scaling
        self.longitude = QuantumAPolynomial.names_to_lattice_coordinate(unordered_longitude,self.gens_dict)
        logger.debug("Ordered Longitude: {}".format(self.longitude))
        logger.info("TEX Longitude: {}".format(QuantumAPolynomial.names_to_str(
            QuantumAPolynomial.lattice_coord_to_dict(self.longitude,self.gens_dict),self.omega_with_q,self.gens_dict
            )))

        ## Choose a custom basis for the invariant sublattice

        # remove the q-power before putting vectors in the quotient.
        strip_q = lambda coord: vector([0] + list(coord[1:]))
        unordered_invariant_sublattice = self.weights_matrix.left_kernel()
        self.invariant_center = (unordered_invariant_sublattice.basis_matrix()*self.omega_with_q).right_kernel().intersection(unordered_invariant_sublattice)
        #logger.debug("invariant lattice: {}".format(unordered_invariant_sublattice))
        #logger.debug("self.invariant_center: {}".format(self.invariant_center))

        new_basis = [self.gens_dict['qrt2'],strip_q(self.meridian),strip_q(self.longitude)] 

        # first thread monodromy, then the other monomial relations, then whatever's needed to fill out the lattice.
        inv_rank = unordered_invariant_sublattice.rank()
        while span(new_basis).rank() < inv_rank:
            for vec in self.thread_monomials.rows() + T_monodromy_lattice_coordinate_list + gluing_relation_coordinate_list + self.invariant_center.basis() + unordered_invariant_sublattice.basis():
            #for vec in self.thread_monomials.rows() + unordered_invariant_sublattice.basis():
                if strip_q(vec) not in span(new_basis):
                    new_basis.append(strip_q(vec))

        logger.debug("invariant lattice basis: {}".format([QuantumAPolynomial.lattice_coord_to_dict(ee,self.gens_dict) for ee in new_basis]))
        for vec in new_basis:
            if matrix(vec)*self.weights_matrix !=0:
                logger.warning("Non-invariant basis element: {}".format(QuantumAPolynomial.lattice_coord_to_dict(vec,self.gens_dict)))

        #logger.debug("T-monodromies: {}".format(self.T_monodromy_variable_names_list))
        self.invariant_sublattice = unordered_invariant_sublattice.submodule_with_basis(new_basis)
        logger.debug("T-invariant lattice is rank: {0}".format(self.invariant_sublattice.rank()))


        logger.debug("central_relations rank: {}".format(self.central_relations.rank()))

        # Take the quotient!
        self.center_of_relations_with_q = self.invariant_center.intersection(self.monomial_relations.row_space())
        self.center_of_relations = self.center_of_relations_with_q.submodule_with_basis([
            g for g in self.center_of_relations_with_q.gens() if g != self.gens_dict['qrt2']
            ])

        logger.debug("center_of_relations: {}".format(self.center_of_relations))

        self.quotient_lattice = self.invariant_sublattice.quotient(self.center_of_relations)
        self.pi = self.quotient_lattice.coerce_map_from(self.quotient_lattice.V())
        logger.debug("quotient: {}".format(self.quotient_lattice))
        #logger.debug("intersection's quotient by central relations {}".format(self.center_of_relations.quotient(self.central_relations.row_space()   )))
        #logger.debug("center's quotient by central relations: {}".format(self.invariant_center.quotient(self.invariant_center.intersection(self.central_relations.row_space()))))
        #logger.debug("center's quotient by center of all monomial relations: {}".format(self.invariant_center.quotient(self.center_of_relations)))



        # checking in on torsion - we don't expect to have any, I think.
        self.invariants = self.quotient_lattice.invariants()
        for g in self.quotient_lattice.gens():
            if g.additive_order() < sage.rings.infinity.PlusInfinity():
                logger.warning("Torison element in the quotient! {}".format(QuantumAPolynomial.lattice_coord_to_dict(g.lift(),self.gens_dict)))


        logger.debug("In the quotient: \n qrt2 = {0}\n M = {1} \n L = {2}".format(self.pi(self.gens_dict['qrt2']), self.pi(self.meridian), self.pi(self.longitude)))
        logger.debug("thread monomials in quotient: {}".format([self.pi(v) for v in self.thread_monomials.rows()]))
        lift_and_name = lambda v : QuantumAPolynomial.lattice_coord_to_dict(v.lift(),self.gens_dict)
        long_edge_filter = lambda elem : reduce(lambda x,y: x or y, map(lambda k : k.count('A') >0, lift_and_name(elem).keys()))
        logger.info("TEX Skeletal variables {}".format([QuantumAPolynomial.names_to_str(lift_and_name(g),self.omega_with_q,self.gens_dict) for g in self.quotient_lattice.gens() if long_edge_filter(g)]))

        logger.info("TEX All Generators of quotient: {}".format("\n"+"\n".join([QuantumAPolynomial.names_to_str(lift_and_name(g),self.omega_with_q,self.gens_dict) for g in self.quotient_lattice.gens()])))
        # make the relations matrix for the quotient
        self.quotient_omega = Matrix([
            [(Matrix(v.lift())*self.omega_with_q*Matrix(w.lift()).transpose())[0,0] for v in self.quotient_lattice.gens()]
            for w in self.quotient_lattice.gens()
            ])


        #################### Make the g algebra! ###########################
        self.H = QuantumAPolynomial.build_g_algebra(self,self.quotient_omega,q_override_value=None)
        # make qrt2 a usable value. Note that self.qrt2 was already defined and is the specialization
        self.H.base().inject_variables()
        logger.info("big g_algebra: {}".format(self.H))
        
        alg_gens = [QuantumAPolynomial.get_alg_element_from_coord(vector(g),self.H,relations=self.quotient_omega) for g in self.quotient_lattice.gens()]
        logger.debug("lattice generators in g algebra: {}".format(alg_gens))
        alg_inv_gens = [QuantumAPolynomial.get_alg_element_from_coord(-vector(g),self.H,relations=self.quotient_omega) for g in self.quotient_lattice.gens()]
        logger.debug("inverse lattice generators in g algebra: {}".format(alg_inv_gens))
    
        self.inverse_relations = [alg_gens[1:][i]*alg_inv_gens[1:][i] -1 for i in range(self.quotient_lattice.ngens()-1)]
        logger.debug("inverse_relations: {}".format(self.inverse_relations))

        self.thread_monomials_in_g_alg = [QuantumAPolynomial.get_alg_element_from_coord(vector(self.pi(mon)),self.H,relations=self.quotient_omega) for mon in self.thread_monomials.rows() ]
        logger.debug("thread monomials in g algebra: {}".format(self.thread_monomials_in_g_alg))
        tmp_thread_variables_list = []
        
        # recreates the behavior of free_algebra's variable function.
        g_alg_mon_variables = lambda mon : [mon.parent()(g) for g in mon.parent().free_algebra()(str(mon.lm())).variables()]
        for mon in self.thread_monomials.rows():
            tmp_thread_variables_list += g_alg_mon_variables(QuantumAPolynomial.get_alg_element_from_coord(vector(self.pi(mon)),self.H,relations=self.quotient_omega))

        self.thread_monodromy_variables = set(tmp_thread_variables_list)
        logger.debug("thread_monodromy_variables: {}".format(self.thread_monodromy_variables))
        self.skeletal_gens = [g for g in self.H.gens() if str(g)[0] == 'Y']


        self.quotient_basis = matrix([vector(g) for g in self.quotient_lattice.gens()]) 
        #logger.debug("quotient basis:\n{}".format(self.quotient_basis))

        self.bulk_relations_in_alg = QuantumAPolynomial.get_bulk_relations_in_g_algebra(self)
        logger.info("bulk relations: {}".format(self.bulk_relations_in_alg))


        if self.finish_variable_elimination == True:
            logger.info("starting variable elimination!")
            ## Eliminate the skeletal variables, then make a g_algebra out of what's left.
            self.partial_elimination_ideal = self.H.ideal(self.inverse_relations+self.bulk_relations_in_alg).elimination_ideal(self.skeletal_gens)
            logger.debug("eliminated {0} : {1}".format(self.skeletal_gens,self.partial_elimination_ideal.gens()))

            tmp_variables_list = []
            for g in self.partial_elimination_ideal.gens():
                for mon in g.monomials():
                    tmp_variables_list += g_alg_mon_variables(mon)

            self.partial_elimination_variables = list(set(tmp_variables_list))
            self.partial_elimination_variables.sort(reverse=True)
            logger.debug("variables after partial elimination: {}".format(self.partial_elimination_variables))

            restricted_F = FreeAlgebra(self.H.base(),[str(g) for g in self.partial_elimination_variables])

            # build the relations dictionary for the subalgebra. 
            restricted_relations = {}
            tmp_F = self.H.free_algebra()
            tmp_relations = self.H.relations()
            import itertools
            for g2,g1 in itertools.combinations(restricted_F.gens(),2):
                tmp_RHS = tmp_relations.get(tmp_F(str(g1*g2)),None)
                if tmp_RHS == None:
                    pass
                else:
                    new_RHS = restricted_F.base()(str(tmp_RHS.lc())) * restricted_F(str(tmp_RHS.lm()))
                    restricted_relations[g1*g2] = new_RHS

            logger.debug("restricted_relations: {}".format(restricted_relations))

            self.restricted_algebra = restricted_F.g_algebra(restricted_relations)
            logger.debug("restricted algebra: {}".format(self.restricted_algebra))

            # find the right factors to set the thread relations to:
            thread_subs_dict = {}
            for g in self.thread_monomials_in_g_alg:
                thread_subs_dict[QuantumAPolynomial.restrict_to_subalgebra(g.lm(),self.restricted_algebra)] = self.restricted_algebra.base()(str(g.lc()))^-1

            logger.debug("thread subs dict: {}".format(thread_subs_dict))

            # import the partial elimination ideal into the restricted algebra:
            self.restricted_ideal = self.restricted_algebra.ideal([
                QuantumAPolynomial.restrict_to_subalgebra(g,self.restricted_algebra).subs(thread_subs_dict)
                for g in self.partial_elimination_ideal.gens()
                ])


            self.quantum_A_poly_ideal = self.restricted_ideal.elimination_ideal([self.restricted_algebra(restricted_F('Mi')), self.restricted_algebra(restricted_F('Li'))])
            logger.info("Quantum A-polynomial is in this ideal:\n {}\n".format(self.quantum_A_poly_ideal))
            logger.info("Compare with the literature classical A-polynomial: {}".format(self.ref_A_poly))
            logger.debug("Quantum A-poly ideal has {0} generators.".format(self.quantum_A_poly_ideal.ngens()))
        else:
            logger.info("Skipping final eliminations step.")

