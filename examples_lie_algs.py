from sage.all import QQ, LieAlgebra, libgap, FunctionField


def lie_algebra_upper_triangular_matrices(n, strict=False, field=QQ):

    # the dimension of the algebra
    dimL = n*(n-1)/2 if strict else n*(n+1)//2

    # this function calculates the product of two elementary matrices.
    # the first matrix has 1 at position i, j; the second has 1 at position k,l
    # this lambda function returns the position where the product has one;
    # of no such position exists, then the output is ()
    def prod_func(i, j, k, m):
        return (j == k)*(i, m)

    # the standard basis is labelled as xij where i, j indicates the
    # non-zero entry
    basis_labels = [(i+1, j+1) for i in range(n)
                    for j in range(i+1 if strict else i, n)]

    # we sort the basis labels so that they correspond to the derived series
    # of L this helps in the invariant field computation
    basis_labels.sort(key=lambda x: x[1]-x[0], reverse=True)

    # prod_dict is the dictionary for the bracket on L
    prod_dict = {}
    # this is just the list of strings to name the basis elements of L
    basis_names = ['x'+str(x[0])+str(x[1]) for x in basis_labels]

    # we run i, j over indices 0<= i < j <= dimL-1
    for i in range(dimL):
        # get the basis label for the first basis element b1
        x, y = basis_labels[i]
        for j in range(i+1, dimL):
            # basis label for the second basis element b2
            u, v = basis_labels[j]
            # we calculate b1*b2 and b2*b1 using the lambda expression above
            pr1, pr2 = prod_func(x, y, u, v), prod_func(u, v, x, y)

            # we calculate the entry to be added to the product dictionary
            prod_dict_entry = {}

            if pr1 != ():
                prod_dict_entry['x'+str(pr1[0])+str(pr1[1])] = 1
            if pr2 != ():
                prod_dict_entry['x'+str(pr2[0])+str(pr2[1])] = -1

            # add the entry to the prod dictionary
            if prod_dict_entry != {}:
                prod_dict[('x'+str(x)+str(y), 'x'+str(u)+str(v))
                          ] = prod_dict_entry

    # return the Lie algebra
    return LieAlgebra(field, prod_dict, basis_names)


def standard_filiform_lie_algebra(n, field=QQ):

    vars = ['y'+str(k) for k in range(n-1)] + ['x']
    rel_dict = {('x', 'y'+str(k)): {'y'+str(k-1): 1} for k in range(1, n-1)}
    return LieAlgebra(field, rel_dict, vars)


def lie_alg_from_libgap_to_sage(lie_alg):
    '''
        INPUT:

            - a Lie algebra written in Gap via libgap

        OUTPUT:

            - a Lie algebra written in Sage.

        COMMENT:

            This function converts a Lie algebra, written in the Gap language,
            into a Lie algebra written in the Sage language.

        EXAMPLES:

            sage: l = libgap.NilpotentLieAlgebra(QQ, [6,14])

            sage: L = lie_alg_from_gap_to_sage(l)

            sage: L

            Lie algebra on 6 generators (x0, x1, x2, x3, x4, x5) over Rational
            Field

    '''
    dimL = libgap.Dimension(lie_alg).sage()
    bL = lie_alg.Basis()
    var = [f"x{i}" for i in range(dimL)]  # ['x0', ..., 'xn']
    dicio = {}  # dicionário para definir nova álgebra de Lie
    # escrever o dicionário dicio {('xi','xj'): {'x0':num0, 'x1':num1, ...},
    # ... }
    for i in range(dimL):
        for j in range(i, dimL):
            dicioAux = {}
            v = bL.Coefficients(bL[i]*bL[j])
            # escrever [xi,xj] na base especial
            for k in range(dimL):
                dicioAux[var[k]] = v[k]
            dicio[(var[i], var[j])] = dicioAux
    return LieAlgebra(QQ, dicio, names=var)


def nilpotent_lie_algebra(F, args, standard_basis=False):
    '''
        INPUT:

            - a field, F
            - a list with arguments, args
            - an optional Boolean value

        OUTPUT:

            - a nilpotent Lie algebra

        COMMENT:

            This function uses the gap.NilpotentLieAlgebra function to return
            a nilpotent Lie algebra in Gap. After that, it uses the
            lie_alg_from_gap_to_sage function to convert that algebra, written
            in the Gap language, to an algebra in the Sage language. If the
            optional value of the function is True, then the function changes
            the basis of the algebra to the basis given by the
            triangular_basis_lie_algebra function.

        EXAMPLES:

            sage: L = nilpotent_lie_algebra(QQ, [6,16])

            sage: L

            Lie algebra on 6 generators (x0, x1, x2, x3, x4, x5) over Rational
            Field

    '''
    if not hasattr(libgap, 'NilpotentLieAlgebra'):
        libgap.LoadPackage('liealgdb')
    L = lie_alg_from_libgap_to_sage(libgap.NilpotentLieAlgebra(F, args))
    return L


def lie_alg_family_6_19(F):

    K = FunctionField(F, 'ε')
    ε = K.gens()[0]
    mult_table = {('x2', 'x5'): {'x1': -1},
                  ('x3', 'x4'): {'x1': -1},
                  ('x3', 'x6'): {'x1': -ε},
                  ('x4', 'x5'): {'x2': 1},
                  ('x4', 'x6'): {'x3': 1}}

    return LieAlgebra(K, mult_table, names='x1,x2,x3,x4,x5,x6')


def lie_alg_family_6_21(F):

    K = FunctionField(F, 'ε')
    ε = K.gens()[0]
    mult_table = {('x2', 'x5'): {'x1': -1},
                  ('x3', 'x6'): {'x1': -ε},
                  ('x4', 'x5'): {'x2': -1},
                  ('x4', 'x6'): {'x3': -1},
                  ('x5', 'x6'): {'x4': 1}}

    return LieAlgebra(K, mult_table, names='x1,x2,x3,x4,x5,x6')


def lie_alg_family_6_22(F):

    K = FunctionField(F, 'ε')
    ε = K.gens()[0]
    mult_table = {('x3', 'x4'): {'x1': 1},
                  ('x3', 'x5'): {'x2': 1},
                  ('x4', 'x6'): {'x2': ε},
                  ('x5', 'x6'): {'x1': 1}}

    return LieAlgebra(K, mult_table, names='x1,x2,x3,x4,x5,x6')


def lie_alg_family_6_24(F):

    K = FunctionField(F, 'ε')
    ε = K.gens()[0]
    mult_table = {('x3', 'x5'): {'x1': -1},
                  ('x3', 'x6'): {'x2': -1},
                  ('x4', 'x5'): {'x2': -ε},
                  ('x4', 'x6'): {'x1': -1},
                  ('x5', 'x6'): {'x3': 1}}

    return LieAlgebra(K, mult_table, names='x1,x2,x3,x4,x5,x6')
