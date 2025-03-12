"""
This module provides functions to create and manipulate various Lie algebras.
It includes functions to generate upper triangular matrices, standard filiform
Lie algebras, and specific families of Lie algebras. Additionally, it provides
utilities to convert Lie algebras from the GAP system to SageMath.

Functions:
- lie_algebra_upper_triangular_matrices: Generates a Lie algebra of upper
  triangular matrices.
- standard_filiform_lie_algebra: Generates a standard filiform Lie algebra.
- lie_alg_from_libgap_to_sage: Converts a Lie algebra from GAP to SageMath.
- nilpotent_lie_algebra: Generates a nilpotent Lie algebra using GAP and
  converts it to SageMath.
- lie_alg_family_6_19: Generates a specific family of Lie algebra.
- lie_alg_family_6_21: Generates a specific family of Lie algebra.
- lie_alg_family_6_22: Generates a specific family of Lie algebra.
- lie_alg_family_6_24: Generates a specific family of Lie algebra.

The Lie algebras generated  by the last four functions are the Lie algebras with 
dimension 6 and index 19, 21, 22, 24 in the paper by Cicalo-de Graaf-Schneder
"6-dimensional nilpotent Lie algebras."
"""

from sage.all import QQ, LieAlgebra, libgap, FunctionField


def lie_algebra_upper_triangular_matrices(n, strict=False, field=QQ):
    """
    Generates a Lie algebra of upper triangular matrices.

    Args:
        n (int): The size of the matrices.
        strict (bool, optional): If True, generates strictly upper triangular 
        matrices. Defaults to False.
        field (Field, optional): The field over which the Lie algebra is 
        defined. Defaults to QQ.

    Returns:
        LieAlgebra: The Lie algebra of upper triangular matrices.

    Examples:
        sage: L = lie_algebra_upper_triangular_matrices(3)
        sage: L
        Lie algebra on 6 generators (x12, x13, x23, x11, x22, x33) over Rational Field
    """
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
            if prod_dict_entry:
                prod_dict[('x'+str(x)+str(y), 'x'+str(u)+str(v))
                          ] = prod_dict_entry

    # return the Lie algebra
    return LieAlgebra(field, prod_dict, basis_names)


def standard_filiform_lie_algebra(n, field=QQ):
    """
    Generates a standard filiform Lie algebra.

    Args:
        n (int): The dimension of the Lie algebra.
        field (Field, optional): The field over which the Lie algebra is
        defined. Defaults to QQ.

    Returns:
        LieAlgebra: The standard filiform Lie algebra.

    Examples:
        sage: L = standard_filiform_lie_algebra(5)
        sage: L
        Lie algebra on 5 generators (y0, y1, y2, y3, x) over Rational Field
    """
    bas_strings = ['y'+str(k) for k in range(n-1)] + ['x']
    rel_dict = {('x', 'y'+str(k)): {'y'+str(k-1): 1} for k in range(1, n-1)}
    return LieAlgebra(field, rel_dict, bas_strings)


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


def nilpotent_lie_algebra(F, args):
    """
    Generates a nilpotent Lie algebra using GAP and converts it to SageMath.

    Args:
        F (Field): The field over which the Lie algebra is defined.
        args (list): A list of arguments for the GAP NilpotentLieAlgebra function.
        change_basis (bool, optional): If True, changes the basis of the algebra 
        to the basis given by the triangular_basis_lie_algebra function. Defaults to False.

    Returns:
        LieAlgebra: The nilpotent Lie algebra.

    Examples:
        sage: L = nilpotent_lie_algebra(QQ, [6,16])
        sage: L
        Lie algebra on 6 generators (x0, x1, x2, x3, x4, x5) over Rational Field
    """
    if not hasattr(libgap, 'NilpotentLieAlgebra'):
        libgap.LoadPackage('liealgdb')
    L = lie_alg_from_libgap_to_sage(libgap.NilpotentLieAlgebra(F, args))
    return L


def lie_alg_family_6_19(F):
    """
    Generates a specific family of Lie algebra with dimension 6 and index 19.

    Args:
        F (Field): The field over which the Lie algebra is defined.

    Returns:
        LieAlgebra: The Lie algebra with dimension 6 and index 19.

    Examples:
        sage: L = lie_alg_family_6_19(QQ)
        sage: L
        Lie algebra on 6 generators (x1, x2, x3, x4, x5, x6) over Rational Field
    """
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
