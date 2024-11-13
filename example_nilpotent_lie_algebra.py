#-------------
def lie_algebra_strict_upper_triangular_matrices(n):
    dimL = (n*(n-1)/2).numerator()
    basis_matrices = [[]]*n
    for i in range(n):
        for j in range(n):
            if j <= i:
                basis_matrices[i] = basis_matrices[i] + [0]
            else:
                elem_matrix = matrix(QQ, n)
                elem_matrix[i,j] = 1
                basis_matrices[i] = basis_matrices[i] + [elem_matrix]
    basis_lie_alg = [[]]*n
    for i in range(n):
        for j in range(n):
            if j <= i:
                basis_lie_alg[i] = basis_lie_alg[i] + [0]
            else:
                basis_lie_alg[i] = basis_lie_alg[i] + ["E" + str(i+1) + str(j+1)]
    dict_matr_lie = {}
    for i in range(n):
        for j in range(n):
            dict_matr_lie[basis_lie_alg[i][j]] = basis_matrices[i][j]
    keys = [0]*dimL
    aux_number = 0
    for i in range(n):
        for j in range(n):
            if basis_lie_alg[i][j] != 0:
                keys[aux_number] = basis_lie_alg[i][j]
                aux_number = aux_number + 1
    value = [0]*dimL
    aux_number = 0
    for i in range(n):
        for j in range(n):
            if basis_matrices[i][j] != 0:
                value[aux_number] = basis_matrices[i][j]
                aux_number = aux_number + 1
    dict_lie = {}
    for i in range(dimL):
        for j in range(i,dimL):
            if get_keys_from_value(dict_matr_lie,value[i]*value[j] - value[j]*value[i])[0] != 0:
                dict_lie[(keys[i], keys[j])] = {get_keys_from_value(dict_matr_lie,value[i]*value[j] - value[j]*value[i])[0]:1}
    #return dict_lie
    L = LieAlgebra(QQ, dict_lie, names = keys)
    return L
#-------------

#-------------
# the solvable Lie algebra of upper triangular matrices
def lie_algebra_upper_triangular_matrices( n ):

    # the dimension of the algebra
    dimL = n*(n+1)//2

    # this function calculates the product of two elementary matrices.
    # the first matrix has 1 at position i, j; the second has 1 at position k,l
    # this lambda function returns the position where the product has one;
    # of no such position exists, then the output is ()
    prod_func = lambda i,j,k,l: (j==k)*(i,l)

    # the standard basis is labelled as xij where i, j indicates the non-zero entry
    basis_labels = [(i+1,j+1) for i in range( n ) for j in range( i, n )]

    # we sort the basis labels so that they correspond to the derived series of L
    # this helps in the invariant field computation
    basis_labels.sort( key = lambda l : l[1]-l[0], reverse = True  )

    # prod_dict is the dictionary for the bracket on L
    prod_dict = {}
    # this is just the list of strings to name the basis elements of L
    basis_names = [ 'x'+str(l[0])+str(l[1]) for l in basis_labels ]

    # we run i, j over indices 0<= i < j <= dimL-1
    for i in range( dimL ):
        # get the basis label for the first basis element b1
        x, y = basis_labels[i]
        for j in range( i+1, dimL ):
            # basis label for the second basis element b2
            u, v = basis_labels[j]
            # we calculate b1*b2 and b2*b1 using the lambda expression above
            pr1, pr2 = prod_func( x, y, u, v ), prod_func( u, v, x, y )

            # we calculate the entry to be added to the product dictionary
            prod_dict_entry = {}

            if pr1 != ():
                prod_dict_entry['x'+str(pr1[0])+str(pr1[1])] = 1
            if pr2 != ():
                prod_dict_entry['x'+str(pr2[0])+str(pr2[1])] = -1

            # add the entry to the prod dictionary
            if prod_dict_entry != {}:
                prod_dict[('x'+str(x)+str(y), 'x'+str(u)+str(v))] = prod_dict_entry

    # return the Lie algebra
    return LieAlgebra(QQ, prod_dict, basis_names )
#-------------

#-------------
def standard_filiform_lie_algebra( n ):

    vars = [ 'y'+str(k) for k in range( n-1 ) ] + ['x']
    rel_dict = { ('x','y'+str(k)): {'y'+str(k-1): 1} for k in range(1,n-1) }
    return LieAlgebra( QQ, rel_dict, vars )
#-------------
