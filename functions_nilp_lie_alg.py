#------------- 
def lie_alg_from_gap_to_sage(l):
    '''
        INPUT:
        
            - a Lie algebra written in Gap.
            
        OUTPUT:
            
            - a Lie algebra written in Sage.
        
        COMMENT:
        
            This function converts a Lie algebra, written in the Gap language, into a Lie algebra written in the Sage language.
        
        EXAMPLES:
        
            sage: l = gap.NilpotentLieAlgebra(QQ, [6,14])

            sage: L = lie_alg_from_gap_to_sage(l)

            sage: L

            Lie algebra on 6 generators (x0, x1, x2, x3, x4, x5) over Rational Field
        
    '''
    dimL = int(gap.Dimension(l))
    P = PolynomialRing(QQ, dimL, "x")
    bL = gap.Basis(l)
    x = list(P.gens())
    xstr = [str(i) for i in x] # ['x0', ..., 'xn']
    var = ', '.join(str(i) for i in x) # 'x0, ..., xn'
    dicio = {} # dicionário para definir nova álgebra de Lie
    # escrever o dicionário dicio {('xi','xj'): {'x0':num0, 'x1':num1, ...}, ... }
    for i in range(dimL):
        for j in range(i, dimL):
            dicioAux = {}
            v = gap.Coefficients(bL, bL[i + 1]*bL[j + 1]) # escrever [xi,xj] na base especial
            for k in range(dimL):
                dicioAux[xstr[k]] = v[k+1]
            dicio[(xstr[i], xstr[j])] = dicioAux
    L = LieAlgebra(QQ, dicio, names = var)
    return L
#-------------

#-------------
def satisfies_jacobi(L):
    '''
        INPUT:
        
            - a Lie algebra constructed from a dictionary.
            
        OUTPUT:
            
            - `True` if the identity is satisfied or list with an example of three elements where the identity fails.
        
        COMMENT:
        
            This function tests the Jacobi identity for the base elements of an algebra and returns `True` if the identity is satisfied for all elements. Otherwise, it returns a list with an example of three elements where the identity fails.
        
        EXAMPLES:
        
            sage: dicio = { ('x0','x1'):{'x3':1}, ('x0','x2'):{'x4':1}, ('x0','x3'):{'x5':1}, ('x0','x4'):{'x5':1}, ('x1','x2'):{'x5':1}, ('x1','x4'):{'x6':1}, ('x2','x3'):{'x6':1} }

            sage: L = LieAlgebra(QQ, dicio, names='x0,x1,x2,x3,x4,x5,x6')

            sage: satisfies_jacobi(L)

            True
        
    '''
    bL = list(L.basis())
    dimL = len(bL)
    for i in range(dimL):
        for j in range(i, dimL):
            for k in range(j, dimL):
                a = L.bracket(L.bracket(bL[i],bL[j]),bL[k])
                b = L.bracket(L.bracket(bL[j],bL[k]),bL[i])
                c = L.bracket(L.bracket(bL[k],bL[i]),bL[j])
                if a + b + c != 0:
                    return [bL[i], bL[j], bL[k]]
    return True
#-------------

#-------------
def special_basis_nilp_lie_alg(L):
    r'''
        INPUT:
        
            - a nilpotent lie algebra, L.
            
        OUTPUT:
            
            - a base of the Lie algebra in list format, `{x_{1}, \dots, x_{n}}`, such that for every `x \in L`,
             
            .. MATH::
            
            [x, x_{i}] = \sum_{j = 1}^{i - 1} \lambda_{j}x_{j}.
        
        COMMENT:
        
            The base is calculated by the union of the bases of the ideals below.

            .. MATH::

            0 \subseteq I_{1} \subseteq I_{2} \subseteq \cdots \subseteq I_{\ell} = L,

            where `I_{1}` is the center of `L`, Z(L), and `I_{j}` is the ideal of `L` such that

            .. MATH::

            \faktor{I_{j}}{I_{j-1}} = Z(\faktor{L}{I_{j-1}})
        
        EXAMPLES:
        
            sage: L = nilpotent_lie_algebra(QQ, [6,16], True)
            
            sage: special_basis_nilp_lie_alg(L) 

            [x0, x1, x2, x3, x4, x5]

    '''
    if L.is_nilpotent() == False:
        return print("Essa álgebra não é nilponte.")
    baseL = list(L.basis()) # base de L no formato de lista
    dimL = len(baseL)
    Z = L.center() # centro de L
    baseZ = list(Z.basis())
    dimZ = len(baseZ) 
    baseLnova = [0]*dimL # base que vamos devolver ao final da computação
    baseLnovaTrun = [0]*dimZ 
    for i in range(dimZ):
        baseLnova[i] = baseZ[i]
        baseLnovaTrun[i] = baseZ[i]
    Q = L.quotient(Z) # quociente L/Z
    ZQ = Q.center() # centro de L/Z
    baseZQ = list(ZQ.basis()) # base de Z(L/Z)
    dimZQ = len(baseZQ)
    preImBaseZQ = [0]*dimZQ
    for i in range(dimZ, dimZ + dimZQ):
        baseLnova[i] = Q.lift(baseZQ[i - dimZ]) 
        preImBaseZQ[i - dimZ] = Q.lift(baseZQ[i - dimZ]) # pré-imagem da base de Z(L/Z)
    baseLnovaTrun = baseLnovaTrun + preImBaseZQ
    cont = len(baseLnovaTrun)
    while cont < dimL:
        I = L.ideal(baseLnovaTrun)
        Q = L.quotient(I)
        ZQ = Q.center()
        baseZQ = list(ZQ.basis())
        dimZQ = len(baseZQ)
        preImBaseZQ = [0]*dimZQ
        for i in range(cont, cont+dimZQ):
            baseLnova[i] = Q.lift(baseZQ[i - cont])
            preImBaseZQ[i - cont] = Q.lift(baseZQ[i - cont])
        cont = cont + dimZQ
        baseLnovaTrun = baseLnovaTrun + preImBaseZQ
    return baseLnova
#-------------

#-------------
def coord_base(L, bLdada, x):
    r'''
        INPUT:
        
            - a Lie algebra
            - a basis for Lie algebra in list format
            - an element of Lie algebra

        OUTPUT:
            
            - coordinates of the element, according to the given base.
        
        COMMENT:
        
            If `B = {x_{1}, \dots, x_{n}}` is a basis of `L` and `x = \sum_{i=1}^{n}\lambda_{i}x_{i}`, then the function will return the list `[\lambda_{1}, \dots, \lambda_{n}]`.
        
        EXAMPLES:
        
            sage: L = nilpotent_lie_algebra(QQ, [6,16], True)

            sage: bEspL = special_basis_nilp_lie_alg(L)

            sage: x = L.random_element()

            sage: x

            -2*x1

            sage: coord_base(L, bEspL, x)

            [0, -2, 0, 0, 0, 0]
        
    '''
    bLstr = list(L.basis().keys())
    dimL = len(bLdada)
    coordParbLdada = [0]*dimL
    coordbLdada = [[0]*dimL for i in range(dimL)]
    for i in range(dimL):
        coordParbLdada[i] = [list(j) for j in bLdada[i]]
    for i in range(dimL):
        for j in range(len(coordParbLdada[i])):
            for k in range(dimL):
                if bLstr[k] == coordParbLdada[i][j][0]:
                    coordbLdada[i][k] = coordParbLdada[i][j][1]
                    k = dimL
    coordParx = [list(i) for i in x]
    coordx = [0]*dimL # coordenada de x na base bL
    # preenche a variável coord
    for i in range(len(coordParx)):
        for j in range(len(bLdada)):
            if bLstr[j] == coordParx[i][0]:
                coordx[j] = coordParx[i][1]
                j = len(bLdada)
    V = VectorSpace(QQ, dimL)
    bLdadaVec = [vector(QQ, coordbLdada[i]) for i in range(dimL)]
    W = V.subspace_with_basis(bLdadaVec)
    xvec = vector(QQ, coordx)
    coord = W.coordinates(xvec)
    return coord
#-------------

#-------------
def base_change_nilp_lie_alg(L, bL):
    '''
        INPUT:
        
            - a Lie algebra
            - a basis for Lie algebra in list format
            
        OUTPUT:
            
            - the Lie algebra whose basis is the one given
        
        COMMENT:
        
            This function returns the same Lie algebra given as input, but with the basis changed to the one provided as input.
        
        EXAMPLES:
        
            sage: L = nilpotent_lie_algebra(QQ, [6,16])

            sage: bEspL = special_basis_nilp_lie_alg(L)

            sage: Lesp = base_change_nilp_lie_alg(L, bEspL)

    '''
    dimL = len(bL)
    P = PolynomialRing(QQ, dimL, "x")
    x = list(P.gens())
    xstr = [str(i) for i in x] # ['x0', ..., 'xn']
    var = ', '.join(str(i) for i in x) # 'x0, ..., xn'
    dicio = {} # dicionário para definir nova álgebra de Lie
    # escrever o dicionário dicio {('xi','xj'): {'x0':num0, 'x1':num1, ...}, ... }
    for i in range(dimL):
        for j in range(i, dimL):
            dicioAux = {}
            v = coord_base(L, bL, L.bracket(bL[i],bL[j])) # escrever [xi,xj] na base especial
            for k in range(dimL):
                dicioAux[xstr[k]] = v[k]
            dicio[(xstr[i], xstr[j])] = dicioAux
    LnovaBase = LieAlgebra(QQ, dicio, names = var)
    return LnovaBase
#-------------

#-------------
def nilpotent_lie_algebra( F, args, standard_basis = False ):
    '''
        INPUT:
        
            - a field, F
            - a list with arguments, args
            - an optional Boolean value
            
        OUTPUT:
            
            - a nilpotent Lie algebra 
        
        COMMENT:
        
            This function uses the gap.NilpotentLieAlgebra function to return a nilpotent Lie algebra in Gap. After that, it uses the lie_alg_from_gap_to_sage function to convert that algebra, written in the Gap language, to an algebra in the Sage language. If the optional value of the function is True, then the function changes the basis of the algebra to the basis given by the special_basis_nilp_lie_alg function.
        
        EXAMPLES:
        
            sage: L = nilpotent_lie_algebra(QQ, [6,16])

            sage: L

            Lie algebra on 6 generators (x0, x1, x2, x3, x4, x5) over Rational Field
        
    '''
    L = lie_alg_from_gap_to_sage( gap.NilpotentLieAlgebra( F, args ))
    if not standard_basis:
        return L
    else:
        return base_change_nilp_lie_alg( L, special_basis_nilp_lie_alg( L ))
#-------------

#-------------
def solvable_lie_algebra( F, args, standard_basis = False ):
    '''
        INPUT:
        
            - a field, F
            - a list with arguments, args
            - an optional Boolean value
            
        OUTPUT:
            
            - a solvable Lie algebra 
        
        COMMENT:
        
            This function uses the gap.SolvableLieAlgebra function to return a solvable Lie algebra in Gap. After that, it uses the lie_alg_from_gap_to_sage function to convert that algebra, written in the Gap language, to an algebra in the Sage language. If the optional value of the function is True, then the function changes the basis of the algebra to the basis given by the special_basis_nilp_lie_alg function.
        
        EXAMPLES:
        
            sage: L = solvable_lie_algebra(QQ, [3,2])

            sage: L

            Lie algebra on 3 generators (x0, x1, x2) over Rational Field
        
    '''
    L = lie_alg_from_gap_to_sage( gap.SolvableLieAlgebra( F, args ))
    return L
    
#-------------


#-------------
def isomorphic_random_nilp_lie_alg(L):
    '''
        INPUT:
        
            - a Lie algebra
            
        OUTPUT:
            
            - a Lie algebra
        
        COMMENT:
        
            This function returns a Lie algebra isomorphic to the given Lie algebra, with a basis formed by random elements.
        
        EXAMPLES:
        
            sage: L = nilpotent_lie_algebra(QQ, [6,16])

            sage: isomorphic_random_nilp_lie_alg(L)

            Lie algebra on 6 generators (x0, x1, x2, x3, x4, x5) over Rational Field
        
    '''
    bL = L.basis().list()
    dimL = len(bL)
    V = VectorSpace(QQ, dimL)
    bLnova = [0]*dimL
    bV = [0]*dimL
    logAux = True
    while logAux:
        for i in range(dimL):
            bV[i] = V.random_element()
        if V.are_linearly_dependent(bV) == False:
            logAux = False
    for i in range(dimL):
        for j in range(dimL):
            bLnova[i] = bLnova[i] + bV[i][j]*bL[j]
    LnovaBase = base_change_nilp_lie_alg(L, bLnova)
    return LnovaBase
#-------------

#-------------
def structure_constants(L, bLdada):
    '''
        INPUT:
        
            - a Lie algebra
            - a basis for the Lie algebra
            
        OUTPUT:
            
            - the matrix with the structure constants of the Lie algebra
        
        COMMENT:
        
            Given a Lie algebra L and a basis `{x_{1}, ..., x_{n}}` of `L`, this function returns the matrix of structure constants, `M`, where `[M]_{ij} = [x{i}, x_{j}]`.
        
        EXAMPLES:
        
            sage: L = nilpotent_lie_algebra(QQ, [6,16], True)

            sage: bEspL = special_basis_nilp_lie_alg(L)

            sage: structure_constants(L, bEspL)

            [  0   0   0   0   0   0]
            [  0   0   0   0   0 -x0]
            [  0   0   0  x0 -x1   0]
            [  0   0 -x0   0 -x2   0]
            [  0   0  x1  x2   0  x3]
            [  0  x0   0   0 -x3   0]
        
    '''
    bL = list(L.basis())
    dimL = len(bL)
    R = PolynomialRing(QQ, dimL, L.basis().keys().list())
    FracR = FractionField(R)
    x = R.gens()
    M = matrix(FracR, dimL)
    for j in range(dimL):
        for i in range(dimL):
            for k in range(dimL):
                M[i,j] = M[i,j] + coord_base(L, bL, L.bracket(bL[j],bL[i]))[k]*x[k]
    P = matrix(QQ, dimL)
    for i in range(dimL):
        P[:,i] = vector(coord_base(L, bL, bLdada[i]))
    Pinv = P.inverse()
    E = Pinv*M*P
    return E
#-------------

#-------------
def adjoint_matrix_element(L, x):
    '''
        INPUT:
        
            - a Lie algebra, L
            - an element of L, x
            
        OUTPUT:
            
            - matrix of the adjoint operator of x with respect to the special basis
        
        COMMENT:
        
            Given a Lie algebra `L` and an element `x` from `L`, this function returns the matrix of the `ad_{x}` operator according to the basis given by the special_basis_nilp_lie_alg function.
        
        EXAMPLES:
        
            sage: L = nilpotent_lie_algebra(QQ, [6,16], True) 

            sage: bEspL = special_basis_nilp_lie_alg(L)

            sage: x = L.random_element()

            sage: x

            -x5

            sage: adjoint_matrix_element(L,x)
            [ 0 -1  0  0  0  0]
            [ 0  0  0  0  0  0]
            [ 0  0  0  0  0  0]
            [ 0  0  0  0  1  0]
            [ 0  0  0  0  0  0]
            [ 0  0  0  0  0  0]
        
    '''
    bEspL = special_basis_nilp_lie_alg(L)
    dimL = len(bEspL)
    M = matrix(QQ, dimL)
    colchXbase = [0]*dimL
    for j in range(dimL):
        colchXbase[j] = L.bracket(x,bEspL[j])
    for j in range(dimL):
        for i in range(dimL):
            M[i, j] = coord_base(L, bEspL, colchXbase[j])[i]
    return M    
#-------------

#-------------
def get_keys_from_value(d, val):
    return [k for k, v in d.items() if v == val]
#-------------

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


def standard_filiform_lie_algebra( n ):

    vars = [ 'y'+str(k) for k in range( n-1 ) ] + ['x']
    rel_dict = { ('x','y'+str(k)): {'y'+str(k-1): 1} for k in range(1,n-1) }
    return LieAlgebra( QQ, rel_dict, vars )
#-------------

#-------------
def solvable_lie_triangular_base(L):
    dimL = L.dimension()
    triagBase = [0]*dimL
    D = L.derived_series()
    I = D[len(D)-2]
    dimI = I.dimension()
    triagBaseTrunc = [0]*dimI
    for i in range(dimI):
        triagBase[i] = I.gens()[i]
        triagBaseTrunc[i] = I.gens()[i]
    while triagBase[dimL-1] == 0:
        Q = L.quotient(I)
        D = Q.derived_series()
        Iaux = D[len(D)-2]
        preImaIaux = [0]*(Iaux.dimension())
        for i in range(Iaux.dimension()):
            preImaIaux[i] = Q.lift(Iaux.gens()[i])
            triagBase[i+dimI] = preImaIaux[i]
        triagBaseTrunc = triagBaseTrunc + preImaIaux
        I = L.ideal(triagBaseTrunc)
        dimI = I.dimension()
    return triagBase
#-------------

#-------------
def matrix_adjoint_element(L, bL, x):
    dimL = L.dimension()
    M = Matrix(QQ,dimL)
    for j in range(dimL):
        for i in range(dimL):
             M[i,j] = M[i,j] + coord_base(L, bL, L.bracket(x,bL[j]))[i]
    return M
#-------------