from sage.all import QQ, LieAlgebra, libgap, PolynomialRing, gap, ZZ, matrix, vector, FractionField, Matrix, VectorSpace

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
    return LieAlgebra(QQ, dicio, names = var)


# libgap version

def lie_alg_from_libgap_to_sage(l):
    '''
        INPUT:

            - a Lie algebra written in Gap via libgap

        OUTPUT:

            - a Lie algebra written in Sage.

        COMMENT:

            This function converts a Lie algebra, written in the Gap language, into a Lie algebra written in the Sage language.

        EXAMPLES:

            sage: l = libgap.NilpotentLieAlgebra(QQ, [6,14])

            sage: L = lie_alg_from_gap_to_sage(l)

            sage: L

            Lie algebra on 6 generators (x0, x1, x2, x3, x4, x5) over Rational Field

    '''
    dimL = libgap.Dimension(l).sage()
    bL = l.Basis()
    var = [f"x{i}" for i in range(dimL)]  # ['x0', ..., 'xn']
    dicio = {}  # dicionário para definir nova álgebra de Lie
    # escrever o dicionário dicio {('xi','xj'): {'x0':num0, 'x1':num1, ...}, ... }
    for i in range(dimL):
        for j in range(i, dimL):
            dicioAux = {}
            v = bL.Coefficients(bL[i] * bL[j])  # escrever [xi,xj] na base especial
            for k in range(dimL):
                dicioAux[var[k]] = v[k]
            dicio[(var[i], var[j])] = dicioAux
    return LieAlgebra(QQ, dicio, names=var)


#-------------

#-------------
def jacobi_satisfies(L):
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

            sage: jacobi_satisfies(L)

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

            sage: bEspL = triangular_basis_lie_algebra(L)

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
    F = L.base()
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
    V = VectorSpace(F, dimL)
    bLdadaVec = [vector(F, coordbLdada[i]) for i in range(dimL)]
    W = V.subspace_with_basis(bLdadaVec)
    xvec = vector(F, coordx)
    coord = W.coordinates(xvec)
    return coord
#-------------

#-------------
def base_change_lie_algebra(L, bL):
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

            sage: bEspL = triangular_basis_lie_algebra(L)

            sage: Lesp = base_change_lie_algebra(L, bEspL)

    '''
    dimL = len(bL)
    F = L.base()
    P = PolynomialRing(F, dimL, "x")
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
    LnovaBase = LieAlgebra(F, dicio, names = var)
    return LnovaBase
#-------------

#-------------
def nilpotent_lie_algebra2( F, args, standard_basis = False ):
    '''
        INPUT:

            - a field, F
            - a list with arguments, args
            - an optional Boolean value

        OUTPUT:

            - a nilpotent Lie algebra

        COMMENT:

            This function uses the gap.NilpotentLieAlgebra function to return a nilpotent Lie algebra in Gap. After that, it uses the lie_alg_from_gap_to_sage function to convert that algebra, written in the Gap language, to an algebra in the Sage language. If the optional value of the function is True, then the function changes the basis of the algebra to the basis given by the triangular_basis_lie_algebra function.

        EXAMPLES:

            sage: L = nilpotent_lie_algebra(QQ, [6,16])

            sage: L

            Lie algebra on 6 generators (x0, x1, x2, x3, x4, x5) over Rational Field

    '''
    L = lie_alg_from_libgap_to_sage( libgap.NilpotentLieAlgebra( F, args ))
    if not standard_basis:
        return L
    else:
        return base_change_lie_algebra( L, triangular_basis_lie_algebra( L ))
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

            This function uses the gap.SolvableLieAlgebra function to return a solvable Lie algebra in Gap. After that, it uses the lie_alg_from_gap_to_sage function to convert that algebra, written in the Gap language, to an algebra in the Sage language. If the optional value of the function is True, then the function changes the basis of the algebra to the basis given by the triangular_basis_lie_algebra function.

        EXAMPLES:

            sage: L = solvable_lie_algebra(QQ, [3,2])

            sage: L

            Lie algebra on 3 generators (x0, x1, x2) over Rational Field

    '''
    return lie_alg_from_libgap_to_sage( libgap.SolvableLieAlgebra( F, args ))

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
            a = [0]*dimL
            for j in range(dimL):
                a[j] = ZZ.random_element(-3,3)
            bV[i] = V(a)
        if not V.are_linearly_dependent(bV):
            logAux = False
    for i in range(dimL):
        for j in range(dimL):
            bLnova[i] = bLnova[i] + bV[i][j]*bL[j]
    LnovaBase = base_change_lie_algebra(L, bLnova)
    return LnovaBase
#-------------

#-------------
def structure_constants(L, bLdada = False):
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

            sage: bEspL = triangular_basis_lie_algebra(L)

            sage: structure_constants(L, bEspL)

            [  0   0   0   0   0   0]
            [  0   0   0   0   0 -x0]
            [  0   0   0  x0 -x1   0]
            [  0   0 -x0   0 -x2   0]
            [  0   0  x1  x2   0  x3]
            [  0  x0   0   0 -x3   0]

    '''
    if bLdada == False:
        bL = list(L.basis())
    else:
        bL = bLdada
    dimL = len(bL)
    F = L.base()
    R = PolynomialRing(F, dimL, L.basis().keys().list())
    FracR = FractionField(R)
    x = R.gens()
    M = matrix(FracR, dimL)
    for i in range(dimL):
        for j in range(dimL):
            for k in range(dimL):
                M[i,j] = M[i,j] + L.bracket(bL[i],bL[j]).to_vector()[k]*x[k]
    return M
#-------------

#-------------
def get_keys_from_value(d, val):
    return [k for k, v in d.items() if v == val]
#-------------

#-------------
def triangular_basis_nilpotent_lie_algebra(L):
    r'''
        INPUT:

            -

        OUTPUT:

            -

        COMMENT:



        EXAMPLES:



    '''
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
def triangular_basis_solvable_lie_algebra(L):
    r'''
        INPUT:

            -

        OUTPUT:

            -

        COMMENT:



        EXAMPLES:



    '''
    dimL = L.dimension()
    bL = L.basis().list()
    V = VectorSpace(QQ,dimL)
    serie = L.derived_series()
    n = len(serie)
    serieV = [0]*n
    for i in range(n):
        bS = serie[i].gens()
        serieV[i] = V.subspace_with_basis([V(coord_base(L,bL,bS[j])) for j in range(len(bS))])
    Q = serieV[n-2].quotient(serieV[n-1])
    bVserie = [V(Q.lift(Q.basis()[i])) for i in range(len(Q.basis()))]
    aux = n-2
    while aux !=0:
        aux = aux - 1
        Q = serieV[aux].quotient(serieV[aux+1])
        bVserie = bVserie + [V(Q.lift(Q.basis()[i])) for i in range(len(Q.basis()))]
    baseTriang = []
    for i in range(len(bVserie)):
        a = 0
        for j in range(dimL):
            a = a + bVserie[i][j]*bL[j]
        baseTriang = baseTriang + [a]
    return baseTriang
#-------------

#-------------
def triangular_basis_lie_algebra(L):
    r'''
        INPUT:

            -

        OUTPUT:

            -

        COMMENT:



        EXAMPLES:



    '''
    if L.is_nilpotent():
        return triangular_basis_nilpotent_lie_algebra(L)
    if L.is_solvable():
        return triangular_basis_solvable_lie_algebra(L)
#-------------

#-------------
def adjoint_matrix_element(L, bL, x):
    dimL = L.dimension()
    M = Matrix(QQ,dimL)
    for j in range(dimL):
        for i in range(dimL):
            M[i,j] = M[i,j] + coord_base(L, bL, L.bracket(x,bL[j]))[i]
    return M
#-------------
