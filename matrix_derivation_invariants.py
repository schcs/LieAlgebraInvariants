from sage.rings.number_field.splitting_field import SplittingFieldAbort

#-------------
def matrix_for_derivation(M):
    K = M.base_ring()
    n = M.nrows()
    P = PolynomialRing(K, "x", n)
    F = P.fraction_field()
    gensF = F.gens()
    D = F.derivation_module()
    bD = D.basis().list()
    diff = D.zero()
    for j in range(n):
        for i in range(n):
            diff = diff + M[i,j]*gensF[i]*bD[j]
    return diff
#-------------

#-------------
def derivation_for_matrix(diff):
    diff = diff.extend_to_fraction_field()
    D = diff.parent()
    F = D.base()
    gensF = F.gens()
    n = len(gensF)
    P = F.base()
    gensP = P.gens()
    M = Matrix(F.base_ring(), n)
    for i in range(n):
        pol = P(diff(gensF[i]))
        mon = pol.monomials()
        coeff = pol.coefficients()
        for j in range(len(mon)):
            varBool = True
            cont = 0
            while varBool:
                if mon[j] == gensP[cont]:
                    M[cont,i] = coeff[j]
                    varBool = False
                else:
                    cont = cont + 1
                if cont == n:
                    varBool = False
    return M
#-------------

#-------------
def invariants_zero_jordan_block(gensF):
    n = len(gensF)
    if n == 1:
        return []
    a = [0]*(n-1)
    a[0] = gensF[0]
    for i in range(2,n):
        fact = 1
        for j in range(i+1):
            fact = fact*j
            if fact == 0:
                fact = 1
            invFact = QQ(1/fact)
            num = ((-1)**j)*invFact
            a[i-1] = a[i-1] + num*(gensF[0]**(i-1-j))*(gensF[1]**j)*(gensF[i-j])
        a[i-1] = a[i-1].numerator()
    return a
#-------------

#-------------
def invariants_eigenvalue_jordan_block(gensF):
    v = invariants_zero_jordan_block(gensF)
    if len(v) == 0:
        return []
    a = [0]*(len(v)-1)
    for i in range(len(a)):
        a[i] = v[i+1]/v[0]**(i+2)
    return a
#-------------

#-------------
def invariants_diagonal(diag, gensF):
    n = len(diag)
    a = []
    varBool = True
    cont = 0
    while varBool:
        for i in range(cont+1,n):
            if diag[cont]/diag[i] in QQ:
                quoc = diag[cont]/diag[i]
                den = quoc.denominator()
                num = quoc*den
                a = a + [(gensF[cont]**(den))*(gensF[i]**(num)).inverse()]
        cont = cont + 1
        if len(a) == n-1 or cont == n-1:
            varBool = False
    return a
#-------------

#-------------
def invariants_matrix_derivation(diff):
    r'''
        INPUT:
        
            - Derivation over a field of fractions defined from a matrix.
            
        OUTPUT:
            
            - List of algebraically independent generators of the rational invariant algebra of the derivation.
            
        COMMENT:

        Given a square matrix M of size n x n, we can define a derivation D_{M} in the polynomial algebra in the variables x_{1}, \cdots, x_{n}, by setting D_{M} = \sum_{j=1}^{n}\sum_{i=1}^{n}D_{M}[i,j]D_{x_{i}}. From the Jordan matrix of D_{M}​, we determine the invariants of the rational function algebra of D_{M}​.

        EXAMPLES:

        sage: j1 = jordan_block(2,1)
        sage: j2 = jordan_block(3,1)
        sage: j3 = jordan_block(5,1)
        sage: j4 = jordan_block(4,2)
        sage: j5 = jordan_block(0,1)
        sage: j6 = jordan_block(0,3)
        sage: j7 = jordan_block(1,3)
        sage: j8 = jordan_block(1,4)
        sage: j9 = jordan_block(0,3)
        sage: M = block_diagonal_matrix(j1,j4,j5,j6,j7,j8,j2,j3,j9)
        sage: diff = matrix_for_derivation(M)
        sage: inv = invariants_matrix_derivation(diff)

    '''
    # Extending the matrix to a field where its Jordan form can be taken
    M = derivation_for_matrix(diff)
    n = M.nrows()
    f = M.characteristic_polynomial()
    K = f.splitting_field("a")
    M = M.base_extend(K)
    J, P = M.jordan_form(transformation=True)
    diff = matrix_for_derivation(M)
    D = diff.parent()
    F = D.base()
    gensF = list(F.gens())
    # Change of basis
    gensAlt = [0]*len(gensF)
    for j in range(len(gensF)):
        for i in range(len(gensF)):
            gensAlt[j] = gensAlt[j] + P[i,j]*gensF[i]
    # Non-repeated eigenvalues
    eigenVR = J.eigenvalues()
    eigenV = [J.eigenvalues()[0]]
    cont = 0
    for i in range(1,len(eigenVR)):
        if eigenVR[i] != eigenV[cont]:
            eigenV = eigenV + [eigenVR[i]]
            cont = cont + 1
    # Number of Jordan blocks
    numBlocks = 0
    for i in range(len(eigenV)):
        numBlocks = numBlocks + J.left_eigenspaces()[i][1].dimension()
    # Jordan blocks
    blocks = []
    lenBlocks = []
    cont = 0
    for i in range(numBlocks):
        blocks = blocks + [J.subdivision(i,i)]
        lenBlocks = lenBlocks + [J.subdivision(i,i).nrows()]
        cont = cont + 1
    # Variables corresponding to each block
    acumLenBlocks = [0]*len(lenBlocks)
    acumLenBlocks[0] = 0
    for i in range(1,len(lenBlocks)):
        acumLenBlocks[i] = acumLenBlocks[i-1] + lenBlocks[i-1]
    # Eigenvalues corresponding to 1 x 1 blocks
    cont = 0
    diag = []
    notDiag = []
    gensDiag = []
    for i in range(numBlocks):
        if blocks[i].nrows() == 1 and blocks[i][0] != 0:
            diag = diag + [blocks[i].list()[0]]
            gensDiag = gensDiag + [gensAlt[cont]]
        else:
            notDiag = notDiag + [i]
        cont = cont + blocks[i].nrows()
    diagExist = 0
    if len(diag) != 0:
        diagExist = 1
    # List with the data used to calculate the invariants. Each position of this list is another list with 2 components: the first one, with the eigenvalue; the second one, with the variables corresponding to its block in Jordan form
    dateInv = [0]*(len(notDiag)+diagExist)
    if diagExist != 0:
        dateInv[0] = [diag, gensDiag]
        for i in range(len(notDiag)):
            gensBlock = []
            for j in range(blocks[notDiag[i]].nrows()):
                        gensBlock = gensBlock + [gensAlt[acumLenBlocks[notDiag[i]]+j]]
            dateInv[diagExist+i] = [[blocks[notDiag[i]][0,0]], gensBlock]
    else:
        for i in range(len(notDiag)):
            gensBlock = []
            for j in range(blocks[notDiag[i]].nrows()):
                        gensBlock = gensBlock + [gensAlt[acumLenBlocks[i]+j]]
            dateInv[i] = [[blocks[notDiag[i]][0,0]], gensBlock]
    # Calculation of the invariants
    inv = []
    for i in range(len(dateInv)):
        if len(dateInv[i][0]) != 1:
            inv = inv + invariants_diagonal(dateInv[i][0], dateInv[i][1])
        else:
            if dateInv[i][0][0] != 0 and dateInv[i][0][0] != 1 and dateInv[i][0][0] != -1:
                inv = inv + invariants_eigenvalue_jordan_block(dateInv[i][1])
            if dateInv[i][0][0] == 1:
                inv = inv + invariants_eigenvalue_jordan_block(dateInv[i][1])
                for j in range(len(dateInv)):
                    for k in range(len(dateInv[j][0])):
                        if j != i and dateInv[j][0][k] in QQ and dateInv[j][0][k] != 0:
                            den = dateInv[j][0][k].denominator()
                            num = dateInv[j][0][k]*den
                            inv = inv + [(dateInv[j][1][k]**(den))*(dateInv[i][1][0]**(num)).inverse()]
            if dateInv[i][0][0] == -1:
                inv = inv + invariants_eigenvalue_jordan_block(dateInv[i][1])
                for j in range(len(dateInv)):
                    if len(dateInv[j][0]) != 1:
                        for k in range(len(dateInv[j][0])):
                            if dateInv[j][0][k] in QQ:
                                den = dateInv[j][0][k].denominator()
                                num = dateInv[j][0][k]*den
                                inv = inv + [(dateInv[j][1][k]**(den))*(dateInv[i][1][0]**(num))]
                    else:
                        for k in range(len(dateInv[j][0])):
                            if j != i and dateInv[j][0][k] in QQ and dateInv[j][0][k] != 0 and dateInv[j][0][k] != 1 and dateInv[j][0][k] != -1:
                                den = dateInv[j][0][k].denominator()
                                num = dateInv[j][0][k]*den
                                inv = inv + [(dateInv[j][1][k]**(den))*(dateInv[i][1][0]**(num))]
            if dateInv[i][0][0] == 0:
                if len(dateInv[i][1]) == 1:
                    inv = inv + dateInv[i][1]
                else:
                    inv = inv + invariants_zero_jordan_block(dateInv[i][1])
    return inv
#-------------
