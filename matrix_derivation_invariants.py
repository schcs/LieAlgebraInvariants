from sage.all import ZZ, Matrix, QQ, lcm, PolynomialRing, vector
# from sage.rings.number_field.splitting_field import SplittingFieldAbort

#-------------
def derivation_of_matrix(M):
    K = M.base_ring()
    n = M.nrows()
    F = PolynomialRing( K, "x", n ).fraction_field()
    gensF = F.gens()
    D = F.derivation_module()
    bD = D.basis().list()

    return sum( M[i,j]*gensF[i]*bD[j] for i in range( n ) for j in range( n ))

#-------------

#-------------
def matrix_of_derivation(diff):

    P = diff.parent().base()
    gensP = P.gens()
    n = len(gensP)
    M = Matrix(P.base_ring(), n)
    mc = diff.monomial_coefficients()

    for i in range(n):
        pol = mc[i].numerator() if i in mc.keys() else P.zero().numerator()
        mon = pol.monomials()
        coeff = pol.coefficients()
        mon_coeff_dict = dict( zip( mon, coeff ))
        for j in range(len(mon)):
            cont = gensP.index( mon[j] )
            M[cont,i] = mon_coeff_dict[mon[j]]

    return M
#-------------

#-------------

# the function implements the computation of the invariants given
# in Lemma 3.1 of Snobl and Winternitz for the Jordan block with
# zero eigenvaliue

def invariants_nilpotent_jordan_block_lemma_3_1(gensF):
    n = len(gensF)
    if n == 1:
        return []
    a = [0]*(n-1)
    a[0] = gensF[0]
    for k in range(2,n):
        fact = 1
        for j in range(k+1):
            fact = fact*j if j > 0 else 1
            a[k-1] = a[k-1] + ((-1)**j)*(gensF[0]**(k-1-j))*(gensF[1]**j)*(gensF[k-j])/fact
        #a[k-1] = a[k-1].numerator()
    return a
#-------------

#-------------
def invariants_eigenvalue_jordan_block(gensF):

    v = invariants_nilpotent_jordan_block_lemma_3_1(gensF)
    return [ v[i+1]/v[0]**(i+2) for i in range(len( v )-1)]
#-------------

#-------------
def invariants_diagonal(diag, gensF, K):
    len_diag = len(diag)
    degree_K = K.polynomial().degree()
    for i in range(len_diag):
        diag[i] = K(diag[i])
    coeff = [0]*len_diag
    for i in range(len_diag):
        coeff_diag = list(diag[i].polynomial())
        if len(coeff_diag) < degree_K:
            for j in range(degree_K-len(coeff_diag)):
                coeff_diag = coeff_diag + [0]
        coeff[i] = coeff_diag
    sist_qq = Matrix(QQ, degree_K, len_diag)
    for i in range(len_diag):
        sist_qq[:,i] = vector(coeff[i])
    for i in range(degree_K):
        den_line = []
        for j in range(len_diag):
            den_line = den_line + [sist_qq[i,j].denominator()]
        lcm_diag = lcm(den_line)
        sist_qq[i,:] = sist_qq[i,:]*lcm_diag
    sist_zz = Matrix(ZZ, degree_K, len_diag)
    for i in range(degree_K):
        sist_zz[i,:] = sist_qq[i,:]
    sol = sist_zz.right_kernel()
    sol_basis = sol.basis()
    if len(sol_basis) == 0:
        return []
    matrix_sol = Matrix(ZZ,len(sol_basis),len_diag)
    for i in range(len(sol_basis)):
        matrix_sol[i,:] = sol_basis[i]
    matrix_sol_hermite = matrix_sol.hermite_form()
    inv = []
    for i in range(len(sol_basis)):
        inv_comp = 1
        for j in range(len_diag):
            inv_comp = inv_comp*gensF[j]**matrix_sol_hermite[i,j]
        inv = inv + [inv_comp]
    return inv
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
        sage: diff = derivation_of_matrix(M)
        sage: inv = invariants_matrix_derivation(diff)

    '''
    # Extending the matrix to a field where its Jordan form can be taken
    M = matrix_of_derivation(diff)
    f = M.characteristic_polynomial()
    K = f.splitting_field("a")
    M = M.change_ring(K)
    J, P = M.jordan_form(transformation=True)
    diff = derivation_of_matrix(M)
    D = diff.parent()
    F = D.base()
    gensF = F.gens()
    # Change of basis
    gensAlt = [0]*len(gensF)
    for j in range(len(gensF)):
        for i in range(len(gensF)):
            gensAlt[j] = gensAlt[j] + P[i,j]*gensF[i]
    # Non-repeated eigenvalues
    eigenVR = J.eigenvalues()
    eigenV = [J.eigenvalues()[0]]
    cont = 0
    for i in range(1, len(eigenVR)):
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
    # List with the data used to calculate the invariants. Each position of this list is another list with 2 components: the first one, with the eigenvalue; the second one, with the variables corresponding to its block in Jordan form
    dateInv = [0]*(numBlocks)
    for i in range(len(dateInv)):
        gensBlock = []
        for j in range(blocks[i].nrows()):
            gensBlock = gensBlock + [gensAlt[acumLenBlocks[i]+j]]
        dateInv[i] = [[blocks[i][0,0]], gensBlock]
    # Interaction between blocks
    diag = [] # first entries of each block (except the nilpotent one of size 1)
    gensDiag = [] # variables of diag
    blockNil = [] # position in dateInv at which there is a nilpotent block with size greater than 1
    blockNotNil = [] # position in dateInv at which there is a non-nilpotent block with size greater than 1
    for i in range(len(dateInv)):
        if dateInv[i][0][0] != 0:
            diag = diag + [dateInv[i][0][0]]
            gensDiag = gensDiag + [dateInv[i][1][0]]
            if len(dateInv[i][1]) != 1:
                blockNotNil = blockNotNil + [i]
        else:
            if len(dateInv[i][1]) == 1:
                diag = diag + [dateInv[i][0][0]]
                gensDiag = gensDiag + [dateInv[i][1][0]]
            else:
                blockNil = blockNil + [i]
    # Calculation of the invariants
    inv = []
    for i in range(len(dateInv)):
        if len(dateInv[i][1]) != 1:
            if dateInv[i][0][0] == 0:
                inv = inv +  invariants_nilpotent_jordan_block_lemma_3_1(dateInv[i][1])
            else:
                inv = inv + invariants_eigenvalue_jordan_block(dateInv[i][1])
    if len(diag) >= 1:
        inv = inv + invariants_diagonal(diag, gensDiag, K)
    if len(blockNotNil) == 0:
        if len(blockNil) > 1:
            for i in range(1,len(blockNil)):
                inv = inv + [dateInv[blockNil[0]][1][1]/dateInv[blockNil[0]][1][0] - dateInv[blockNil[i]][1][1]/dateInv[blockNil[i]][1][0]]
    else:
        for j in range(len(blockNil)):
            inv = inv + [dateInv[blockNotNil[0]][1][1]/dateInv[blockNotNil[0]][1][0] - dateInv[blockNil[j]][1][1]/dateInv[blockNil[j]][1][0]]
        for i in range(1,len(blockNotNil)):
            inv = inv + [dateInv[blockNotNil[i]][1][1]/dateInv[blockNotNil[i]][1][0] - dateInv[blockNil[0]][1][1]/dateInv[blockNil[0]][1][0]]
    return inv
#-------------
