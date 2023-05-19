#------------- FUNÇÃO algLieGAPparaSAGE -------------
def algLieGAPparaSAGE(l):
    '''Dada uma álgebra de Lie em GAP, entrega uma álgebra de Lie em SAGE'''
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

#------------- FUNÇÃO satisIdentJacobi 
def satisIdentJacobi(L):
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

#------------- FUNÇÃO baseEspAlgNil -------------
def baseEspAlgNil(L):
    '''Dada uma álgebra de Lie nilpotente L, a função devolve uma base especial de L'''
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

#------------- FUNÇÃO coordBase -------------
def coordBase(x, L, bLdada):
    '''Devolve as coordenadas de x na base bLdada'''
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

#------------- FUNÇÃO algLieNilMudBase -------------
def algLieNilMudBase(L, bL):
    '''Dado uma álgebra de Lie nilpotente e uma para para ela, a função devolve uma álgebra de Lie nilpotente isomorfa, com a base dada'''
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
            v = coordBase(L.bracket(bL[i],bL[j]), L, bL) # escrever [xi,xj] na base especial
            for k in range(dimL):
                dicioAux[xstr[k]] = v[k]
            dicio[(xstr[i], xstr[j])] = dicioAux
    LnovaBase = LieAlgebra(QQ, dicio, names = var)
    return LnovaBase
#-------------

#------------- FUNÇÃO algLieIsoAlea -------------
def algLieIsoAlea(L):
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
    LnovaBase = algLieNilMudBase(L, bLnova)
    return LnovaBase
#-------------

#------------- FUNÇÃO consEst -------------
def consEst(L, bLdada):
    bL = list(L.basis())
    dimL = len(bL)
    R = PolynomialRing(QQ, dimL, "x")
    FracR = FractionField(R)
    x = R.gens()
    M = matrix(FracR, dimL)
    for i in range(dimL):
        for j in range(dimL):
            for k in range(dimL):
                M[i,j] = M[i,j] + coordBase(L.bracket(bL[i],bL[j]), L, bL)[k]*x[k]
    P = matrix(QQ, dimL)
    for i in range(dimL):
        P[:,i] = vector(coordBase(bLdada[i], L, bL))
    Pinv = P.inverse()
    E = Pinv*M*P
    return E
#-------------

#------------- FUNÇÃO matrizDer -------------
def matrizDer(x, L):
    bLesp = baseEspAlgNil(L)
    dimL = len(bLesp)
    M = matrix(QQ, dimL)
    colchXbase = [0]*dimL
    for i in range(dimL):
        colchXbase[i] = L.bracket(x,bLesp[i])
    for i in range(dimL):
        for j in range(dimL):
            M[i, j] = coordBase(colchXbase[i], L, bLesp)[j]
    return M    
#-------------