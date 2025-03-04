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

def lie_alg_family_6_19(F):

    K = FunctionField(F, 'ε')
    ε = K.gens()[0]
    mult_table = {('x2','x5') : {'x1' : -1},
                  ('x3','x4') : {'x1' : -1},
                  ('x3','x6') : {'x1' : -ε},
                  ('x4','x5') : {'x2' : 1},
                  ('x4','x6') : {'x3' : 1}}

    return LieAlgebra(K, mult_table, names = 'x1,x2,x3,x4,x5,x6')

def lie_alg_family_6_21(F):

    K = FunctionField(F, 'ε')
    ε = K.gens()[0]
    mult_table = {('x2','x5') : {'x1' : -1},
                  ('x3','x6') : {'x1' : -ε},
                  ('x4','x5') : {'x2' : -1},
                  ('x4','x6') : {'x3' : -1},
                  ('x5','x6') : {'x4' : 1}}

    return LieAlgebra(K, mult_table, names = 'x1,x2,x3,x4,x5,x6')

def lie_alg_family_6_22(F):

    K = FunctionField(F, 'ε')
    ε = K.gens()[0]
    mult_table = {('x3','x4') : {'x1' : 1},
                  ('x3','x5') : {'x2' : 1},
                  ('x4','x6') : {'x2' : ε},
                  ('x5','x6') : {'x1' : 1}}

    return LieAlgebra(K, mult_table, names = 'x1,x2,x3,x4,x5,x6')

def lie_alg_family_6_24(F):

    K = FunctionField(F, 'ε')
    ε = K.gens()[0]
    mult_table = {('x3','x5') : {'x1' : -1},
                  ('x3','x6') : {'x2' : -1},
                  ('x4','x5') : {'x2' : -ε},
                  ('x4','x6') : {'x1' : -1},
                  ('x5','x6') : {'x3' : 1}}

    return LieAlgebra(K, mult_table, names = 'x1,x2,x3,x4,x5,x6')


def nilpotent_6_19(F, trian = False):
    P = PolynomialRing(F, "a")
    a = list(P.gens())[0]
    K = P.fraction_field()
    d = {('x0','x1'): {'x3':1}, ('x0','x2'): {'x4':1}, ('x0','x4'): {'x5':1}, ('x1','x3'): {'x5':1}, ('x2','x4'): {'x5':K(a)}}
    L = LieAlgebra(K, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L

def nilpotent_6_21(F, trian = False):
    P = PolynomialRing(F, "a")
    a = list(P.gens())[0]
    K = P.fraction_field()
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x3':1}, ('x1','x2'): {'x4':1}, ('x0','x3'): {'x5':1}, ('x1','x4'): {'x5':a}}
    L = LieAlgebra(K, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_22(F, trian = False):
    P = PolynomialRing(F, "a")
    a = list(P.gens())[0]
    K = P.fraction_field()
    d = {('x0','x1'): {'x4':1}, ('x0','x2'): {'x5':1}, ('x1','x3'): {'x5':a}, ('x2','x3'): {'x4':1}}
    L = LieAlgebra(K, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

def nilpotent_6_24(F, trian = False):
    P = PolynomialRing(F, "a")
    a = list(P.gens())[0]
    K = P.fraction_field()
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x4':1}, ('x0','x3'): {'x5':a}, ('x1','x2'): {'x5':1}, ('x1','x3'): {'x4':1}}
    L = LieAlgebra(K, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
