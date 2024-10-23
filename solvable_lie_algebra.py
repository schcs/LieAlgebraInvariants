#-------------
def solvable_1_1(F, trian = False):
    L = lie_algebras.abelian(F, names='x0')
    return L
#-------------

#-------------
def solvable_2_1(F, trian = False):
    L = lie_algebras.abelian(F, names='x0, x1')
    return L
#-------------

#-------------
def solvable_2_2(F, trian = False):
    d = {('x0','x1'): {'x0':-1}}
    L = LieAlgebra(F, d, names='x0, x1')
    if trian == True:
        bTL = triangular_basis_solvable_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def solvable_3_1(F, trian = False):
    L = lie_algebras.abelian(F, names='x0, x1, x2')
    return L
#-------------

#-------------
def solvable_3_2(F, trian = False):
    d = {('x0','x2'): {'x0':-1}, ('x1','x2'): {'x1':-1}}
    L = LieAlgebra(F, d, names='x0, x1, x2')
    if trian == True:
        bTL = triangular_basis_solvable_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def solvable_3_3(F, trian = False):
    P = PolynomialRing(F, "t")
    t = list(P.gens())[0]
    f = t**2 + 1
    K = NumberField(f, 'a')
    a = list(K.gens())[0]
    d = {('x0','x2'): {'x1':-1}, ('x1','x2'): {'x0':-a, 'x1':-1}}
    L = LieAlgebra(K, d, names='x0, x1, x2')
    if trian == True:
        bTL = triangular_basis_solvable_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def solvable_3_4(F, trian = False):
    P = PolynomialRing(F, "t")
    t = list(P.gens())[0]
    f = t**2 + 1
    K = NumberField(f, 'a')
    a = list(K.gens())[0]
    d = {('x0','x2'): {'x1':-1}, ('x1','x2'): {'x0':-a}}
    L = LieAlgebra(K, d, names='x0, x1, x2')
    if trian == True:
        bTL = triangular_basis_solvable_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def solvable_4_1(F, trian = False):
    L = lie_algebras.abelian(F, names='x0, x1, x2, x3')
    return L
#-------------

#-------------
def solvable_4_2(F, trian = False):
    d = {('x0','x3'): {'x0':-1}, ('x1','x3'): {'x1':-1}, ('x2','x3'): {'x2':-1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3')
    if trian == True:
        bTL = triangular_basis_solvable_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def solvable_4_3(F, trian = False):
    P = PolynomialRing(F, "t")
    t = list(P.gens())[0]
    f = t**2 + 1
    K = NumberField(f, 'a')
    a = list(K.gens())[0]
    d = {('x0','x3'): {'x0':-1}, ('x1','x3'): {'x2':-1}, ('x2','x3'): {'x1':a, 'x2':-a-1}}
    L = LieAlgebra(K, d, names='x0, x1, x2, x3')
    if trian == True:
        bTL = triangular_basis_solvable_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def solvable_4_4(F, trian = False):
    d = {('x1','x3'): {'x2':-1}, ('x2','x3'): {'x2':-1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3')
    if trian == True:
        bTL = triangular_basis_solvable_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def solvable_4_5(F, trian = False):
    d = {('x1','x3'): {'x2':-1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3')
    if trian == True:
        bTL = triangular_basis_solvable_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def solvable_4_6(F, trian = False):
    P1 = PolynomialRing(F, "t1")
    t1 = list(P1.gens())[0]
    f1 = t1**2 - 2
    K1 = NumberField(f1, 'a1')
    a1 = list(K1.gens())[0]
    P2 = PolynomialRing(K1, "t2")
    t2 = list(P2.gens())[0]
    f2 = t2**2 + 1
    K2 = NumberField(f2, 'a2')
    a2 = list(K2.gens())[0]
    d = {('x0','x3'): {'x1':-1}, ('x1','x3'): {'x2':-1}, ('x2','x3'): {'x0':-a1, 'x1':-a2, 'x2':-1}}
    L = LieAlgebra(K2, d, names='x0, x1, x2, x3')
    if trian == True:
        bTL = triangular_basis_solvable_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def solvable_4_7(F, trian = False):
    P1 = PolynomialRing(F, "t1")
    t1 = list(P1.gens())[0]
    f1 = t1**2 - 2
    K1 = NumberField(f1, 'a1')
    a1 = list(K1.gens())[0]
    P2 = PolynomialRing(K1, "t2")
    t2 = list(P2.gens())[0]
    f2 = t2**2 + 1
    K2 = NumberField(f2, 'a2')
    a2 = list(K2.gens())[0]
    d = {('x0','x3'): {'x1':-1}, ('x1','x3'): {'x2':-1}, ('x2','x3'): {'x0':-a1, 'x1':-a2}}
    L = LieAlgebra(K2, d, names='x0, x1, x2, x3')
    if trian == True:
        bTL = triangular_basis_solvable_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def solvable_4_8(F, trian = False):
    d = {('x0','x1'): {'x1':1}, ('x2','x3'): {'x3':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3')
    if trian == True:
        bTL = triangular_basis_solvable_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def solvable_4_9(F, trian = False):
    P = PolynomialRing(F, "t")
    t = list(P.gens())[0]
    f = t**2 + 1
    K = NumberField(f, 'a')
    a = list(K.gens())[0]
    d = {('x0','x2'): {'x0':-1}, ('x0','x3'): {'x0':-1, 'x1':-a}, ('x1','x2'): {'x1':-1}, ('x1','x3'): {'x0':-1}}
    L = LieAlgebra(K, d, names='x0, x1, x2, x3')
    if trian == True:
        bTL = triangular_basis_solvable_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def solvable_4_10(F, trian = False):
    P = PolynomialRing(F, "t")
    t = list(P.gens())[0]
    f = t**2 + 1
    K = NumberField(f, 'a')
    a = list(K.gens())[0]
    d = {('x0','x2'): {'x0':-1}, ('x0','x3'): {'x1':-1}, ('x1','x2'): {'x1':-1}, ('x1','x3'): {'x0':-a}}
    L = LieAlgebra(K, d, names='x0, x1, x2, x3')
    if trian == True:
        bTL = triangular_basis_solvable_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def solvable_4_11(F, trian = False):
    char = F.characteristic()
    if char != 2:
        raise ValueError("The field needs to have characteristic 2.")
    P1 = PolynomialRing(F, "t1")
    t1 = list(P1.gens())[0]
    f1 = t1**2 - 2
    K1 = NumberField(f1, 'a1')
    a1 = list(K1.gens())[0]
    P2 = PolynomialRing(K1, "t2")
    t2 = list(P2.gens())[0]
    f2 = t2**2 + 1
    K2 = NumberField(f2, 'a2')
    a2 = list(K2.gens())[0]
    d = {('x3','x0'): {'x0':1}, ('x3','x1'): {'x1':a2}, ('x3','x2'): {'x2':1+a2}, ('x2','x0'): {'x1':1}, ('x2','x1'): {'x0':a1}}
    L = LieAlgebra(K2, d, names='x0, x1, x2, x3')
    if trian == True:
        bTL = triangular_basis_solvable_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def solvable_4_12(F, trian = False):
    d = {('x3','x0'): {'x0':1}, ('x3','x1'): {'x1':2}, ('x3','x2'): {'x2':1}, ('x2','x0'): {'x1':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3')
    if trian == True:
        bTL = triangular_basis_solvable_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def solvable_4_13(F, trian = False):
    P = PolynomialRing(F, "t")
    t = list(P.gens())[0]
    f = t**2 + 1
    K = NumberField(f, 'a')
    a = list(K.gens())[0]
    d = {('x3','x0'): {'x0':1, 'x2':a}, ('x3','x1'): {'x1':1}, ('x3','x2'): {'x0':1}, ('x2','x0'): {'x1':1}}
    L = LieAlgebra(K, d, names='x0, x1, x2, x3')
    if trian == True:
        bTL = triangular_basis_solvable_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def solvable_4_14(F, trian = False):
    P = PolynomialRing(F, "t")
    t = list(P.gens())[0]
    f = t**2 + 1
    K = NumberField(f, 'a')
    a = list(K.gens())[0]
    d = {('x3','x0'): {'x2':a}, ('x3','x2'): {'x0':1}, ('x2','x0'): {'x1':1}}
    L = LieAlgebra(K, d, names='x0, x1, x2, x3')
    if trian == True:
        bTL = triangular_basis_solvable_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def solvable_lie_algebra(F, param, trian = False):
    e = str(param[0]) + '_' + str(param[1])
    return globals()['solvable_' + e](F, trian)
#-------------

def comp_solv(p1,p2):
    l = solvable_lie_algebra(QQ,p1)
    l2 = solvable_lie_algebra2(QQ,p2)
    bl = l.basis().list()
    bl2 = l2.basis().list()
    m1 = structure_constants(l,bl)
    m2 = structure_constants(l2,bl2)
    sys.displayhook(m1)
    print("------")
    sys.displayhook(m2)
