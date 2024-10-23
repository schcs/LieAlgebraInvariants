#-------------
def nilpotent_1_1(F, trian = False):
    L = lie_algebras.abelian(F, names='x0')
    return L
#-------------

#-------------
def nilpotent_2_1(F, trian = False):
    L = lie_algebras.abelian(F, names='x0, x1')
    return L
#-------------

#-------------
def nilpotent_3_1(F, trian = False):
    L = lie_algebras.abelian(F, names='x0, x1, x2')
    return L
#-------------

#-------------
def nilpotent_3_2(F, trian = False):
    d = {('x0','x1'): {'x2':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_4_1(F, trian = False):
    L = lie_algebras.abelian(F, names='x0, x1, x2, x3')
    return L
#-------------

#-------------
def nilpotent_4_2(F, trian = False):
    d = {('x0','x1'): {'x2':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_4_3(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x3':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_5_1(F, trian = False):
    L = lie_algebras.abelian(F, names='x0, x1, x2, x3, x4')
    return L
#-------------

#-------------
def nilpotent_5_2(F, trian = False):
    d = {('x0','x1'): {'x2':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_5_3(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x3':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_5_4(F, trian = False):
    d = {('x0','x1'): {'x4':1}, ('x2','x3'): {'x4':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_5_5(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x4':1}, ('x1','x3'): {'x4':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_5_6(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x3':1}, ('x0','x3'): {'x4':1}, ('x1','x2'): {'x4':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_5_7(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x3':1}, ('x0','x3'): {'x4':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_5_8(F, trian = False):
    d = {('x0','x1'): {'x3':1}, ('x0','x2'): {'x4':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_5_9(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x3':1}, ('x1','x2'): {'x4':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_1(F, trian = False):
    L = lie_algebras.abelian(F, names='x0, x1, x2, x3, x4, x5')
    return L
#-------------

#-------------
def nilpotent_6_2(F, trian = False):
    d = {('x0','x1'): {'x2':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_3(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x3':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_4(F, trian = False):
    d = {('x0','x1'): {'x4':1}, ('x2','x3'): {'x4':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_5(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x4':1}, ('x1','x3'): {'x4':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_6(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x3':1}, ('x0','x3'): {'x4':1}, ('x1','x2'): {'x4':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_7(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x3':1}, ('x0','x3'): {'x4':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_8(F, trian = False):
    d = {('x0','x1'): {'x3':1}, ('x0','x2'): {'x4':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_9(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x3':1}, ('x1','x2'): {'x4':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_10(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x5':1}, ('x3','x4'): {'x5':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_11(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x3':1}, ('x0','x3'): {'x5':1}, ('x1','x2'): {'x5':1}, ('x1','x4'): {'x5':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_12(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x3':1}, ('x0','x3'): {'x5':1}, ('x1','x4'): {'x5':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_13(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x4':1}, ('x1','x3'): {'x4':1}, ('x0','x4'): {'x5':1}, ('x2','x3'): {'x5':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_14(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x3':1}, ('x0','x3'): {'x4':1}, ('x1','x2'): {'x4':1}, ('x1','x4'): {'x5':1}, ('x2','x3'): {'x5':-1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_15(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x3':1}, ('x0','x3'): {'x4':1}, ('x1','x2'): {'x4':1}, ('x0','x4'): {'x5':1}, ('x1','x3'): {'x5':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_16(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x3':1}, ('x0','x3'): {'x4':1}, ('x1','x4'): {'x5':1}, ('x2','x3'): {'x5':-1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_17(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x3':1}, ('x0','x3'): {'x4':1}, ('x0','x4'): {'x5':1}, ('x1','x2'): {'x5':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_18(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x3':1}, ('x0','x3'): {'x4':1}, ('x0','x4'): {'x5':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_19(F, trian = False):
    P = PolynomialRing(F, "t")
    t = list(P.gens())[0]
    f = t**2 + 1
    K = NumberField(f, 'a')
    a = list(K.gens())[0]
    d = {('x0','x1'): {'x3':1}, ('x0','x2'): {'x4':1}, ('x0','x4'): {'x5':1}, ('x1','x3'): {'x5':1}, ('x2','x4'): {'x5':a}}
    L = LieAlgebra(K, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_20(F, trian = False):
    d = {('x0','x1'): {'x3':1}, ('x0','x2'): {'x4':1}, ('x0','x4'): {'x5':1}, ('x1','x3'): {'x5':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_21(F, trian = False):
    P = PolynomialRing(F, "t")
    t = list(P.gens())[0]
    f = t**2 + 1
    K = NumberField(f, 'a')
    a = list(K.gens())[0]
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
    P = PolynomialRing(F, "t")
    t = list(P.gens())[0]
    f = t**2 + 1
    K = NumberField(f, 'a')
    a = list(K.gens())[0]
    d = {('x0','x1'): {'x4':1}, ('x0','x2'): {'x5':1}, ('x1','x3'): {'x5':a}, ('x2','x3'): {'x4':1}}
    L = LieAlgebra(K, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_23(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x4':1}, ('x0','x3'): {'x5':1}, ('x1','x3'): {'x4':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_24(F, trian = False):
    P = PolynomialRing(F, "t")
    t = list(P.gens())[0]
    f = t**2 + 1
    K = NumberField(f, 'a')
    a = list(K.gens())[0]
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x4':1}, ('x0','x3'): {'x5':a}, ('x1','x2'): {'x5':1}, ('x1','x3'): {'x4':1}}
    L = LieAlgebra(K, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_25(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x4':1}, ('x0','x3'): {'x5':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_6_26(F, trian = False):
    d = {('x0','x1'): {'x2':1}, ('x0','x2'): {'x4':1}, ('x1','x3'): {'x5':1}}
    L = LieAlgebra(F, d, names='x0, x1, x2, x3, x4, x5')
    if trian == True:
        bTL = triangular_basis_nilpotent_lie_algebra(L)
        l_iso = base_change_lie_algebra(L, bTL)
        return l_iso
    return L
#-------------

#-------------
def nilpotent_lie_algebra(F, param, trian = False):
    e = str(param[0]) + '_' + str(param[1])
    return globals()['nilpotent_' + e](F, trian)
#-------------

def comp_nilp(p1,p2):
    l = nilpotent_lie_algebra(QQ,p1)
    l2 = nilpotent_lie_algebra(QQ,p2)
    bl = l.basis().list()
    bl2 = l2.basis().list()
    m1 = structure_constants(l,bl)
    m2 = structure_constants(l2,bl2)
    sys.displayhook(m1)
    print("------")
    sys.displayhook(m2)
