import sys

#-------------
def derivations_associated_with_base(L):
    bL = L.basis().list()
    S = structure_constants(L,bL)
    F = S[0,0].parent()
    D = F.derivation_module()
    bD_comp = D.basis().list()
    bD = bD_comp[len(bD_comp)-len(bL):len(bD_comp)] # Elimina os possíveis elementos do corpo base de L
    der = []
    for i in range(len(bL)):
        z = D.zero()
        for j in range(len(bL)):
            z = z + S[i,j]*bD[j]
        der = der + [z]
    return der
#-------------

#-------------
def matrices_associated_with_base(L):
    der = derivations_associated_with_base(L)
    mat = []
    for i in range(len(der)):
        mat = mat + [matrix_of_derivation(der[i])]
    return mat
#-------------

#-------------
def is_linear_derivation(diff):
    coeff = diff.list()
    for i in range(len(coeff)):
        if coeff[i].denominator() != 1 or coeff[i].numerator().degree() not in [-1,0,1]:
            return False
    return True
#-------------

#-------------
def is_local_nilpotent_derivation(diff):
    D = diff.parent()
    bD = list(D.gens())
    F = D.base()
    bF = list(F.gens())
    P = F.base()
    R = P.base()
    c = 0
    if len(bF) != len(bD):
        c = len(bD) - len(bF)
    diff_l = diff.list()
    for i in range(c):
        diff_l.pop(i)
    if diff == D.zero():
        return True, [i for i in range(len(diff_l))], [-1,-1]
    list_org = []
    for i in range(len(diff_l)):
        if diff_l[i] not in P:
            return False, [], [0,0]
        else:
            diff_l[i] = P(diff_l[i])
    for i in range(len(diff_l)):
        if diff_l[i] == 0:
            list_org = list_org + [i]
    firt_non_zero = len(list_org)
    for i in range(len(diff_l)):
        if diff_l[i] in R and diff_l[i] != 0:
            list_org = list_org + [i]
    firt_non_const = len(list_org)
    if list_org == []:
        return False, [], [0,0]
    val_bool = True
    while val_bool:
        val_bool = False
        list_aux = []
        for i in range(len(diff_l)):
            if i in list_org:
                continue
            var = list(diff_l[i].variables())
            var_test = 0
            for j in range(len(var)):
                for k in range(len(list_org)):
                    if var[j] == bF[list_org[k]]:
                        var_test = var_test + 1
            if var_test == len(var):
                list_aux = list_aux + [i]
        if len(list_org) < len(diff_l) and list_aux == []:
            return False, [], [0,0]
        if list_aux != []:
            list_org = list_org + list_aux
            val_bool = True
    return True, list_org, [firt_non_zero, firt_non_const]
#-------------

#-------------
def get_derivation_on_generator(diff, param):
    coeff = [diff(param[i]) for i in range(len(param))]
    new_coeff = []
    for i in range(len(coeff)):
        v, der = is_element_of_subalgebra(param,coeff[i])
        if v == False:
            return False
        new_coeff = new_coeff + [der[0]]
    F = FractionField(new_coeff[0].parent())
    D = F.derivation_module()
    diff_alt = sum(new_coeff[i]*D.gens()[i] for i in range(len(new_coeff)))
    dict_param = {F.gens()[i]:param[i] for i in range(len(param))}
    return diff_alt, dict_param
#-------------

#-------------
def simplify_invariants(inv):
    if len(inv) in [0,1]:
        return inv
    inv_simp = []
    with_den = []
    for i in range(len(inv)):
        if inv[i].denominator() != 1:
            with_den = with_den + [inv[i]]
        else:
            inv_simp = inv_simp + [inv[i]]
    if inv_simp == []:
        return inv
    var_bool = True
    while var_bool:
        var_bool = False
        rem = []
        for i in range(len(with_den)):
            v, der = is_element_of_subalgebra(inv_simp, with_den[i].denominator())
            if v == True:
                var_bool = True
                inv_simp = inv_simp + [with_den[i].numerator()]
                rem = rem + [with_den[i]]
        with_den = [x for x in with_den if x not in rem]
    inv_simp = inv_simp + with_den
    for i in range(len(inv_simp)):
        coeff = inv_simp[i].numerator().coefficients()
        if len(coeff) == 1:
            inv_simp[i] = inv_simp[i]/coeff[0]
        else:
            var = 0
            for j in range(1,len(coeff)):
                if coeff[j] not in [coeff[0], -coeff[0]]:
                    var = 1
            if var == 0:
                inv_simp[i] = inv_simp[i]/coeff[0]
    return inv_simp
#-------------

#-------------
def generators_algebra_rational_invariants_nilpotent(l):
    if l.is_nilpotent() == False:
        raise ValueError("The algebra is not nilpotent.")
    bl_tri = triangular_basis_nilpotent_lie_algebra(l)
    l = base_change_lie_algebra(l, bl_tri)
    der = derivations_associated_with_base(l)
    F = der[0].parent().base()
    gensF = list(F.gens())
    center = [gensF[0]]
    for i in range(1,len(der)):
        if der[i] == 0:
            center = center + [gensF[i]]
        else:
            break
    if len(center) == len(der):
        return gensF
    inv = method_characteristics_local_nilpotent(der[len(center)])
    for i in range(len(center)+1,len(der)):
        diff, dict_diff = get_derivation_on_generator(der[i], inv)
        var_bool, list_org, [first_not_zero, first_not_const] = is_local_nilpotent_derivation(diff)
        if var_bool == False:
            raise ValueError("The derivation is neither linear nor locally nilpotent.")
        else:
            inv_aux = method_characteristics_local_nilpotent(diff)
            inv = [inv_aux[j].subs(dict_diff) for j in range(len(inv_aux))]
            for j in range(len(inv)):
                fact = inv[j].factor()
                if len(fact) > 1:
                    for k in range(len(fact)):
                        var_aux = 0
                        for l in range(len(fact[k][0].variables())):
                            if fact[k][0].variables()[l] not in center:
                                var_aux = 1
                        if var_aux == 0:
                            inv[j] = inv[j]/fact[k][0]
    return inv
#-------------

#-------------
def generators_algebra_rational_invariants2(L):
    bEspL = triangular_basis_lie_algebra(L)
    l = base_change_lie_algebra(L, bEspL)
    der = derivations_associated_with_base(l)
    inv = invariants_matrix_derivation(der[0])
    inv = simplify_invariants(inv)
    for i in range(1,len(der)):
        diff, dict_diff = get_derivation_on_generator(der[i], inv)
        if len(diff.parent().gens()) == 1 and len(diff.list()) == 1 and diff.list()[0] != 0:
            return []
        if is_linear_derivation(diff) == True:
            inv_aux = invariants_matrix_derivation(diff)
            inv = [inv_aux[j].subs(dict_diff) for j in range(len(inv_aux))]
            inv = simplify_invariants(inv)
        else:
            var_bool, list_org, [first_not_zero, first_not_const] = is_local_nilpotent_derivation(diff)
            if var_bool == False:
                raise ValueError("The derivation is neither linear nor locally nilpotent.")
            else:
                inv_aux = method_characteristics_local_nilpotent(diff)
                inv = [inv_aux[j].subs(dict_diff) for j in range(len(inv_aux))]
                inv = simplify_invariants(inv)
    return inv
#-------------

#-------------
def print_solvable():
    list_solv = [
    [2,1],
    [2,2],
    [3,1],
    [3,2],
    [3,3,2],
    [3,4,2],
    [4,1],
    [4,2],
    [4,3,2],
    [4,4],
    [4,5],
    [4,6,2,3],
    [4,7,2,3],
    [4,8],
    [4,9,2],
    [4,10,2],
    [4,12],
    [4,13,2],
    [4,14,2]
    ]
    for i in range(len(list_solv)):
        print("Parâmentro")
        print(list_solv[i],"\n")
        l = solvable_lie_algebra(QQ,list_solv[i])
        btri = triangular_basis_solvable_lie_algebra(l)
        L = base_change_lie_algebra(l, btri)
        print("Matrizes")
        mat = matrices_associated_with_base(L)
        sys.displayhook(mat)
        print("\nDerivações")
        der = derivations_associated_with_base(L)
        sys.displayhook(der)
        print("\nInvariantes de cada elemento da base")
        inv = [invariants_matrix_derivation(der[i]) for i in range(len(der))]
        sys.displayhook(inv)
        inv_field = [[2,2], [4,8], [4,9,2], [4,10,2], [4,12], [4,13,2]]
        inv_dif = [[4,14,2]]
        inv_todos = [[2,1], [3,1], [4,1]]
        if list_solv[i] in inv_field:
            print("\nInvariantes são o corpo")
        if list_solv[i] in inv_dif:
            print("\nPrecisa de mais conta")
        if list_solv[i] in inv_todos:
            print("\nTodo os racionais são invariantes")
        if list_solv[i] not in inv_field and list_solv[i] not in inv_dif and list_solv[i] not in inv_todos:
            M = mat[L.dimension()-1].submatrix(0, 0, L.dimension()-1, L.dimension()-1)
            print("\nMatriz a ser resolvida")
            print(M)
            print()
            f = M.characteristic_polynomial()
            #K = extension_field_roots(f)
            K = f.splitting_field("a")
            M = Matrix(K,M)
            J, P = M.jordan_form(transformation=True)
            print("Jordan")
            print(J)
            print()
            print("Invariantes da álgebra de Lie")
            diff = derivation_of_matrix(M)
            sol = invariants_matrix_derivation(diff)
            sys.displayhook(sol)
        print("------------\n")
#-------------

#-------------
def print_solvable2():
    list_solv = [
    [2,1],
    [2,2],
    [3,1],
    [3,2],
    [3,3],
    [3,4],
    [4,1],
    [4,2],
    [4,3],
    [4,4],
    [4,5],
    [4,6],
    [4,7],
    [4,8],
    [4,9],
    [4,10],
    [4,12],
    [4,13],
    [4,14]
    ]
    for i in range(len(list_solv)):
        print("Parâmentro")
        print(list_solv[i],"\n")
        l = solvable_lie_algebra2(QQ,list_solv[i])
        btri = triangular_basis_solvable_lie_algebra(l)
        L = base_change_lie_algebra(l, btri)
        print("Matrizes")
        mat = matrices_associated_with_base(L)
        sys.displayhook(mat)
        print("\nDerivações")
        der = derivations_associated_with_base(L)
        sys.displayhook(der)
        print("\nInvariantes de cada elemento da base")
        inv = [invariants_matrix_derivation(der[i]) for i in range(len(der))]
        sys.displayhook(inv)
        inv_field = [[2,2], [4,8], [4,9], [4,10], [4,12], [4,13]]
        inv_dif = [[4,14]]
        inv_todos = [[2,1], [3,1], [4,1]]
        if list_solv[i] in inv_field:
            print("\nInvariantes são o corpo")
        if list_solv[i] in inv_dif:
            print("\nPrecisa de mais conta")
        if list_solv[i] in inv_todos:
            print("\nTodo os racionais são invariantes")
        if list_solv[i] not in inv_field and list_solv[i] not in inv_dif and list_solv[i] not in inv_todos:
            M = mat[L.dimension()-1].submatrix(0, 0, L.dimension()-1, L.dimension()-1)
            print("\nMatriz a ser resolvida")
            print(M)
            print()
            f = M.characteristic_polynomial()
            K = extension_field_roots(f)
            #K = f.splitting_field("a")
            M = Matrix(K,M)
            J, P = M.jordan_form(transformation=True)
            print("Jordan")
            print(J)
            print()
            print("Invariantes da álgebra de Lie")
            diff = derivation_of_matrix(M)
            sol = invariants_matrix_derivation(diff)
            sys.displayhook(sol)
        print("------------\n")
#-------------
