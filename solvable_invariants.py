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
    for j in range(len(bL)):
        z = D.zero()
        for i in range(len(bL)):
            z = z + S[i,j]*bD[i]
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
    F = D.base()
    P = F.base()
    R = P.base()
    bF = list(F.gens())
    diff_l = diff.list()
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
def generators_algebra_rational_invariants2(L):
    bEspL = triangular_basis_lie_algebra(L)
    l = base_change_lie_algebra(L, bEspL)
    der = derivations_associated_with_base(l)
    inv = invariants_matrix_derivation(der[0])
    for i in range(1,len(der)):
        diff, dict_diff = get_derivation_on_generator(der[i], inv)
        if is_linear_derivation(diff) == True:
            inv_aux = invariants_matrix_derivation(diff)
            inv = [inv_aux[j].subs(dict_diff) for j in range(len(inv_aux))]
        else:
            var_bool, list_org, [first_not_zero, first_not_const] = is_local_nilpotent_derivation(diff)
            if var_bool == False:
                raise ValueError("The derivation is neither linear nor locally nilpotent.")
            else:
                inv_aux = method_characteristics_local_nilpotent(diff)
                inv = [inv_aux[j].subs(dict_diff) for j in range(len(inv_aux))]
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
