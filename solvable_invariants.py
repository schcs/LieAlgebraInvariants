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
