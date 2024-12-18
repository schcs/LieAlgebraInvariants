#-------------
def differential_operator( L, x ):
    r'''
        INPUT:
        
        - 
            
        OUTPUT:
            
        - 
            
        COMMENT:

        Constructs the differential operators on the polynomial algebra generated by the L.basis()
        that corresponds to the adjoint action of the element `x` on `L`.
        If `L` is a Lie algebra with basis `x_1,\ldots,x_n`, and `d_{x_i}` is the partial derivative by 
        `x_i`, then the differential operator is defined as
       
        .. MATH:: 
       
            \sum_{i=1}^{\dim L} [x,x_i]d_{x_1}

        EXAMPLES:



    '''
    
    F = L.base_ring()
    d = L.dimension()
    
    if hasattr( L, "polynomialRing" ):
        P = L.polynomialRing
        F = L.fractionField
    else: 
        P = PolynomialRing( F, d, list(L.basis().keys()) )
        F = P.fraction_field()
        L.polynomialRing = P
        L.fractionField = F

    D = F.derivation_module()
    op = D.zero()
    bas = L.basis()
    
    mc = [ F(L.bracket( x, y ))*F.derivation( F( str( y ))) for y in bas ]
    return sum( mc )

    for i in range( d ):
        prod = L.bracket( x, bas[i] )
        coeffs = prod.monomial_coefficients()
        
        d_coeff = sum( coeffs[i]*F.gens()[i] for i in range( d ))
        op += d_coeff*D.gens()[i]

    return op
#-------------

#-------------
def differential_operator_from_coeffs( P, coeffs ):

    D = P.derivation_module()
    op = D.zero()
    
    for k in range( len( coeffs )): 
        op += coeffs[k]*D.gens()[k]

    return op
#-------------

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
    D = diff.parent()
    P = D.base()
    gensD = D.gens()
    gensP = P.gens()
    n = len(gensP)
    c = 0
    if n != len(gensD):
        c = len(gensD) - n
    M = Matrix(P.base_ring(), n)
    mc = diff.monomial_coefficients()

    for i in range(c,n+c):
        pol = mc[i].numerator() if i in mc.keys() else P.zero().numerator()
        mon = pol.monomials()
        coeff = pol.coefficients()
        mon_coeff_dict = dict( zip( mon, coeff ))
        for j in range(len(mon)):
            cont = gensP.index( mon[j] )
            M[cont,i-c] = mon_coeff_dict[mon[j]]

    return M
#-------------

#-------------
def derivations_associated_with_base(L):
    bL = L.basis().list()
    der = [0]*(len(bL))
    for i in range(len(bL)):
        der[i] = differential_operator(L, bL[i])
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
        if coeff[i].denominator() != 1 or coeff[i].numerator().degree() not in [-1,1]:
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
        return True, [i for i in range(len(diff_l))], -1
    list_org = []
    for i in range(len(diff_l)):
        if diff_l[i] not in P:
            return False, [], 0
        else:
            diff_l[i] = P(diff_l[i])
    for i in range(len(diff_l)):
        if diff_l[i] == 0:
            list_org = list_org + [i]
    firt_non_zero = len(list_org)
    for i in range(len(diff_l)):
        if diff_l[i] in R and diff_l[i] != 0:
            list_org = list_org + [i]
    #firt_non_const = len(list_org)
    if list_org == []:
        return False, [], 0
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
            return False, [], 0
        if list_aux != []:
            list_org = list_org + list_aux
            val_bool = True
    return True, list_org, firt_non_zero #[firt_non_zero, firt_non_const]
#-------------

#-------------
def get_derivation_on_generator(diff, param):
    coeff = [diff(param[i]) for i in range(len(param))]
    non_zero = []
    for i in range(len(coeff)):
        if coeff[i] != 0:
            non_zero = non_zero + [i]
    if non_zero == []:
        R = param[0].parent().base_ring()
        P = PolynomialRing(R, len(param), "t")
        F = FractionField(P)
        D = F.derivation_module()
        c = 0
        if len(F.gens()) != len(D.gens()):
            c = len(D.gens()) - len(F.gens())
        diff = D.zero()
        dict_param = {F.gens()[i]:param[i] for i in range(len(param))}
        return diff, dict_param
    if len(non_zero) == 1:
        print(coeff[non_zero[0]])
        R = param[0].parent().base_ring()
        P = PolynomialRing(R, len(param), "t")
        F = FractionField(P)
        D = F.derivation_module()
        c = 0
        if len(F.gens()) != len(D.gens()):
            c = len(D.gens()) - len(F.gens())
        diff = D.gens()[c+non_zero[0]] 
        dict_param = {F.gens()[i]:param[i] for i in range(len(param))}
        return diff, dict_param
    new_coeff = [0]*len(coeff)
    for i in range(len(non_zero)):
        v, der = is_element_of_subalgebra(param,coeff[non_zero[i]])
        if v == False:
            return False
        new_coeff[non_zero[i]] = der[0]
    F = FractionField(new_coeff[non_zero[0]].parent())
    for i in range(len(new_coeff)):
        new_coeff[i] = F(new_coeff[i])
    D = F.derivation_module()
    c = 0
    if len(F.gens()) != len(D.gens()):
        c = len(D.gens()) - len(F.gens()) 
    diff_alt = sum(new_coeff[i]*D.gens()[c+i] for i in range(len(new_coeff)))
    dict_param = {F.gens()[i]:param[i] for i in range(len(param))}
    return diff_alt, dict_param
#-------------
