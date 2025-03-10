#-------------
def invar_test_nilp_lie_alg(L, phi):
    r'''
        INPUT:
        
            - 
            
        OUTPUT:
            
            - 

        COMMENT:


            
        EXAMPLES:
        
            

    '''
    bL = L.basis().list()
    M = structure_constants(L, bL)
    rank_M = M.rank()
    dimL = len(bL)
    print("")
    print("A dimensão da álgebra de Lie L é ", dimL, ".", sep="")
    print("O posto da matriz das constantes estruturas de L é ", rank_M, ".", sep="")
    print("Assim, o número de invariantes algebricamente independentes (?) é", end=": ")
    print(dimL, " - ", rank_M ," = ", dimL - rank_M, ".", sep="")
    gens_domain_phi = phi.domain().gens()
    print("Tais invariantes são:")
    for i in range(len(gens_domain_phi)):
        print(phi(gens_domain_phi[i]))
    print("")
    dimL = len(bL)
    d = [0]*dimL
    print("As derivações são", end=": \n")
    for i in range(dimL):
        d[i] = differential_operator(L, bL[i])
        print("d", i, " = ", d[i], sep="")
    print("")
    table = [[]]*(dimL + 1)
    table[0] = [" "]
    for i in range(1, dimL + 1):
        table[i] = ["d"+str(i-1)]
    for i in range(len(gens_domain_phi)):
        table[0] = table[0] + [phi(gens_domain_phi[i])]
    for i in range(len(gens_domain_phi)):
        for j in range(1, dimL + 1):
            table[j] = table[j] + [d[j-1](phi(gens_domain_phi[i]))]
    print("Os valores das derivações aplicadas nos invariantes são apresentados na seguinte tabela:")
    for i in range(dimL + 1):
        for j in range(len(gens_domain_phi) + 1):
            print(table[i][j], end="\t ")
        print("")
#-------------

#-------------
def polynomial_integral(pol):
    P = pol.parent()
    keys_pol = list(pol.dict().keys())
    values_pol = list(pol.dict().values())
    for i in range(len(keys_pol)):
        keys_pol[i] = keys_pol[i] + 1
    for i in range(len(keys_pol)):
        values_pol[i] = values_pol[i]/keys_pol[i]
    dict_pol = {}
    for i in range(len(keys_pol)):
        dict_pol[keys_pol[i]] = values_pol[i]
    int_pol = P(dict_pol)
    return int_pol
#-------------

#-------------
def method_characteristics_local_nilpotent(diff):
    var_bool, list_org, first_not_zero = is_local_nilpotent_derivation(diff)
    if var_bool == False:
        raise ValueError("The derivation is not locally nilpotent.")
    D = diff.parent()
    bD = list(D.gens())
    F = FractionField(D.base())
    bF = list(F.gens())
    P = F.base()
    R = P.base()
    c = 0
    if len(bF) != len(bD):
        c = len(bD) - len(bF)
    diff_l = diff.list()
    for i in range(c):
        diff_l.pop(i)
    if first_not_zero == -1:
        return bF
    for i in range(len(bF)):
        bF[i] = P(bF[i])
    curve = [0]*len(bF)
    T = PolynomialRing(F, "t")
    t = T.gens()[0]
    #center = []
    #for i in range(first_not_zero):
    #    center = center + [bF[list_org[i]]]
    for i in range(first_not_zero):
        curve[list_org[i]] = diff_l[list_org[i]]*t + T(bF[list_org[i]])
    q = -bF[list_org[first_not_zero]]/diff_l[list_org[first_not_zero]]
    for i in range(first_not_zero,len(bF)):
        aux_c = diff_l[list_org[i]]
        der_aux_c = T(aux_c.subs({aux_c.parent().gens()[j] : curve[j] for j in range(len(bF))}))
        curve[list_org[i]] = polynomial_integral(der_aux_c) + bF[list_org[i]]
    inv = [0]*(len(bF)-1)
    curve.pop(list_org[first_not_zero])
    for i in range(len(bF)-1):
        d = curve[i].degree()
        inv[i] = P(curve[i].subs({t:q})*diff_l[list_org[first_not_zero]]**d)
        inv[i] = inv[i]/inv[i].content()
        #fact = inv[i].factor()
        #if len(fact) > 1:
        #    for j in range(len(fact)):
        #        var_aux = 0
        #        for k in range(len(fact[j][0].variables())):
        #            if fact[j][0].variables()[k] not in center:
        #                var_aux = 1
        #        if var_aux == 0:
        #            inv[i] = inv[i]/fact[j][0]
        inv[i] = F(inv[i])
    return inv
#-------------

#-------------
def method_characteristics_nilpotent(d, phi = 0):
    domain_d = d.domain()
    frac_domain_d = FractionField(domain_d)
    ring_domain_d = frac_domain_d.ring()
    emb = frac_domain_d.coerce_map_from(ring_domain_d).section()
    gens = [0]*len(frac_domain_d.gens())
    for i in range(len(gens)):
        gens[i] = emb(frac_domain_d.gens()[i])
        #gens[i] = domain_d.gens()[i]
    if phi != 0:
        gens_domain_phi = phi.domain().gens()
        gens = [0]*len(gens_domain_phi)
        for i in range(len(gens_domain_phi)):
            f = phi(gens_domain_phi[i])
            f.reduce()
            gens[i] = emb(f.numerator())
            #gens[i] = f
    len_gens = len(gens)
    pols = [0]*len_gens
    curve = [0]*len_gens
    inicial_value = [0]*(len_gens - 1)
    S = PolynomialRing(QQ, len_gens - 1, "y")
    FracS = FractionField(S)
    HomFracSR = Hom(FracS,frac_domain_d)
    for i in range(len_gens):
        f = d(gens[i])
        f.reduce()
        pols[i] = f
    print(pols)
    if pols[0] != 0:
        return False
    first_not_zero = 0
    for i in range(len_gens):
        if pols[i] == 0:
            first_not_zero = first_not_zero + 1
        else:
            break
    if first_not_zero == len_gens:
        S = PolynomialRing(QQ, len_gens, "y")
        FracS = FractionField(S)
        HomFracSR = Hom(FracS,frac_domain_d)
        phi = HomFracSR(gens)
        return phi
    for i in range(first_not_zero):
        curve[i] = gens[i]
        inicial_value[i] = curve[i]
    if first_not_zero == len_gens - 1:
        phi = HomFracSR(inicial_value)
        return phi
    P = PolynomialRing(frac_domain_d, "t")
    t = P.gens()[0]
    curve[first_not_zero] = pols[first_not_zero]*t
    q = gens[first_not_zero]/pols[first_not_zero]
    for i in range(first_not_zero + 1, len_gens):
        print(i)
        if pols[i] == 0:
            curve[i] = gens[i]
            inicial_value[i-1] = curve[i]
        else:
            list_pols = [0]*i
            for j in range(i):
                list_pols[j] = gens[j]
            print( "gens are ", list_pols, "pol is", pols[i] )
            is_el = is_element_of_subalgebra(list_pols,pols[i])
            if is_el[0] == False:
                return False
            aux_c = is_el[1][0]
            der_aux_c = P(aux_c.subs({aux_c.parent().gens()[j] : curve[j] for j in range(i)}))
            curve_without_const = polynomial_integral(der_aux_c)
            inicial_value[i-1] = gens[i] - curve_without_const.subs({t:q})
            curve[i] = curve_without_const + inicial_value[i-1]
    for i in range(len_gens - 1):
        inicial_value[i] = inicial_value[i].numerator()
        inicial_value[i] = inicial_value[i]*inicial_value[i].denominator()
        inicial_value[i] = inicial_value[i]/1
    for i in range(1,len_gens - 1):
        for j in range(i):
                if inicial_value[j].divides(inicial_value[i]):
                    inicial_value[i] = inicial_value[i]/inicial_value[j]
                    inicial_value[i] = inicial_value[i].numerator()
                    inicial_value[i] = inicial_value[i]*inicial_value[i].denominator()
    phi = HomFracSR(inicial_value)
    return phi
#-------------

#-------------
def method_characteristics_diagonal(d, phi = 0):
    domain_d = d.domain()
    frac_domain_d = FractionField(domain_d)
    ring_domain_d = frac_domain_d.ring()
    gens = frac_domain_d.gens()
    if phi != 0:
        gens = phi
    coeff = [d(gens[i]) for i in range(len(gens))]
    first_not_zero = 0
    aux_bool = True
    while aux_bool:
        if coeff[first_not_zero] == 0:
            first_not_zero = first_not_zero + 1
        else:
            aux_bool = False
        if first_not_zero == len(coeff):
            return gens
    a = ring_domain_d(coeff[first_not_zero].numerator()).coefficients()[0]
    asign = a.sign()
    aabs = a.abs()
    s = gens[first_not_zero]**asign
    inv = [gens[i] for i in range(first_not_zero)]
    for i in range(first_not_zero + 1, len(gens)):
        if coeff[i] == 0:
            inv = inv + [gens[i]]
        else:
            a = ring_domain_d(coeff[i].numerator()).coefficients()[0]
            inv = inv + [gens[i]/s**(a/aabs)]
    return inv
#-------------

#-------------
def generators_algebra_rational_invariants2(L, needs_basis_change=True):
    bL = L.basis().list()

    if needs_basis_change:
        bEspL = triangular_basis_lie_algebra(L)
        Lesp = base_change_lie_algebra(L, bEspL)
    else:
        bEspL = L.basis().list()
        Lesp = L 

    #bEspLesp = triangular_basis_lie_algebra(Lesp)
    bEspLesp = Lesp.basis().list()
    dimL = len(bL)
    F = Lesp.base_ring()
    P = PolynomialRing( F, dimL, list(Lesp.basis().keys()) )
    Lesp.polynomialRing = P
    Lesp.fractionField = P.fraction_field()
    x = Lesp.polynomialRing.gens()

    if Lesp.center().dimension() == dimL:
        HomFracSR = Hom(Lesp.fractionField,Lesp.fractionField)
        phi = HomFracSR.identity()
        return [phi(phi.domain().gens()[i]) for i in range(len(phi.domain().gens()))]

    first_not_center = 0
    while Lesp.gens()[first_not_center] in Lesp.center():
        first_not_center = first_not_center + 1

    d = differential_operator(Lesp, bEspLesp[first_not_center])
    phi = method_characteristics_nilpotent(d)
    if phi == False:
        return False
    for i in range(first_not_center+1, dimL):
        d = differential_operator(Lesp, bEspLesp[i])
        print( i, d )
        #phi0 = phi
        phi = method_characteristics_nilpotent(d, phi)
        #print( "succeed", d )
        if phi == False:
            return False#, phi0, d

    FracS = phi.domain()
    gens_domain_phi = phi.domain().gens()
    gens_codomain_phi = phi.codomain().gens()
    HomFracSR = Hom(FracS,Lesp.fractionField)
    HR = [0]*(len(gens_domain_phi))
    for i in range(len(gens_domain_phi)):
        HR[i] = phi(gens_domain_phi[i]).subs({gens_codomain_phi[k] : x[k] for k in range(len(gens_codomain_phi))})
    
    y = [0]*dimL
    for i in range(dimL):
        for j in range(dimL):
            y[i] = y[i] + coord_base(L, bL, bEspL[i])[j]*x[j]
    
    for i in range(len(gens_domain_phi)):
        HR[i] = HR[i].subs({P.gens()[k] : y[k] for k in range(dimL)})
    
    
    L2 = base_change_lie_algebra(L, bL)
    bEspL2 = triangular_basis_lie_algebra(L2)
    xstr1 = [str(i) for i in bEspL]
    xstr2 = [str(i) for i in bEspL2]
    Pol1 = PolynomialRing( F, dimL, xstr1 )
    Pol2 = PolynomialRing( F, dimL, xstr2 )
    HomFracPol2Pol1 = Hom(FractionField(Pol2),FractionField(Pol1))
    HR2 = [0]*(dimL)
    for i in range(dimL):
        HR2[i] = Pol1.gens()[i]
    phi = HomFracPol2Pol1(HR2)
    for i in range(len(gens_domain_phi)):
        HR[i] = phi(HR[i])
    
    P3 = PolynomialRing( F, dimL, list(L.basis().keys()))
    HomFracSR = Hom(FracS,FractionField(P3))
    alpha = HomFracSR(HR)
    return [alpha(alpha.domain().gens()[i]) for i in range(len(alpha.domain().gens()))]
#-------------

#-------------
def generators_algebra_rational_invariants_nilpotent(l):
    bl_tri = triangular_basis_nilpotent_lie_algebra(l)
    l = base_change_lie_algebra(l, bl_tri)
    der = derivations_associated_with_base(l)
    F = der[0].parent().base()
    P = F.base()
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
    #print("inv")
    #print(inv)
    ideal_center = P.ideal(center)
    for i in range(len(center)+1,len(der)):
        #print("início")
        #print(i)
        #print("der")
        #print(der[i])
        diff, dicio = get_derivation_on_generator(der[i],inv)
        #print("diff")
        #print(diff)
        #print("first")
        if diff == 0:
            print(0)
        else:
            first = 0
            dl = diff.list()
            for j in range(len(dl)):
                if dl[j] != 0:
                    first = dl[j].subs(dicio)
                    break
            #print(first)
        inv_aux = method_characteristics_local_nilpotent(diff)
        inv = [inv_aux[j].subs(dicio) for j in range(len(inv_aux))]
        for j in range(len(inv)):
            if P(inv[j].denominator()) in ideal_center:
                inv[j] = inv[j].numerator()
            fact = inv[j].factor()
            if len(fact) > 1:
                for k in range(len(fact)):
                    if P(fact[k][0]) in ideal_center:
                        inv[j] = inv[j]/fact[k][0]
        #print("inv")
        #print(inv)
        #print("fim")
        #print("-----")
    return inv
#-------------

#-------------
def generators_algebra_rational_invariants(L):
    if L.is_nilpotent() == True:
        return generators_algebra_rational_invariants_nilpotent(L)
    #bEspL = triangular_basis_solvable_lie_algebra(L)
    #l = base_change_lie_algebra(L, bEspL)
    l=L
    cons_str_rref = structure_constants(l).rref()
    ident = identity_matrix(cons_str_rref[0,0].parent(), cons_str_rref.ncols())
    if cons_str_rref == ident:
        return []
    der = derivations_associated_with_base(l)
    inv = invariants_matrix_derivation(der[0])
    #inv = simplify_invariants(inv)
    for i in range(1,len(der)):
        diff, dict_diff = get_derivation_on_generator(der[i], inv)
        if len(diff.parent().gens()) == 1 and len(diff.list()) == 1 and diff.list()[0] != 0:
            return []
        if is_linear_derivation(diff) == True:
            inv_aux = invariants_matrix_derivation(diff)
            inv = [inv_aux[j].subs(dict_diff) for j in range(len(inv_aux))]
            #inv = simplify_invariants(inv)
        else:
            var_bool, list_org, first_not_zero= is_local_nilpotent_derivation(diff)
            if var_bool == False:
                raise ValueError("The derivation is neither linear nor locally nilpotent.")
            else:
                inv_aux = method_characteristics_local_nilpotent(diff)
                inv = [inv_aux[j].subs(dict_diff) for j in range(len(inv_aux))]
                #inv = simplify_invariants(inv)
    return inv
#-------------

#-------------
def invariant_field_isomorphism( L, needs_basis_change=True):
    return generators_algebra_rational_invariants2( L, 
                                needs_basis_change = needs_basis_change )    
#-------------
