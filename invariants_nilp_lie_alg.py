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
        P = PolynomialRing( F, d, 'x' )
        F = P.fraction_field()
        L.polynomialRing = P
        L.fractionField = F

    D = F.derivation_module()
    op = D.zero()
    bas = L.basis().list()

    for i in range( d ):
        coeffs = L.bracket( x, bas[i] ).dense_coefficient_list()
        print( coeffs )
        d_coeff = sum( coeffs[i]*F.gens()[i] for i in range( d ))
        op += d_coeff*D.gens()[i]

    return op
#-------------

def differential_operator_from_coeffs( P, coeffs ):


    D = P.derivation_module()
    op = D.zero()
    
    for k in range( len( coeffs )): 
        op += coeffs[k]*D.gens()[k]

    return op
#-------------


# Seja L uma álgebra de Lie nilpotente com base {x1, ..., xn} tal que, para todo x em L, [x, xi] é combinação linear de x1 até xi-1. Considere a derivação

# x \cdot f = \sum_{i=1}^{n}[x, xi] (derivada parcial de f com respeito a xi)

# Vamos escrever

# [x,xi] = gamma_{1i}x_{1} + ... + gamma_{i-1,i}x_{i-1}

# Defina a matrix Sigma como sendo

# 0 gamma_{12} . . . gamma_{1n} 
# 0 0 . . .          gamma_{2n}
# 0 0 . . .          gamma_{3n}
# . . . . . 
# 0 0 . . .          gamma_{n-1,n}
# 0 0 . . .          0

# A função abaixo tem como entrada a matrix Sigma e tem como saída um homomorfismo dos invariantes da derivação de x.
#-------------
def invar_nilp_matrix_jordan(Sigma):
    r'''
        INPUT:
        
            - a jordan matrix that represents a nilpotent linear transformation
            
        OUTPUT:
            
            - a homomorphism between fraction fields of polynomial rings, whose image consists of invariants of the extended derivation operator of the nilpotent linear transformation
        
        COMMENT:
        
            Given a nilpotent linear transformation `T`, consider a basis `{x_{1}, \cdots, x_{n}}` such that the matrix representing `T` in this basis is a Jordan matrix. We can extend `T` to the field of fractions of `F[x_{1}, \cdots, x_{n}]` by linearity and the Leibniz rule to obtain a derivation. Considering the matrix of `T` with respect to the fixed basis, this function yields a homomorphism whose image consists of the invariants under the extended derivation `T`.
        
        EXAMPLES:
        
            sage: L = nilpotent_lie_algebra(QQ, [6,16], True)

            sage: x = L.random_element()
            
            sage: x
            
            13*x4
            
            sage: Sigma = adjoint_matrix_element(L,x)

            sage: invar_nilp_matrix_jordan(Sigma)
            
            Ring morphism:
            From: Fraction Field of Multivariate Polynomial Ring in y0, y1, y2, y3, y4 over Rational Field
            To:   Fraction Field of Multivariate Polynomial Ring in x0, x1, x2, x3, x4, x5 over Rational Field
            Defn: y0 |--> x0
                  y1 |--> x1
                  y2 |--> (-x2^2 + 2*x1*x3)/(2*x1)
                  y3 |--> x4
                  y4 |--> (x2^3 - 3*x1*x2*x3 + 3*x1^2*x5)/(3*x1^2)

    '''
    n = Sigma.nrows() # dimensão de L
    R = PolynomialRing(QQ, n, "x") # anel de polinômios em n variáveis
    x = R.gens() # variáveis de R
    X = matrix(R, n, 1) # matriz coluna formada de zeros
    # atribuir a X as variáveis de R
    for i in range(n):
        X[i] = x[i]
    S = PolynomialRing(QQ, n - 1, "y") # anel de polinômios em n - 1 variáveis
    y = S.gens() # variáveis de R
    FracR = FractionField(R) # corpo de frações de R
    FracS = FractionField(S) # corpo de frações de S
    logaux = True # variável lógica auxiliar
    vColNZero = 0 # coluna na qual aparece o primeiro elemento não nulo de Sigma
    # determinar quem é vColNZero
    while logaux:
        vColNZero = vColNZero + 1
        for i in range(vColNZero):
            if Sigma[i,vColNZero] != 0:
                logaux = False
                i = vColNZero
    qDen = 0
    for i in range(vColNZero):
        qDen = qDen + Sigma[i,vColNZero]*x[i]
    q = x[vColNZero]/qDen
    M = -q*Sigma
    E = M.exp() # exponencial da matrix M 
    idM = matrix.identity(n)
    HR = [0]*(n - 1) # lista de zeros de tamanho n - 1
    # preenchendo HR com os valores de H. Note que HR é o H, sem o valor da posição vColNZero
    for i in range(n - 1):
        if i < vColNZero:
            for j in range(n):
                HR[i] = HR[i] + (E*idM[:,i])[j]*X[j]
        else:
            for j in range(n):
                HR[i] = HR[i] + (E*idM[:,i+1])[j]*X[j]
    HomFracSR = Hom(FracS, FracR) # conjunto dos homomorfismo de FracS para FracR
    psi = HomFracSR(HR) # homomorfismo desejado
    return psi
#-------------

#-------------
def invar_nilp_matrix(Sigma):
    r'''
        INPUT:
        
            - a matrix that represents a nilpotent linear transformation

        OUTPUT:
            
            - a homomorphism between fraction fields of polynomial rings, whose image consists of invariants of the extended derivation operator of the nilpotent linear transformation

        COMMENT:
        
            Given a nilpotent linear transformation `T` and a fixed basis of the vector space `{x_{1}, \cdots, x_{n}}`, we can extend `T` to the field of fractions of `F[x_{1}, \cdots, x_{n}]` by linearity and the Leibniz rule to obtain a derivation. Considering the matrix of `T` with respect to the fixed basis, this function yields a homomorphism whose image consists of the invariants under the extended derivation `T`.

        EXAMPLES:
        
            sage: M = matrix([[2,2,2,-3], [6,1,1,-4], [1,6,1,-4], [1,1,6,-4]])

            sage: invar_nilp_matrix(M)

            Ring morphism:
                
                From: Fraction Field of Multivariate Polynomial Ring in y0, y1, y2 over Rational Field
                
                To:   Fraction Field of Multivariate Polynomial Ring in x0, x1, x2, x3 over Rational Field
  
                Defn: y0 |--> 100*x0 + 100*x1 + 100*x2 + 200*x3
                      y1 |--> (175*x0^2 + 1150*x0*x1 + 975*x1^2 - 450*x0*x2 + 350*x1*x2 - 1025*x2^2 + 700*x0*x3 + 2300*x1*x3 - 100*x2*x3 + 300*x3^2)/(200*x0 + 200*x1 + 200*x2 + 400*x3)
                      y2 |--> (24375*x0^3 + 25125*x0^2*x1 - 22875*x0*x1^2 - 23625*x1^3 + 49125*x0^2*x2 - 21750*x0*x1*x2 - 70875*x1^2*x2 + 49125*x0*x2^2 - 22875*x1*x2^2 + 32375*x2^3 + 98250*x0^2*x3 + 28500*x0*x1*x3 - 69750*x1^2*x3 + 76500*x0*x2*x3 - 139500*x1*x2*x3 + 2250*x2^2*x3 + 100500*x0*x3^2 - 43500*x1*x3^2 - 19500*x2*x3^2 - 5000*x3^3)/(30000*x0^2 + 60000*x0*x1 + 30000*x1^2 + 60000*x0*x2 + 60000*x1*x2 + 30000*x2^2 + 120000*x0*x3 + 120000*x1*x3 + 120000*x2*x3 + 120000*x3^2)
        
    '''
    n = Sigma.nrows()
    log_value_aux = True
    for i in range(n+1):
        if Sigma**n == 0:
            i = n
            log_value_aux = False
    if log_value_aux == True:
        return print("The matrix is not nilpotent.")
    J, P = Sigma.jordan_form(transformation=True, subdivide=False)
    phi = invar_nilp_matrix_jordan(J)
    HomFracSR = phi.parent()
    FracS = phi.domain()
    FracR = phi.codomain()
    x = FracR.gens() 
    X = matrix(FracR, n, 1)
    for i in range(n):
        X[i] = x[i]
    HR = [0]*(n - 1)
    y = [0]*n
    for i in range(n):
        for j in range(n):
            y[i] = y[i] + P[:,i][j]*X[j]
    for i in range(n - 1):
        HR[i] = phi(FracS.gens()[i]).subs({ FracR.gens()[j]: y[j] for j in range(n)})
    psi = HomFracSR(HR) 
    return psi
#-------------

#-------------
def invar_nilp_lie_alg( L ):
    #L = nilpotent_lie_algebra(F, args, True)
    dimL = len(L.basis())
    F = L.base_ring()
    P = PolynomialRing( F, dimL, 'x' )
    F = P.fraction_field()
    L.polynomialRing = P
    L.fractionField = F
    x = L.polynomialRing.gens()
    first_not_center = 0
    while L.gens()[first_not_center] in L.center():
        first_not_center = first_not_center + 1
    Sigma = adjoint_matrix_element(L, L.gens()[first_not_center])
    phi = invar_nilp_matrix_jordan(Sigma)
    gens_domain_phi = phi.domain().gens()
    gens_codomain_phi = phi.codomain().gens()
    for i in range(first_not_center + 1, dimL):
        d = differential_operator(L,L.gens()[i])
        list_pols = [0]*(len(gens_domain_phi))
        for j in range(len(gens_domain_phi)):
            list_pols[j] = phi(gens_domain_phi[j]).subs({gens_codomain_phi[k] : x[k] for k in range(len(gens_codomain_phi))})

        Sigma = matrix(QQ, len(gens_domain_phi))
        for j in range(len(gens_domain_phi)):
            pol = is_element_of_subalgebra(list_pols,d(list_pols[j]))[1][0]
            if pol.degree() > 1:
                return("Problema com a combinação linear de polinômios")
            for k in range(len(gens_domain_phi)):
                Sigma[k,j] =  pol.coefficient(pol.parent().gens()[k])
        if Sigma != 0:
            x = [0]*(len(list_pols))
            for j in range(len(list_pols)):
                x[j] = list_pols[j]
            phi = invar_nilp_matrix_jordan(Sigma)
            gens_domain_phi = phi.domain().gens()
            gens_codomain_phi = phi.codomain().gens()
    FracS = phi.domain()
    HomFracSR = Hom(FracS,Lesp.fractionField)
    HR = [0]*(len(gens_domain_phi))
    for i in range(len(gens_domain_phi)):
        HR[i] = phi(gens_domain_phi[i]).subs({gens_codomain_phi[k] : x[k] for k in range(len(gens_codomain_phi))})
    for i in range(len(gens_domain_phi)):
        HR[i] = HR[i].subs({HR[0].parent().gens()[k] : y[k] for k in range(len(HR[0].parent().gens()))})
    alpha = HomFracSR(HR)
    return alpha
#-------------

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
def method_characteristics_simple(d, phi = 0):
    domain_d = d.domain()
    ring_domain_d = domain_d.ring()
    emb = domain_d.coerce_map_from(ring_domain_d).section()
    frac_domain_d = FractionField(domain_d)
    gens = [0]*len(domain_d.gens())
    for i in range(len(gens)):
        gens[i] = emb(domain_d.gens()[i])
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
        if pols[i] == 0:
            curve[i] = gens[i]
            inicial_value[i-1] = curve[i]
        else:
            list_pols = [0]*i
            for j in range(i):
                list_pols[j] = gens[j]
            if is_element_of_subalgebra(list_pols,pols[i])[0] == False:
                return False
            aux_c = is_element_of_subalgebra(list_pols,pols[i])[1][0]
            der_aux_c = P(aux_c.subs({aux_c.parent().gens()[j] : curve[j] for j in range(i)}))
            curve_without_const = polynomial_integral(der_aux_c)
            inicial_value[i-1] = (gens[i] - curve_without_const.subs({t:q})).numerator()
            inicial_value[i-1] = inicial_value[i-1]*inicial_value[i-1].denominator()
            for j in range(i-1):
                if inicial_value[j].divides(inicial_value[i-1]):
                    inicial_value[i-1] = inicial_value[i-1]/inicial_value[j]
                    inicial_value[i-1] = inicial_value[i-1].numerator()
                    inicial_value[i-1] = inicial_value[i-1]*inicial_value[i-1].denominator()
            curve[i] = curve_without_const + inicial_value[i-1]
    phi = HomFracSR(inicial_value)
    return phi
#-------------

#-------------
def invar_nilp_lie_alg_via_method_characteristics_simple(L, needs_basis_change = true ):
    bL = L.basis().list()

    if needs_basis_change:
        bEspL = special_basis_nilp_lie_alg(L)
        Lesp = base_change_nilp_lie_alg(L, bEspL)
    else:
        bEspL = L.basis().list()
        Lesp = L 

    #bEspLesp = special_basis_nilp_lie_alg(Lesp)
    bEspLesp = Lesp.basis().list()
    dimL = len(bL)
    F = Lesp.base_ring()
    P = PolynomialRing( F, dimL, Lesp.basis().keys().list() )
    Lesp.polynomialRing = P
    Lesp.fractionField = P.fraction_field()
    x = Lesp.polynomialRing.gens()

    if Lesp.center().dimension() == dimL:
        HomFracSR = Hom(Lesp.fractionField,Lesp.fractionField)
        phi = HomFracSR.identity()
        return phi

    first_not_center = 0
    while Lesp.gens()[first_not_center] in Lesp.center():
        first_not_center = first_not_center + 1

    d = differential_operator(Lesp, bEspLesp[first_not_center])
    phi = method_characteristics_simple(d)
    if phi == False:
        return False
    for i in range(first_not_center + 1, dimL):
        d = differential_operator(Lesp, bEspLesp[i])
        print( i, d )
        #phi0 = phi
        phi = method_characteristics_simple(d, phi)
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
    
    
    L2 = base_change_nilp_lie_alg(L, bL)
    bEspL2 = special_basis_nilp_lie_alg(L2)
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
    
    P3 = PolynomialRing( F, dimL, L.basis().keys().list() )
    HomFracSR = Hom(FracS,FractionField(P3))
    alpha = HomFracSR(HR)
    return alpha

def invariant_field_isomorphism( L, needs_basis_change = true ):
    return invar_nilp_lie_alg_via_method_characteristics_simple( L, 
                                needs_basis_change = needs_basis_change )    
#-------------

#-x0*x3*x4*x5 + x0*x1*x5*x7 + x0*x2*x3*x8 - x0^2*x7*x8 - x0*x1*x2*x12 + x0^2*x4*x12
