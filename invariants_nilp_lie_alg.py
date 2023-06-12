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
        d_coeff = sum( coeffs[i]*F.gens()[i] for i in range( d ))
        op += d_coeff*D.gens()[i]

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
def invar_nilp_lie_alg(L):
    r'''
        INPUT:
        
            - 
            
        OUTPUT:
            
            - 

        COMMENT:


            
        EXAMPLES:
        
            

    '''
    bL = L.basis().list()
    bEspL = special_basis_nilp_lie_alg(L)
    Lesp = base_change_nilp_lie_alg(L, bEspL)
    bEspLesp = special_basis_nilp_lie_alg(Lesp)
    dimL = len(bL)
    F = Lesp.base_ring()
    P = PolynomialRing( F, dimL, 'x' )
    Lesp.polynomialRing = P
    Lesp.fractionField = P.fraction_field()
    if Lesp.center().dimension() == dimL:
        HomFracSR = Hom(Lesp.fractionField,Lesp.fractionField)
        phi = HomFracSR.identity()
        return phi
    x = Lesp.polynomialRing.gens()
    y = [0]*dimL
    for i in range(dimL):
        for j in range(dimL):
            y[i] = y[i] + coord_base(L, bL, bEspL[i])[j]*x[j]
    first_not_center = 0
    while Lesp.gens()[first_not_center] in Lesp.center():
        first_not_center = first_not_center + 1
    Sigma = adjoint_matrix_element(Lesp, Lesp.gens()[first_not_center])
    phi = invar_nilp_matrix(Sigma)
    gens_domain_phi = phi.domain().gens()
    gens_codomain_phi = phi.codomain().gens()
    for i in range(first_not_center + 1, dimL):
        d = differential_operator(Lesp,Lesp.gens()[i])
        list_pols = [0]*(len(gens_domain_phi))
        for j in range(len(gens_domain_phi)):
            list_pols[j] = phi(gens_domain_phi[j]).subs({gens_codomain_phi[k] : x[k] for k in range(len(gens_codomain_phi))})
        Sigma = matrix(QQ, len(gens_domain_phi))
        for j in range(len(gens_domain_phi)):
            pol = is_element_of_subalgebra(list_pols,d(list_pols[j]))[1][0]
            if pol.degree() > 1:
                return print("Problema com a combinação linear de polinômios")
            for k in range(len(gens_domain_phi)):
                Sigma[k,j] =  pol.coefficient(pol.parent().gens()[k])
        if Sigma != 0:
            x = [0]*(len(list_pols))
            for j in range(len(list_pols)):
                x[j] = list_pols[j]
            phi = invar_nilp_matrix(Sigma)
            gens_domain_phi = phi.domain().gens()
            gens_codomain_phi = phi.codomain().gens()
    FracS = phi.domain()
    HomFracSR = Hom(FracS,Lesp.fractionField)
    HR = [0]*(len(gens_domain_phi))
    for i in range(len(gens_domain_phi)):
        HR[i] = phi(gens_domain_phi[i]).subs({gens_codomain_phi[k] : x[k] for k in range(len(gens_codomain_phi))})
    for i in range(len(gens_domain_phi)):
        HR[i] = HR[i].subs({HR[0].parent().gens()[k] : y[k] for k in range(len(HR[0].parent().gens()))})
    phi = HomFracSR(HR)
    return phi
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
