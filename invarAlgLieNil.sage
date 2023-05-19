# Seja L uma álgebra de Lie nilpotente com base {x1, ..., xn} tal que, para todo x em L, [x, xi] é combinação linear de x1 até xi-1. Considere a derivação

# x \cdot f = \sum_{i=1}^{n}[x, xi] (derivada parcial de f com respeito a xi)

# Vamos escrever

# [x,xi] = gamma_{i}^{1}x1 + ... + gamma_{i}^{i-1}xi-1

# Defina a matrix Sigma como sendo

# 0 0 0 0 ... 0 0
# gamma_{21} 0 0 0 ... 0 0
# gamma_{31} gamma_{32} 0 0 ... 0 0
# . . . . . .
# gamma_{n1} gamma_{n2} gamma_{n3} gamma_{n4} ... gamma_{n,n-1} 0

# A função abaixo tem como entrada a matrix Sigma e tem como saída um homomorfismo dos invariantes da derivação de x.

def invar(Sigma):
    '''Dada uma matriz de derivação de um elemento em uma álgebra de Lie nilpotente, retorna um homomorfismo com seus invariantes'''
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
    vLinNZero = 0 # linha na qual aparece o primeiro elemento não nulo de Sigma
    # determinar quem é vLinNZero
    while logaux:
        vLinNZero = vLinNZero + 1
        for i in range(vLinNZero):
            if Sigma[vLinNZero][i] != 0:
                logaux = False
                i = vLinNZero
    qDen = 0
    for i in range(vLinNZero):
        qDen = qDen + Sigma[vLinNZero][i]*x[i]
    q = x[vLinNZero]/qDen
    M = -q*Sigma
    E = M.exp() # exponencial da matrix M 
    H = (E*X).list() # lista da multiplicação de E com X
    HR = [0]*(n - 1) # lista de zeros de tamanho n - 1
    # preenchendo HR com os valores de H. Note que HR é o H, sem o valor da posição vLinNZero
    for i in range(n - 1):
        if i < vLinNZero:
            HR[i] = H[i]
        else:
            HR[i] = H[i + 1]
    HomFracSR = Hom(FracS, FracR) # conjunto dos homomorfismo de FracS para FracR
    psi = HomFracSR(HR) # homomorfismo desejado
    return psi
