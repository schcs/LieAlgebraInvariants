#     MOTIVAÇÃO DO PROGRAMA

# Seja L uma F-álgebra de Lie nilpotente com base {x_{1}, ..., x_{n}} tal que [x_{i}, x_{j}] é combinação linear de x_{1}, ..., x_{i-1}, se i < j. Seja F(x_{1}, ..., x_{n}) o corpo de frações do anel F[x_{1}, ..., x_{n}]. Considere a ação de L sobre F(x_{1}, ..., x_{n}) dada por
# x \cdot f = r_{21} x_{1} der(f,x_{2}) + (r_{31} x_{1} + r_{32} x_{2}) der(f, x_{3}) + ... + (r_{n1} x_{1} + r_{n2} x_{2} + ... + r_{n,n-1} x_{n-1}) der(f, x_{n-1})
# onde r_{ij} pertence a F e der(f, x_{i}) é a derivada parcial de f com respeito a x_{i}. Queremos determinar quais são os f tais que x \cdot f = 0.

# Usando o Método das Características, as frações f tais que x \cdot f = 0, são determinadas a partir de um dado f_{0} pertencente a F(y_{1}, ..., y_{n-1}). Mais precisamente, as frações f são dadas por f_{0}(k_{1}, k_{3}, ..., k_{n}), onde
# k_{i} = x_{i} - \sum_{j=0}^{i-1}\zeta_{ij}(x_{2}/r_{21} x_{1})^{j}/j!
# e
# \zeta_{ij} = \sum_{l_{1}=1}^{i-j}\sum_{l_{2}=l_{1}+1}^{i-j+1}r_{l_{2}l_{1}}\sum_{l_{3}=l_{2}+1}^{i-j+2}r_{l_{3}l_{2}}...\sum_{l_{j}=l_{j-1} + 1}^{i-1}r_{l_{j},l_{j}-1} r_{il_{j}} k_{l_{1}}

# O programa tem como entrada uma matriz M
# 0       0       0       0      ...  0         0
# r_{21}  0       0       0      ...  0         0
# r_{31}  r_{32}  0       0      ...  0         0
# .       .       .       .           .         .
# .       .       .       .           .         .
# .       .       .       .           .         .
# r_{n1}  r_{n2}  r_{n3}  r_{n4} ...  r_{n,n-1} 0

# e uma fração f_{0}. A saída é uma fração f.

# Função auxiliar para calcular \zeta_{ni}

def partZeta(M, n, i, j, a):
    s = 0
    if n - i + a == n - 1:
        for l in range(j + 1, n - i + a):
            s = s + M[l, j]*M[n - 1, l]
        return s
    else:
        for l in range(j + 1, n - i + a):
            s = s + M[l, j]*partZeta(M, n, i, l, a + 1)
        return s

# Função para calcular \zeta_{ni}
def zeta(k, M, n, i):
    s = 0
    if i == 1:
        for j in range(n - i):
            s = s + M[n - 1,j]*k[j]
    else:
        for j in range(n - i):
            s = s + partZeta(M, n, i, j, 1)*k[j]
    return s
    
# Função principal para calcular a função f
def invar(M, f0):
    n = M.nrows() # dimensão da F-álgebra de Lie L
    R = PolynomialRing(QQ, n, "x") 
    x = R.gens()
    FracR = FractionField(R) 
    k = zero_vector(FracR, n) # vetor de zeros para guardar os valores de k_{i}
    k[0] = x[0]
    t = x[1]/(M[1,0]*x[0])
    Z = matrix(FracR, n) # matriz de zeros para guardar os valores de zeta
#   \zeta_{1,0}   0             0             0
#   \zeta_{2,0}   \zeta_{2,1}                 0
#   \zeta_{3,0}   \zeta_{3,1}   \zeta_{3,2}   0
#   . 
#   . 
#   . 
#   \zeta_{n,0}   \zeta_{n,1}   ...      \zeta_{n,n-1}
    Z[1,1] = M[1,0]*k[0]
    for i in range(1, n):
        for l in range(1, i + 1):
            Z[i,l] = zeta(k, M, i + 1, l)   
        sZ = 0
        for l in range(i + 1):
            sZ = sZ + Z[i,l]*t^(l)/factorial(l)
        k[i] = x[i] - sZ
    q = [0]*(n - 1)
    q[0] = k[0]
    q[1] = k[2]
    for i in range(2, n - 1):
        q[i] = k[i + 1]
    return factor(f0(q))


def mult_table_gap_nilpotent_lie_algebra( F, pars ):

    l = gap.NilpotentLieAlgebra( F, pars )
    d = int(gap.Dimension( l ))
    P = PolynomialRing( F, d, x )
    b = gap.Basis( l )

    mat = matrix( P, d, d )

    for i in range( d ):
        for j in range( d ):
            v = gap.Coefficients( b, b[i+1]*b[j+1] )
            mat[i,j] = sum( v[x+1]*P.gens()[x] for x in range( d ))

    return mat
