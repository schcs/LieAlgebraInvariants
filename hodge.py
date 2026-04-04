from itertools import combinations
from auxfunctions import _polynomial_ring


def omega_form(l):
    n = l.dimension()
    P = _polynomial_ring(l)
    M = EuclideanSpace(n)
    f = M.diff_form(2)
    basis = list(l.basis().keys())
    struct_coeffs = l.structure_coefficients()
    for i in range(n):
        for j in range(i+1,n): 
            if (basis[i],basis[j]) in struct_coeffs.keys():
                pr = struct_coeffs[(basis[i], basis[j])]
                f[i+1,j+1] = sum([ pr.coefficient(z)*P(z) for z in basis])
    return f, M
    

def factors_of_diff_form(M,f):
    n = M.dim()
    index_sets = [tuple(c) for c in combinations(range(1, n + 1), n - 1)]
    P = l.polynomialRing
    F = P.fraction_field()

    mat = zero_matrix(F, n, n)
    for i in range(n):
        ei = zero_vector(n); ei[i] = 1
        di = M.one_form(ei)
        fw = f.wedge(di)
        mat[i,:] = [[fw[ind] for ind in index_sets]]

    ker = mat.left_kernel().basis()
    ker = Matrix(P, [LCM(x.denominator() for x in z)*z for z in ker])    
    return ker



def is_gradient(vec):
    
    P = vec[0].parent()
    DP = P.derivation_module()
    dim = len(vec)

    for i in range(dim):
        for j in range(i+1,dim):
            if DP.gens()[i](vec[j]) != DP.gens()[j](vec[i]):
                return False, i, j

    return True, -1, -1


def polynomial_from_gradient(grad, gens=None, constant=None, check=True):
    """Reconstruct a polynomial from its gradient.

    Parameters
    ----------
    grad : list|tuple
        A list/tuple ``[∂f/∂x_0, …, ∂f/∂x_{n-1}]`` of partial derivatives in a
        Sage polynomial ring ``P``.
    gens : list|tuple|None
        The variables to integrate with respect to, in the same order as
        ``grad``. Defaults to ``P.gens()[:len(grad)]``.
    constant : element|None
        Optional additive constant to return ``f + constant``. Defaults to 0 in
        the same parent ring.
    check : bool
        If True, verify that the reconstructed polynomial differentiates back
        to ``grad``.

    Returns
    -------
    element
        A polynomial ``f`` whose partial derivatives match ``grad`` (up to the
        chosen additive constant).

    Notes
    -----
    This assumes ``grad`` is a conservative vector field, i.e. mixed partials
    agree. If not, a ``ValueError`` is raised when ``check`` is True.
    """


    grad = list(grad)
    P = grad[0].parent()

    if any(g.parent() is not P for g in grad):
        raise ValueError("all entries of grad must lie in the same parent ring")

    if gens is None:
        gens = list(P.gens()[: len(grad)])
    else:
        gens = list(gens)

    if len(gens) != len(grad):
        raise ValueError("gens and grad must have the same length")

    if constant is None:
        constant = P(0)
    else:
        constant = P(constant)

    # Optional pre-check: ensure this is actually a gradient.
    if check:
        res = is_gradient(grad)[0]
        if res is not True:
            # res is expected to be (False, i, j) from is_gradient
            raise ValueError(f"grad is not integrable (mixed partials mismatch): {res}")

    DP = P.derivation_module()
    ders = DP.gens()

    # Reconstruct by successive integration; the residual recovers the
    # 'integration constant' as a function of remaining variables.
    f = grad[0].integral(gens[0])
    for i in range(1, len(grad)):
        residual = grad[i] - ders[i](f)
        f += residual.integral(gens[i])

    f += constant

    if check:
        for i in range(len(grad)):
            if ders[i](f) != grad[i]:
                raise ValueError(f"reconstruction failed: d/d{gens[i]}(f) != grad[{i}]")

    return f
 

def make_it_gradient(vec, gradients):

    # pivots contains the indices of the last nonzero entry in each vector in gradients
    pivots = [next((k for k in range(len(g) - 1, -1, -1) if g[k] != 0), -1) for g in gradients]
    
    P = vec[0].parent()
    DP = P.derivation_module()
    ders = DP.gens()
    dim = len(vec)
    Pgens = P.gens()

    for k in range(dim-1, -1, -1): # backwards from dim-1 to 0

        # compute the k-th partial derivatives of vec[0], ..., vec[k-1]
        ders_vj_dk = [ders[k](vec[j]) for j in range(dim)]
        
        try:
            piv_k = pivots.index(k)
        except ValueError:
            # there is no vector with this pivot
            continue

        grad = gradients[piv_k]
        for j in range(dim):
            if ders_vj_dk[j] == ders[j](vec[k]):
                continue
            
            entry_integral = ders_vj_dk[j].integral(Pgens[j])
            modding_el = entry_integral*grad
            coeff_vec = ders[j](modding_el[k])/ders_vj_dk[j]
            vec = coeff_vec*vec + modding_el
            assert ders[k](vec[j]) == ders[j](vec[k])

    assert is_gradient(vec)[0]  

    return vec

