from sage.all import PolynomialRing 


def _triangular_basis_nilpotent_lie_algebra(lie_alg):
    """
    Compute a triangular basis for a nilpotent Lie algebra.
    This function calculates a basis for a nilpotent Lie algebra using its 
    upper central series. The basis is constructed by lifting elements 
    from the quotient spaces of consecutive terms in the upper central series.
    
    Args:
        l: A Lie algebra object that provides methods for computing the 
           ideal, upper central series, and quotient spaces.
    
    Returns:
        list: A list of elements forming a triangular basis for the nilpotent 
              Lie algebra.
    
    Examples:
        sage: L = nilpotent_lie_algebra(QQ, [6, 14])
        sage: basis = _triangular_basis_nilpotent_lie_algebra(L)
        sage: basis
        [x5, x4, x3, x2, x0, x1]
    """
    
    ucs = [lie_alg.ideal()] + lie_alg.upper_central_series()
    nr_quots = len(ucs)-1
    bas = []
    for i in range(nr_quots):
        quot = ucs[i+1].quotient(ucs[i].gens())
        bas += [quot.lift(x) for x in quot.basis()]
        
    return bas


def _inject_pol_ring(lie_alg):
    """
    Inject a polynomial ring into the Lie algebra.

    This function takes a Lie algebra and injects a polynomial ring over its
    base field with variables corresponding to the basis elements of the Lie
    algebra. It also sets up the fraction field of the polynomial ring.

    Parameters:
    lie_alg (LieAlgebraWithStructureCoefficients): The Lie algebra into which
    the polynomial ring is to be injected.

    Returns:
    True

    Example:
    sage: l = lie_algebras.Heisenberg(QQ,3)
    sage: inject_pol_ring(l)
    sage: l.polynomialRing
    Multivariate Polynomial Ring in p1, p2, p3, q1, q2, q3, z over Rational
    Field
    """
    FF = lie_alg.base_ring()
    lie_alg.polynomialRing = PolynomialRing(FF, lie_alg.dimension(),
                                            list(lie_alg.basis()),
                                            order="invlex")
    lie_alg.fractionField = lie_alg.polynomialRing.fraction_field()
    return True
