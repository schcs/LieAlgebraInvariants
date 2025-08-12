from sage.all import ZZ, factorial

# The image of f under the Dixmier map corresponding to a derivation D,
# a local slice r.
#
# See Freudenburg "The algebraic theory of locally nilpotent derivations",
# Section 1.1.9, page 12.


def dixmier_map(D, r, f):
    """
    Computes the Dixmier map for a given differential operator D, a local slice r, 
    and argument f. A local slice is a function r such that D(r) != 0, but D^2(r) = 0. 
    The Dixmier map is represented as:

        φ(f) = Σ ( (-1)^i / i! ) * ( D^i(f) * r^i ) / ( D(r)^i )
              i>=0

    Parameters:
    D (function): A differential operator that acts on functions.
    r (function): A local slice in the computation of the Dixmier map.
    f (function): The argument function to which the Dixmier map is applied.
    
    Returns:
    result (expression): The image of f under the Dixmier map.
    
    Notes:
    - The computation iteratively applies the differential operator D to f 
      until the derivative becomes zero.
    - The factorial and alternating sign are used in the series expansion.
    """
    i, result = 0, 0
    Df = f
    while True:
        result += ZZ((-1)**i)/factorial(i)*Df*r**i/D(r)**i
        i += 1
        Df = D(Df)
        if Df == 0:
            break

    return result


def generators_of_kernel_with_dixmier_map(D):
    """
    Compute the generators of the kernel of a given
    triangular derivation D using the Dixmier map.

    Parameters:
    D (Derivation): The derivation for which to compute the kernel generators.

    Returns:
    list: A list of generators for the kernel of the derivation D.
    """
    # Get the field and its base field

    F = D.base_ring()
    gens_F = F.gens()
    if D == 0:
        return gens_F

    k = next((k for k, xk in enumerate(gens_F) if D(xk) != 0), 0)
    r = gens_F[k]
    return [dixmier_map(D, r, gens_F[j])
            for j, _ in enumerate(gens_F) if j != k]
