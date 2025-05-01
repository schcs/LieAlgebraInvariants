from sage.all import ZZ

def dixmier_map(D, r, f):

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
    k = next((k for k, xk in enumerate(F.gens()) if D(xk) != 0), 0)
    r = gens_F[k]
    return [dixmier_map(D,r,gens_F[j]) for j, _ in enumerate(gens_F) if j != k ]
