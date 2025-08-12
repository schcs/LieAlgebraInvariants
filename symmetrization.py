from sage.all import Permutation, Permutations, factorial


def _variable_into_uea(uea,v):
    return uea.gens()[str(v)]


def _sign_perm(perm):
    return Permutation([x+1 for x in perm]).signature()


def _symmetrization_mon(uea, mon):
    mon_vars = sum((x[1]*[x[0]] for x in mon.factor()),[])
    k = len(mon_vars)
    F = uea.base_ring()
    perms = Permutations(range(k))
    gens = uea.gens()
    return F(1/factorial(k)) * uea.sum(uea.prod(gens[str(mon_vars[val])] for val in p) for p in perms)


def symmetrization(uea, pol):
    F = uea.base_ring()
    mons = pol.monomials()
    return sum(F(pol.monomial_coefficient(m))*_symmetrization_mon(uea, m) for m in mons)
