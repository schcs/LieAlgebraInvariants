def _variable_into_uea(uea,v):
     return uea.gens()[str(v)]

def _sign_perm(perm):
     return Permutation([x+1 for x in perm]).signature()

def _symmetrization_mon(uea, mon):
     mon_vars = sum((x[1]*[x[0]] for x in mon.factor()),[])
     k = len(mon_vars)
     F = uea.base_ring()
     perms = Permutations(range(k))
     sym_mon = (F(1/factorial(k)))*sum(prod(_variable_into_uea(uea,mon_vars[p[i]]) for i in range(k)) for p in perms)
     return sym_mon

def symmetrization(uea, pol):
    F = uea.base_ring()
    mons = pol.monomials()
    return sum(F(pol.coefficient(m))*_symmetrization_mon(uea, m) for m in mons)
