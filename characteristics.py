from sage.all import PolynomialRing, zero_matrix, zero_vector, prod, gcd
from sage.algebras.lie_algebras.structure_coefficients import (
    LieAlgebraWithStructureCoefficients
)
from sage.rings.derivation import RingDerivationWithoutTwist
from membership_pols import is_element_of_subalgebra

lie_algebra_type = LieAlgebraWithStructureCoefficients
derivation_type = RingDerivationWithoutTwist


def differential_operator_from_coeffs(P, coeffs):
    r"""
    """

    D = P.derivation_module()
    op = D.zero()

    for k, _ in enumerate(coeffs):
        op += coeffs[k]*D.gens()[k]

    return op


def generators_of_kernel_triangular_derivation(d_op):

    # the number of generators of P
    nr_gens = len(d_op.domain().gens())
    P = d_op.domain()
    if d_op == 0:
        return P.gens()

    Pt = PolynomialRing(P, 't')
    t = Pt.gens()[0]

    coeffs = [d_op(x) for x in P.gens()]
    k = next((x for x in range(nr_gens) if coeffs[x]), None)
    zk = P.gens()[k]
    fk = coeffs[k]
    ck = fk*t+P.gens()[k]
    gens = list(P.gens()[:k])
    substitution = dict(zip(P.gens()[:k], gens))
    substitution[P.gens()[k]] = ck
    for i in range(k+1, nr_gens):
        Fi = Pt(coeffs[i].subs(substitution)).integral(t)
        substitution[P.gens()[i]] = Fi+P.gens()[i]
        gens.append(Fi.subs({t: -zk/fk})+P.gens()[i])

    return gens


def inject_pol_ring(lie_alg):

    FF = lie_alg.base_ring()
    lie_alg.polynomialRing = PolynomialRing(FF, lie_alg.dimension(),
                                            list(lie_alg.basis()),
                                            order="invlex")
    lie_alg.fractionField = lie_alg.polynomialRing.fraction_field()


def rational_invariant_field(lie_alg):

    # setting up
    bl, l_dim = list(lie_alg.basis()), lie_alg.dimension()
    inject_pol_ring(lie_alg)
    P = lie_alg.polynomialRing
    Pt, K = P, lie_alg.base_ring()

    # step 0
    # gens contains the generators of invariants after applying
    # the derivation that corresponds to the basis elements
    gens = P.gens()

    # table contains how each derivation acts on the current generating set
    # the initial generating set is the basis of l, so initually, this
    # is just the multiplication table
    table = zero_matrix(Pt, l_dim, l_dim)
    for x in range(l_dim):
        for y in range(l_dim):
            # essentially,
            # table[x, y] = P(l.bracket(bl[x], bl[y]))
            # however, this does not work over more complicated fields
            lprod = lie_alg.bracket(bl[x], bl[y])
            prod_dict = lprod.monomial_coefficients()
            table[x, y] = sum(prod_dict[k]*P(k) for k in prod_dict)
    # set up initial denominator and substitution
    denom, subs = Pt(1), {}

    # start computation here
    for k in range(lie_alg.dimension()):
        # compute the new substitution. subs contains at each step, the
        # expressions for the generators in terms of the original generators
        subs = dict(zip(Pt.gens(), [x.subs(subs) for x in gens]))

        # the current derivation is the derivation of Pt whose coefficients
        # are in the k-th line of table
        d = Pt.derivation(list(table[k]))

        # denom is the first non-zero coefficient of d
        denom = Pt(next((x for x in table[k] if x), Pt(1)))

        # compute the generators of the kernel of d and take their enumerators
        # this is valid, since the denominators are invariants
        gens = [Pt(x.numerator())
                for x in generators_of_kernel_triangular_derivation(d)]

        # construct the new ring Pt whose rank is len(gens)
        # and remember the old Pt as oldPt
        oldPt, Pt = Pt, PolynomialRing(K, len(gens), 't', order='invlex')
        Ft = Pt.fraction_field()

        # recompute the table for the action of the derivations on the new
        # generating set
        newtable = zero_matrix(Pt, l_dim, len(gens))

        # compute the derivations of l as they act on the generators
        l_derivations = {x: oldPt.derivation(
            list(table[x])) for x in range(k, l_dim)}

        # write the denominator as polynomial in the new generators
        # does this always work???
        denom_in_t = is_element_of_subalgebra(gens, denom, 1, Pt=Pt)[1]

        for y, _ in enumerate(gens):
            coeffvec = zero_vector(Ft, l_dim)
            for x in range(k, l_dim):
                # compute the image of gens[y] under l_derivation[x] as
                # expressions in the generators
                dxy_in_t = is_element_of_subalgebra(gens,
                                                    l_derivations[x](gens[y]),
                                                    denom,
                                                    denom_in_t=denom_in_t,
                                                    Pt=Pt)
                # put it into coeffvec
                coeffvec[x] = Ft(dxy_in_t[1])
            # compute the denominators and multiply everything with the lcm
            # of the denominators
            list_denoms = [x.denominator() for x in coeffvec]
            lcm_denom = prod(list_denoms)/gcd(list_denoms)
            if lcm_denom != 1:
                gens[y] *= lcm_denom
                coeffvec *= lcm_denom
            newtable[:, y] = coeffvec
        # update table
        table = newtable

    # return the final generating set under substitution by subs
    return [x.subs(subs) for x in gens]


def reduce_gen_set(gen_set):
    nr_gens = len(gen_set)
    for i in range(nr_gens):
        for j in range(i+1, nr_gens):
            while gen_set[j] % gen_set[i] == 0:
                gen_set[j] //= gen_set[i]

    return gen_set


setattr(derivation_type, "generators_of_kernel",
        generators_of_kernel_triangular_derivation)
setattr(lie_algebra_type, 'rational_invariant_field', rational_invariant_field)
