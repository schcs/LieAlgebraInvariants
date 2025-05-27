"""
This module provides functions to compute rational invariants of nilpotent
Lie algebras. It includes utilities for working with polynomial rings,
derivations, and generating sets of invariants. The main functionality is
to compute algebraically independent generators for the rational invariant
field of a given nilpotent Lie algebra.
"""

from sage.algebras.lie_algebras.structure_coefficients import (
    LieAlgebraWithStructureCoefficients)
from sage.rings.derivation import RingDerivationWithoutTwist
from membership_pols import is_element_of_subalgebra
from auxfunctions import _triangular_basis_nilpotent_lie_algebra, _polynomial_ring
from dixmier import generators_of_kernel_with_dixmier_map
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.fraction_field import FractionField, FractionField_generic


lie_algebra_type = LieAlgebraWithStructureCoefficients
derivation_type = RingDerivationWithoutTwist


class RationalInvariantField(FractionField_generic, UniqueRepresentation):
    def __init__(self, lie_algebra, has_triangular_basis=False):

        # Both of these are used by _compute_symspace_gens()
        self._lie_algebra = lie_algebra
        self._has_triangular_basis = has_triangular_basis

        self._symspace_gens = self._compute_symspace_gens()

        if not self._symspace_gens:
            raise ValueError("the Lie algebra does not have a rational "
                             "invariant field")
        ngens = len(self._symspace_gens)
        R = PolynomialRing(self._lie_algebra.base_ring(), 't', ngens)
        FractionField_generic.__init__(self, R)

        SF = self._symspace_gens[0].parent().fraction_field()
        self.lift = self.hom([SF(g) for g in self._symspace_gens], check=False)
        # Make this act like a method using the generators in the symmetric space
        self.lift.register_as_coercion()

    def to_symmetric_space(self, elt):
        return self.lift(elt)

    def symmetric_space(self):
        return self._symspace_gens[0].parent()
    
    def _repr_(self):
        """
        Return a string representation of the RationalInvariantField object.
        """
        return f"Rational Invariant Field of {self._lie_algebra}"

    def __contains__(self, element):
        """
        Check if ``element`` lies in ``self``.

        INPUT:
        
        - ``element`` -- the rational function to check
        """
        P = self._polynomial_ring
        basis = self._lie_algebra.basis()
        for b in basis:
            der_coeffs = [P(self._lie_algebra.bracket(b, x)) for x in basis]
            der = P.derivation(der_coeffs)
            if der(element):
                return False
        return True

    def _compute_symspace_gens(self):
        """
        Compute the rational invariant field of a nilpotent Lie algebra.
    
        This function computes a set of algebraically independent generators for
        the rational invariant field of a given nilpotent Lie algebra. The
        computation is based on the derivations of the Lie algebra and their
        actions on a polynomial ring injected into the Lie algebra.

        OUTPUT:
        
        A list of algebraically independent generators for the rational
        invariant field of the Lie algebra.
    
        EXAMPLES::
        
            sage: l = lie_algebras.Heisenberg(QQ, 3)
            sage: rational_invariant_field(l)
            [z]
        """
        # setting up
        from sage.matrix.special import zero_matrix
        from sage.modules.free_module_element import zero_vector
        from sage.arith.misc import GCD
        from sage.misc.misc_c import prod
        lie_alg = self._lie_algebra
        l_dim = lie_alg.dimension()
        K = lie_alg.base_ring()
      
        if self._has_triangular_basis:
            bl = list(lie_alg.basis())
            from sage.matrix.special import identity_matrix
            basis_trans_matrix = identity_matrix(K, l_dim)
            basis_trans_matrix_inv = identity_matrix(K, l_dim)
        else:
            from sage.matrix.constructor import Matrix
            bl = _triangular_basis_nilpotent_lie_algebra(lie_alg)
            basis_trans_matrix = Matrix(K, l_dim, l_dim,
                                        [x.to_vector() for x in bl])
            basis_trans_matrix_inv = basis_trans_matrix.inverse()
        # basis trans_matrix contains the matrix of the identity 
        # transformation in the bases bl -> standard basis
    
        # new polynomial ring whose generators correspond to the 
        # element in bl
        P = PolynomialRing(K, l_dim, 'z')
        Pt = P
    
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
                lprod_vect = lie_alg.bracket(bl[x], bl[y]).to_vector()
                # write l_prod vect as l.c. in the triangular basis
                lprod_bl = lprod_vect*basis_trans_matrix_inv
                table[x, y] = sum(lprod_bl[k]*Pt.gens()[k] for k in range(l_dim))
        # set up initial denominator and substitution
        denom = Pt.one()
        subs = {}
    
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
            l_derivations = {x: oldPt.derivation(list(table[x]))
                             for x in range(k, l_dim)}
    
            # write the denominator as polynomial in the new generators
            # does this always work???
            denom_in_t = is_element_of_subalgebra(gens, denom, 1, Pt=Pt)[1]
    
            for y in range(len(gens)):
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
                lcm_denom = prod(list_denoms) / gcd(x.denominator() for x in coeffvec)
                if lcm_denom != 1:
                    gens[y] *= lcm_denom
                    coeffvec *= lcm_denom
                newtable[:, y] = coeffvec
            # update table
            table = newtable
    
        # return the final generating set under substitution by subs
        gens_subs = [x.subs(subs) for x in gens]
        if not self._has_triangular_basis:
            z_subs = dict(zip(P.gens(), [_lie_alg_element_to_pol(x) for x in bl]))
            gens_subs = [x.subs(z_subs) for x in gens_subs] 
    
        return gens_subs

def _lie_alg_element_to_pol(x):
    lie_alg = x.parent()
    P = _polynomial_ring(lie_alg)
    mons_x = x.monomial_coefficients()
    return sum(mons_x[x] * P(x) for x in mons_x)


def generators_of_kernel_triangular_derivation(d_op):
    """
    Compute the generators of the kernel of a triangular derivation.

    This function takes a derivation ``d_op`` and computes the generators of its
    kernel. The derivation is assumed to be triangular, meaning that it has a
    triangular action on the polynomial ring.

    INPUT:

    - ``d_op`` -- RingDerivationWithoutTwist; the derivation for which to
      compute the kernel

    OUTPUT:

    A list of generators of the kernel of the derivation.

    EXAMPLES::

        sage: P = PolynomialRing(QQ, 4, 't')
        sage: P.inject_variables()
        Defining t0, t1, t2, t3
        sage: d = P.derivation( [0,t0,t1,t2] )
        sage: generators_of_kernel_triangular_derivation(d)
        [t0, (-1/2*t1^2 + t0*t2)/t0, (1/3*t1^3 - t0*t1*t2 + t0^2*t3)/t0^2]
    """
    # the number of generators of P
    nr_gens = len(d_op.domain().gens())
    P = d_op.domain()
    if d_op == 0:
        return P.gens()

    Pt = PolynomialRing(P, 't')
    Pgens = P.gens()
    t = Pt.gen(0)

    coeffs = [d_op(x) for x in P.gens()]
    k = next((x for x in range(nr_gens) if coeffs[x]), None)
    zk = Pgens[k]
    fk = coeffs[k]
    ck = fk*t + Pgens[k]
    gens = list(Pgens[:k])
    substitution = dict(zip(gens, gens))
    substitution[Pgens[k]] = ck

    for i in range(k+1, nr_gens):
        Fi = Pt(coeffs[i].subs(substitution)).integral(t)
        substitution[Pgens[i]] = Fi + Pgens[i]
        gens.append(Fi.subs({t: -zk/fk}) + Pgens[i])

    return gens

def reduce_gen_set(gen_set):
    """
    Reduce a generating set of invariants.

    This function takes a set of generators for the invariant field and reduces
    it by removing redundant elements. Specifically, it divides each generator
    by any other generator that divides it, ensuring that the resulting set is
    minimal.

    Parameters:
    gen_set (list): A list of generators for the invariant field.

    Returns:
    list: A reduced list of generators for the invariant field.

    Example:
    """
    nr_gens = len(gen_set)
    for i in range(nr_gens):
        for j in range(i+1, nr_gens):
            while gen_set[j] % gen_set[i] == 0:
                gen_set[j] //= gen_set[i]

    return gen_set

# Monkey patching
setattr(derivation_type, "generators_of_kernel",
        generators_of_kernel_triangular_derivation)
setattr(lie_algebra_type, 'rational_invariant_field', RationalInvariantField)
