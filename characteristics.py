from sage.all import SR, function, diff, PolynomialRing, QQ, desolve_system, var
from sage.algebras.lie_algebras.structure_coefficients import LieAlgebraWithStructureCoefficients
from sage.rings.derivation import RingDerivationWithoutTwist
from invariants_nilp_lie_alg import differential_operator, differential_operator_from_coeffs
from membership_pols import is_element_of_subalgebra

lie_algebra_type = LieAlgebraWithStructureCoefficients
derivation_type = RingDerivationWithoutTwist

"""
using the sage solver:

sage: l = solvable_lie_algebra( QQ, [4,2,1] )
sage: d = differential_operator( l, l.3 );
sage: system4 = equations_from_differential_operator( d ); system4
([diff(c0(t), t) == c0(t),
  diff(c1(t), t) == c1(t),
  diff(c2(t), t) == c2(t),
  diff(c3(t), t) == 0],
 (c0(t), c1(t), c2(t), c3(t)))
sage: desolve_system(*system4, [0,1,2,3,1])
[c0(t) == e^t, c1(t) == 2*e^t, c2(t) == 3*e^t, c3(t) == 1]
"""


def equations_from_differential_operator(d_op):
    '''
        INPUT:

            - a differential operator d_op

        OUTPUT:

            - The differential equations that determine the integral curves
              that correspond to d_op.

        COMMENT:

            Given a differential operator d_op, of the form

            p_0(x)*d/dx0 + ... + p_k(x)*d/dxk,

            a curve is said to be an integral curve if it is tangent to the
            vector field (p_0....,p(k)). These integral curves are determined
            by the system of differential equations

            d/dt(c_i(t)) = p_i(c_1,...,c_k) for i=1,...,k.

            Thus function returns these equations and also the list of
            functions c_i that appear in the equations.
/2*t^2*x1
            The equations are returned as symbolic expressions.

        EXAMPLES::

            sage: l = solvable_lie_algebra( QQ, [4,2,1] )
            sage: d = differential_operator( l, l.3 )
            sage: equations_from_differential_operator( d )
            ([diff(c0(t), t) == c0(t),
              diff(c1(t), t) == c1(t),
              diff(c2(t), t) == c2(t),
              diff(c3(t), t) == 0],
             (c0(t), c1(t), c2(t), c3(t)))

    '''
    # find the coefficients of d_op
    coeffs = d_op.list()
    nr_coeffs = len(coeffs)

    # t is the variable of the curve
    t = SR.var('t')

    # build the list of functions c0,...,ck where k = l-1
    c_funcs = tuple(function(f'c{k}')(t) for k in range(nr_coeffs))

    # the rhs of the equations is obtained by substituting
    # c_funcs into the coeffs of d_opo
    rhs_eqs = [c(c_funcs) for c in coeffs]

    # return c_funcs, rhs_eqs
    return [diff(c_funcs[k], t) == rhs_eqs[k]
            for k in range(nr_coeffs)], c_funcs


def integral_curves(d_op):
    '''
        INPUT:

            - a differential operator d_op
            - a polynomial ring P

        OUTPUT:

            -

        COMMENT:



        EXAMPLES:

    '''
    # the number of generators of P
    nr_gens = len(d_op.domain().gens())
    P0 = d_op.domain()

    # build the list of equations and the list of functions for d_op
    eqs, c_funcs = equations_from_differential_operator(d_op)

    # we need initial condition
    init_cond = [0] + [P0.gens()[k] for k in range(nr_gens)]
    return desolve_system(eqs, c_funcs, init_cond)


def integral_curves_triangular_derivation(d_op):
    pass


def generators_of_kernel(d_op):
    r"""

    """
    curves = [c.rhs() for c in integral_curves(d_op)]
    t = var('t')
    gens_P = d_op.domain().gens()
    rank = len(gens_P)
    # k is the least number such that the coefficient of d_op for the
    # k-th variable is non-zero
    d_op_mon_coeffs = d_op.monomial_coefficients()
    d_op_coeffs = d_op.list()
    k = min(d_op_mon_coeffs.keys())
    xk, fk = gens_P[k], d_op_mon_coeffs[k]
    subs_dict = dict(zip(gens_P, curves))
    Fs = [x.subs(subs_dict).integral(t) for x in d_op_coeffs]
    kernel_gens = [Fs[i].subs(t=-xk/fk)+gens_P[i] for i in range(rank)
                   if i != k]
    return kernel_gens


def inject_pol_ring(lie_alg):
    lie_alg.polynomialRing = PolynomialRing(QQ, lie_alg.dimension(),
                                            list(lie_alg.basis())[::-1],
                                            order="lex")
    lie_alg.fractionField = lie_alg.polynomialRing.fraction_field()


def rational_invariant_field(self):

    # get the dimension of the Lie algebra and construct its
    # polynomial algebra
    d = self.dimension()
    if True or not hasattr(self, "polynomialRing"):
        inject_pol_ring(self)
    P = self.polynomialRing

    # lists of generators for P and list of basis for lie_alg
    gens = P.gens()[::-1]
    bas = [x for x in self.basis()]
    denoms = []

    # do the following computation for each basis element of lie_alg
    for i in range(d):
        # di is the current differential operator
        di = differential_operator(self, bas[i])

        # compute the current differential operator
        # in terms of the current generating set.
        # initialize the coffs with zero-vector
        coeffs = [0]*len(gens)
        Pt = PolynomialRing(QQ, len(gens), ['t'+str(i)
                            for i in range(len(gens))])
        for k in range(len(gens)):
            # apply di to gens[k]
            d_gen = di(gens[k])

            # write d_gen as polynomial in the earlier generators
            if d_gen == 0:
                coeffs[k] = 0
            else:
                v, cs = _is_element_of_subalgebra(gens[0:k], denoms, d_gen)
                assert v
                coeffs[k] = Pt(cs[0].numerator())

        # construct the differential operator in terms of the current gens
        # the indeterminates are gonna be t1,...,tk where k is #gens
        dt = differential_operator_from_coeffs(Pt, coeffs)

        if dt == 0:
            continue

        # compute generators for the kernel of dt
        dt_kernel_gens = dt.generators_of_kernel()
        dt_kernel_gens_enum = [Pt(x.numerator()) for x in dt_kernel_gens]
        dt_kernel_gens_deno = [Pt(x.denominator()) for x in dt_kernel_gens]
        
        # dictionary for substitution
        substitution = {Pt.gens()[i]: gens[i] for i in range(len(gens))}

        gens = [dt_k.subs(substitution) for dt_k in dt_kernel_gens_enum]
        for x in dt_kernel_gens_deno:
            x_subs = x.subs(substitution)
            if x_subs != 1 and x_subs not in denoms:
                denoms.append(x.subs(substitution))

        # gens = reduce_gen_set(gens)
    return gens, denoms


def reduce_gen_set(gen_set):
    nr_gens = len(gen_set)
    for i in range(nr_gens):
        for j in range(i+1, nr_gens):
            while gen_set[j] % gen_set[i] == 0:
                gen_set[j] //= gen_set[i]
    
    return gen_set


setattr(derivation_type, "integral_curves", integral_curves)
setattr(derivation_type, "generators_of_kernel", generators_of_kernel)
setattr(lie_algebra_type, 'rational_invariant_field', rational_invariant_field)
