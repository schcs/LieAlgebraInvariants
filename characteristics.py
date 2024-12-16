from sage.all import SR, function, diff, PolynomialRing, QQ, desolve_system, var
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
    v = d_op.list()
    first_non_zero = next(i for i in range(len(v)) if v[i] != 0)
    P0 = d_op.domain()
    namesP = [str(P0.gens()[x]) for x in range(nr_gens) if
              x != first_non_zero] + ['t']

    P = PolynomialRing(QQ, nr_gens, namesP)
    # build the list of equations and the list of functions for d_op
    eqs, c_funcs = equations_from_differential_operator(d_op)

    # we need initial condition
    init_k = 0

    init_cond = [0] + [P0.gens()[k] for k in range(nr_gens)]
    print( init_cond)
    #init_cond.insert(first_non_zero+1, init_k)
    return desolve_system(eqs, c_funcs, init_cond)


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
    kernel_gens = [Fs[i].subs(t=-xk/fk)+gens_P[i] for i in range(rank) if i != k]
    return kernel_gens


def rational_invariant_field(lie_alg):

    d = lie_alg.dimension()
    P = PolynomialRing(QQ, d, [str(x) for x in lie_alg.basis()])

    gens = [x for x in P.gens()]
    bas = [x for x in lie_alg.basis()]

    for i in range(d):
        d = differential_operator(lie_alg, bas[i])
        coeffs = [0]*len(gens)
        for k in range(len(gens)):
            d_gen = d(gens[k])
            if d_gen == 0:
                coeffs[k] = 0
            else:
                v, cs = is_element_of_subalgebra(gens[0:k],d_gen)
                assert v
                coeffs[k] = cs[0]
        

        Pt = PolynomialRing( QQ, len(gens), ['t'+str(i) for i in range(len(gens))])
        dt = differential_operator_from_coeffs(Pt, coeffs)

        if dt == 0:
            continue

        h  = field_isomorphism_from_differential_operator(dt)

        substitution = { h.codomain().gens()[i]: gens[i] for i in range( len( gens )) }

        gens = [ h(x).subs( substitution ) for x in h.domain().gens() ]
        gens = [ x.numerator() for x in gens ]

    return gens


def rational_invariant_field_2( l ):

    d = l.dimension()
    P = PolynomialRing( QQ, d, [ str( x ) for x in l.basis() ])

    gens = [ x for x in P.gens()]
    bas = [ x for x in l.basis() ]

    m = structure_constants( l, l.basis()).rref()

    for i in range( m.nrow ):

        r = m[i]*lcm( x.denominator() for x in m[i] )
        r = ( P(x) for x in r )

        d = differential_operator_with_coeffs( P, r )
        print( i, d )
        coeffs = [ is_element_of_subalgebra( gens[0:k], d(gens[k]) )[1][0] for k in range( len( gens ))]

        Pt = PolynomialRing( QQ, len( gens ), ['t'+str(i) for i in range( len( gens ))])
        dt = differential_operator_from_coeffs( Pt, coeffs )

        if dt == 0:
            continue

        h  = field_isomorphism_from_differential_operator( dt )

        substitution = { h.codomain().gens()[i]: gens[i] for i in range( len( gens )) }

        gens = [ h(x).subs( substitution ) for x in h.domain().gens() ]
        gens = [ x.numerator() for x in gens ]

    return gens
    p = l.polynomialRing

    for i in range( m.nrows()):
        r = m[i]*lcm( x.denominator() for x in m[i] )
        r = ( P(x) for x in r )


