from sage.all import SR, function, diff
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

def equations_from_differential_operator( d_op ):
    '''
        INPUT:

            - a differential operator d_op

        OUTPUT:

            - The differential equations that determine the integral curves that correspond to d_op.

        COMMENT:

            Given a differential operator d_op, of the form

            p_0(x)*d/dx0 + ... + p_k(x)*d/dxk,

            a curve is said to be an integral curve if it is tangent to the vector field (p_0....,p(k)).
            These integral curves are determined by the system of differential equations

            d/dt(c_i(t)) = p_i(c_1,...,c_k) for i=1,...,k.

            Thus function returns these equations and also the list of functions c_i that appear in the
            equations.

            The equations are returned as symbolic expressions.

        EXAMPLES:

            sage: l = solvable_lie_algebra( QQ, [4,2,1] )

            sage: d = differential_operator( l, l.3 );

            sage: equations_from_differential_operator( d )

            ([diff(c0(t), t) == c0(t),
              diff(c1(t), t) == c1(t),
              diff(c2(t), t) == c2(t),
              diff(c3(t), t) == 0],
             (c0(t), c1(t), c2(t), c3(t)))

    '''
    # find the coefficients of d_op
    coeffs = d_op.list()
    l = len(coeffs)

    # t is the variable of the curve
    t = SR.var('t')

    # build the list of functions c0,...,ck where k = l-1
    c_funcs = tuple(function(f'c{k}')(t) for k in range(l))

    # the rhs of the equations is obtained by substituting
    # c_funcs into the coeffs of d_opo
    rhs_eqs = [c(c_funcs) for c in coeffs]

    # return c_funcs, rhs_eqs
    return [diff(c_funcs[k], t) == rhs_eqs[k] for k in range(l)], c_funcs


def desolve_system_alt( rhss, t, func_c, P ):
    """
    examples !

    sage:       sage: x = var('x')
    ....:       sage: y = function('y')(x)
    ....:       sage: de = diff(y,x) - y == 0
    ....:       sage: desolve(de, y, [0,1])
    e^x

    """
    nr_eqs = len( rhss )
    gensP = P.gens()

    sols = []
    seen_non_zero = False
    kP = 0
    for k in range( nr_eqs ):
        expr = rhss[k]

        if not seen_non_zero and expr != 0:
            add_term = 0
            seen_non_zero = True
        else:
            add_term = gensP[kP]
            kP += 1

        sol = expr.integral(t) + add_term
        sols.append( sol )
        for j in range( k+1, nr_eqs ):
            rhss[j] = rhss[j].subs( func_c[k] == sol )

    return [ P( x ) for x in sols ]


def integral_curves( d_op ):
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
    nr_gens = len( d_op.domain().gens())
    v = d_op.list()
    first_non_zero = next( i for i in range( len( v )) if v[i] != 0 )
    P0 = d_op.domain()
    namesP = [ 'y_'+str(P0.gens()[x]) for x in range( nr_gens ) if x != first_non_zero ] + [ 't' ]

    P = PolynomialRing( QQ, nr_gens, namesP )
    t = P.gens()[P.ngens()-1]
    # build the list of equations and the list of functions for d_op
    eqs, c_funcs = equations_from_differential_operator( d_op )

    # we need initial condition
    if eqs[first_non_zero].rhs().coefficient( c_funcs[first_non_zero] ) == 0:
        # nilpotent case
        init_k = 0
    else:
        # semisimple case
        init_k = 1

    init_cond = [0] + [ P.gens()[k] for k in range( nr_gens -1 )]
    init_cond.insert( first_non_zero+1, init_k )

    # sols = desolve_system_alt( [ x.rhs() for x in eqs ], t, c_funcs, P )

    sols = desolve_system(eqs, c_funcs, init_cond)

    # substitute e^t with t. Check if this is valid!!!
    sols = [ x.rhs().subs( { e**t: t } ) for x in sols ]
    return [ P(x) for x in sols ]


def invert_curves( int_curve, P0 ):

    P = int_curve[0].parent()
    nr_vars = len( P.gens())
    t = P.gens()[-1]

    first_with_t = next( c for c in int_curve if c != 0 and c % t == 0 )
    q = first_with_t / t

    ims = [ is_element_of_subalgebra_localization( int_curve, y, q ) for y in [q] + [ x for x in P.gens() ]]
    ims_new = [ ims[k][0]/ims[0][0]**ims[k][1] for k in range( 1, len( ims ))]
    Pims = ims_new[0].parent()
    substitution = { Pims.gens()[i]: P0.gens()[i] for i in range( nr_vars )}
    return [ x.subs( substitution ) for x in ims_new ]


def field_isomorphism_from_curves( int_curve, P, domainF ):

    nr_gens = len( P.gens())

    t = int_curve[0].parent().gens()[nr_gens-1]
    inv_c = invert_curves( int_curve, P )

    fhom = domainF.hom( inv_c[0:len(inv_c)-1] )

    return fhom

def field_isomorphism_from_differential_operator( d ):

    F0 = d.parent().base_ring()
    nr_gens = len( F0.gens())
    v = d.list()
    first_non_zero = next( i for i in range( len( v )) if v[i] != 0 )
    domainF = PolynomialRing( QQ, nr_gens - 1,
                    [str(F0.gens()[x]) for x in range( nr_gens ) if x != first_non_zero]).fraction_field()

    curves = integral_curves( d )
    return field_isomorphism_from_curves( curves, F0, domainF )


def rational_invariant_field( l ):

    d = l.dimension()
    P = PolynomialRing( QQ, d, [ str( x ) for x in l.basis() ])

    gens = [ x for x in P.gens()]
    bas = [ x for x in l.basis() ]

    for i in range( d ):

        d = differential_operator( l, bas[i] )
        coeffs = [0]*len( gens )
        for k in range( len( gens )):
            d_gen = d(gens[k])
            if d_gen == 0:
                coeffs[k] = 0
            else:
                print( "k is ", k )
                print( "i is", i, "gens are ", gens[0:k], "pol is", d_gen )
                v, cs = is_element_of_subalgebra( gens[0:k], d_gen )
                assert v
                coeffs[k] = cs[0]

        print( coeffs )

        #coeffs = [ is_element_of_subalgebra( gens[0:k], d(gens[k]) )[1][0] for k in range( len( gens ))]

        print( i, d )

        Pt = PolynomialRing( QQ, len( gens ), ['t'+str(i) for i in range( len( gens ))])
        dt = differential_operator_from_coeffs( Pt, coeffs )

        if dt == 0:
            continue

        h  = field_isomorphism_from_differential_operator( dt )

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


