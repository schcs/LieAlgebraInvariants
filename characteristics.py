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
    l = len( coeffs )

    # t is the variable of the curve
    t = var( 't' )

    # build the list of functions c0,...,ck where k = l-1
    c_funcs = tuple( [ function( 'c'+str(k) )(t) for k in range( l )])

    # the rhs of the equations is obtained by substituting c_funcs into the coeffs of d_opo
    rhs_eqs = [ c(c_funcs) for c in coeffs ]

    #return c_funcs, rhs_eqs
    return  [ diff( c_funcs[k], t ) == rhs_eqs[k] for k in range( l )], c_funcs


def integral_curves( d_op, P ):

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
    nr_gens = len( P.gens())
    # build the list of equations and the list of functions for d_op
    eqs, c_funcs = equations_from_differential_operator( d_op )
    
    # we need initial condition
    first_nonzero = next( k for k in range( len( eqs )) if not eqs[k].rhs().is_zero())

    if eqs[first_nonzero].rhs().coefficient( c_funcs[first_nonzero] ) == 0:
        init_k = 0
    else:
        init_k = 1
    
    print( "constant in initial condition is", init_k )

    init_cond = [0] + [ P.gens()[k] for k in range( nr_gens -1 )]
    init_cond.insert( first_nonzero+1, init_k )
    
    sols = [ x.rhs() for x in desolve_system( eqs, c_funcs, init_cond )]
    
    # substitute e^t with t. Check if this is valid!!!
    sols = [ x.subs( { e**t: t } ) for x in sols ]
    return [ P(x) for x in sols ]


def invert_curves( int_curve ):

    P = int_curve[0].parent()
    nr_vars = len( P.gens())
    t = P.gens()[nr_vars-1]

    first_with_t = next( c for c in int_curve if c != 0 and c % t == 0 )
    print( first_with_t )
    q = t #first_with_t/t
    
    ims = [ is_element_of_subalgebra_localization( int_curve, y, q ) for y in [q] + [ x for x in P.gens() ]]
    F = ims[0][0].parent().fraction_field()
    ims_new = [ ims[k][0]/ims[0][0]**ims[k][1] for k in range( 1, len( ims ))]

    return ims_new


def field_isomorphism_from_curves( int_curve ):

    P = int_curve[0].parent()
    nr_gens = len( P.gens())
    F = PolynomialRing( QQ, nr_gens-1, ['y'+str(k) for k in range( nr_gens-1 )]).fraction_field() 
    #int_curve[0].parent().fraction_field()

    t = int_curve[0].parent().gens()[nr_gens-1]
    #int_curve_t0 = [ c.subs( {t:0} ) for c in int_curve ]
    inv_c = invert_curves( int_curve )
    #return inv_c, F
    fhom = F.hom( inv_c[0:len(inv_c)-1] )
    return fhom

def field_isomorphism_from_differential_operator( d ):

    F0 = d.parent().base_ring()
    nr_gens = F0.ngens()
    F = PolynomialRing( QQ, nr_gens, [ 'y'+str(k) for k in range( nr_gens-1 )] + ['t'] )
    curves = integral_curves( d, F )

    return field_isomorphism_from_curves( curves )