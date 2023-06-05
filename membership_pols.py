#-------------
def alg_dependence( gens ):
    r''' 
    INPUT: 
        
        - list of polynomials in multivariate polynomial ring
    
    OUTPUT: 
    
        - Groebner basis for the ideal of polynomials that give zero if the elements of gens are subtituted.
    
    EXAMPLES:
    
        sage: gens[0]*gens[1]+gens[1]*gens[2]-gens[0]

        x0^2*x1 - x0^2*x2 - x0^2 + x1*x2*x3 - x2^2*x3

        sage: gens = [x0^2,x1-x2,x2*x3,x0^2*x1 - x0^2*x2 - x0^2 + x1*x2*x3 - x2^2*x3]

        sage: dep = alg_dependence( gens )

        sage: dep

        [t0*t1 - t0 + t1*t2 - t3]

        sage: dep[0].subs( t0=gens[0],t1=gens[1],t2=gens[2], t3 = gens[3] )

        0
        
    '''
    
    # get the parent pol ring
    P = gens[0].parent()
    P_gens = P.gens()
    nr_gens = len(gens)
    n = len(P.gens())
    F = P.base_ring()
    R = PolynomialRing( F, n+nr_gens, 'z' )
    rgens = R.gens()

    rxgens, rzgens = rgens[0:n], rgens[n:n+nr_gens]
    idgens = [ gens[i].subs( { P_gens[j]: rxgens[j] for j in range(n)}) - 
                        rzgens[i] for i in range( nr_gens )]

    I = R.ideal( idgens )
    gr = I.groebner_basis()
    
    n_zeros = tuple( 0 for _ in range( n ))
    #return gr, n, n_zeros
    cent_gens = [ g for g in gr if g.degrees()[0:n] == n_zeros ]
    #return cent_gens

    R1 = PolynomialRing( F, nr_gens, 't', order = "lex" )
    values = { R.gens()[i]: 0 for i in range( n )}
    values.update( { R.gens()[i+n]: R1.gens()[i] for i in range( nr_gens )})
    new_gens = [ x.subs( values ) for x in cent_gens ]
    return R1.ideal( new_gens ).groebner_basis()
#-------------

#-------------
def is_element_of_subalgebra( gens, p ):
    '''
        INPUT:
        
            - a list of polynomials `gens`
            - a polynomial `p`
            
        OUTPUT:
            
            - true/false depending on whether p lies in the algebra generated by gens
            - in case of true, returns a polynomial r such that r( gens ) == p 

        Finds if the polynomial `p` lies in the algebra generated by gens.
        
        EXAMPLES:
        
            sage: gens = [x0^2,x1-x2,x2*x3]

            sage: p = x0^2*x1 - x0^2*x2 - x0^2 + x1*x2*x3 - x2^2*x3

            sage: v, dep = is_element_of_subalgebra( gens, p )

            sage: v

            True

            sage: dep

            [t0*t1 - t0 + t1*t2]

            sage: dep[0].subs( t0=gens[0],t1=gens[1],t2=gens[2] ) == p

            True
    '''
    
    deps = alg_dependence( gens+[p] )
    d = len( gens )
    
    if len( deps ) == 0:
        return false, _
    
    R = deps[0].parent()
    deps = [ x for x in deps if R.gens()[d] in x.monomials() ]
    coeffs = [ x.coefficient( R.gens()[d] ) for x in deps ]
    
    return len( deps ) > 0, [ deps[i] - coeffs[i]*R.gens()[d] for i in range( len( deps )) ] 
#-------------