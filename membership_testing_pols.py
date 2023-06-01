# mudança
# 
# mudança 2
def alg_dependence( gens ):

    # mudança do Csaba 2
    
    P = gens[0].parent()
    # mudança do Igor
    P_gens = P.gens()
    nr_gens = len(gens)
    n = len(P.gens())
    F = P.base_ring()
    R = PolynomialRing( F, n+nr_gens, 'z' )
    rgens = R.gens()
    print( nr_gens, n )
    rxgens, rzgens = rgens[0:n], rgens[n:n+nr_gens]
    idgens = [ gens[i].subs( { P_gens[j]: rxgens[j] for j in range(n)}) - 
                        rzgens[i] for i in range( nr_gens )]

    I = R.ideal( idgens )
    gr = I.groebner_basis()
    
    n_zeros = tuple( 0 for _ in range( n ))
    #return gr, n, n_zeros
    cent_gens = [ g for g in gr if g.degrees()[0:n] == n_zeros ]
    #return cent_gens

    R1 = PolynomialRing( F, nr_gens, 't' )
    values = { R.gens()[i]: 0 for i in range( n )}
    values.update( { R.gens()[i+n]: R1.gens()[i] for i in range( nr_gens )})
    new_gens = [ x.subs( values ) for x in cent_gens ]
    return R1.ideal( new_gens )




