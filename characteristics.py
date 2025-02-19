from sage.all import SR, function, diff, PolynomialRing, QQ, desolve_system, var, desolve_system 
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
    P = d_op.domain()
    coeffs = [d_op(x) for x in P.gens()]
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
    for i in range(k+1,nr_gens):
        Fi = Pt(coeffs[i].subs(substitution)).integral(t)
        substitution[P.gens()[i]] = Fi+P.gens()[i]
        gens.append(Fi.subs({t:-zk/fk})+P.gens()[i])

    return gens

def generators_of_kernel(d_op):
    r"""

    """
    if d_op == 0:
        return d_op.domain().gens()

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
    
    FF = lie_alg.base_ring()
    lie_alg.polynomialRing = PolynomialRing(FF, lie_alg.dimension(),
                                            list(lie_alg.basis()),
                                            order="invlex")
    lie_alg.fractionField = lie_alg.polynomialRing.fraction_field()


def rational_invariant_field(self):

    # get the dimension of the Lie algebra and construct its
    # polynomial algebra
    d = self.dimension()
    if not hasattr(self, "polynomialRing"):
        inject_pol_ring(self)
    P = self.polynomialRing
    FF = self.base_ring()

    # lists of generators for P and list of basis for lie_alg
    # gens of P in reverse order 
    gens = P.gens()[::-1]
    bas = [x for x in self.basis()]
    
    # the current denominator
    denom = 1

    # do the following computation for each basis element of lie_alg
    for i in range(d):
        # di is the current differential operator
        di = differential_operator(self, bas[i])

        # compute the current differential operator
        # in terms of the current generating set.
        # initialize the coffs with zero-vector
    
        coeffs = [0]*len(gens)
        # the coefficients of di will be expressed in terms of
        # t1,...,t_{len(gens)}, t_{len(gens)+1} (the last indeterminate  
        # corresponds to denom)
        Pt = PolynomialRing(FF, len(gens),
                            ['t'+str(i) for i in range(len(gens))])
        
        # set up the pol ring for cs
        Ptt = PolynomialRing(FF, len(gens)+1, names='t')
        
        # new denom will be calculated
        new_denom = 1
        for k in range(len(gens)):
            # apply di to gens[k]
            d_gen = di(gens[k])

            # write d_gen as polynomial in the earlier generators
            if d_gen == 0:
                coeffs[k] = (0,1)
            else:
                v, cs = _is_element_of_subalgebra(gens, d_gen, denom, Pt=Ptt)
                assert v
                # substitute into the second components of cs
                coeffs[k] = cs
                if new_denom == 1: 
                    new_denom = d_gen
        
        coeffs = [Ptt(x[0]) for x in coeffs]
        
        # construct the differential operator in terms of the current gens
        # the indeterminates are gonna be t1,...,tk where k is #gens

        
        dt = differential_operator_from_coeffs(Pt, [Pt(x) for x in coeffs])
        if dt == 0:
            continue
                
        # compute generators for the kernel of dt
        dt_kernel_gens = dt.generators_of_kernel()
        dt_kernel_gens_enum = [Pt(x.numerator()) for x in dt_kernel_gens]
        
        # dictionary for substitution
        substitution = {Pt.gens()[i]: gens[i] for i in range(len(gens))}
        gens = [dt_k.subs(substitution) for dt_k in dt_kernel_gens_enum]
        print(gens)
        substitution = {Pt.gens()[i]: gens[i] for i in range(len(gens))}
        
        # we need to fix up the generating set 
        Ptt = PolynomialRing(FF, len(gens)+1, names='t')
        for ng in range(len(gens)):
            ng_factor = P(1)
            for j in range(i+1,d):
                dj = differential_operator(l, bas[j])
                dj_g = dj(gens[ng])
                v, cs = _is_element_of_subalgebra(gens[:ng], dj_g, denom, Pt=Ptt, denom_var = Ptt.gens()[len(gens)])
                assert v 
                if cs[1] != 1: 
                    # breakpoint()
                    cs_denom_in_x = Pt(cs[1]).subs({Ptt.gens()[len(gens)]:denom})
                    ng_factor = lcm(ng_factor,cs_denom_in_x)
            if ng_factor != 1:
                # breakpoint()
                # print("multiply with", ng_factor)
                gens[ng] *= ng_factor
                substitution = {Pt.gens()[i]: gens[i] for i in range(len(gens))}
        
        substitution = {Pt.gens()[i]: gens[i] for i in range(len(gens))}
        denom_in_gens = _is_element_of_subalgebra( gens, P(denom), new_denom, Pt=Ptt)
        if denom_in_gens[1][1] != 1:
            breakpoint()
        assert denom_in_gens[0]
        denom_in_gens = Pt(denom_in_gens[1][0]).subs(substitution)
        denom = new_denom*denom_in_gens 
        # print( "denom now is ", denom )
        
        
        # print( dt, denom )
        
    return gens, denom


def rational_invariant_field2(l):

    # setting up
    bl = list(l.basis())
    l_dim = l.dimension()
    inject_pol_ring(l)
    P = l.polynomialRing
    Pt = P #PolynomialRing(QQ, l_dim, 't', order='invlex')
    K = l.base_ring()

    # step 0
    gens = P.gens()
    # aux pol ring
    table = zero_matrix(Pt, l_dim, l_dim)
    for x in range(l_dim):
        for y in range(l_dim):
            lprod = l.bracket(bl[x], bl[y])
            prod_dict = lprod.monomial_coefficients()
            table[x, y] = sum( prod_dict[k]*P(k) for k in prod_dict)
    #subs = dict(zip(gens,Pt.gens()))
    #subs = dict(zip(Pt.gens(),gens))
    denom = Pt(1)
    subs = {}

    for k in range(l.dimension()):
        print(k)
        subs = dict(zip(Pt.gens(), [x.subs(subs) for x in gens]))
        d = Pt.derivation(list(table[k]))
        denom = Pt(next((x for x in table[k] if x), Pt(1)))
        gens = [Pt(x.numerator()) for x in generators_of_kernel_triangular_derivation(d)]        
        oldPt, Pt = Pt, PolynomialRing(K, len(gens), 't', order='invlex')
        Ft = Pt.fraction_field()
        newtable = zero_matrix(Pt, l_dim, len(gens))
        diff_ops = {x:oldPt.derivation(list(table[x])) for x in range(k,l_dim)}
        denom_in_t = _is_element_of_subalgebra(gens, denom, 1, Pt=Pt)[1]
        for y in range(len(gens)):
            coeffvec = zero_vector(Ft,l_dim)
            for x in range(k, l_dim):
                dxy_in_t = _is_element_of_subalgebra(gens, diff_ops[x](gens[y]), denom, 
                                                     denom_in_t = denom_in_t, Pt=Pt)
                coeffvec[x] = Ft(dxy_in_t[1])
            list_denoms = [x.denominator() for x in coeffvec]
            lcm_denom = prod(list_denoms)/gcd(list_denoms)
            if lcm_denom != 1:
                gens[y] *= lcm_denom
                coeffvec *= lcm_denom
            newtable[:, y] = coeffvec
        table = newtable

    return [x.subs(subs) for x in gens]




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
