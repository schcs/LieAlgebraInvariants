"""
This module provides functionality to determine if a given polynomial lies
inside a subalgebra generated by a set of polynomials localized by a given
polynomial. This is not a general purpose module, the generator polynomials
have very specific form since they come from applying the method of integral
curves to a triangular derivation.


Functions:
- is_element_of_subalgebra: Determines if a polynomial is an element of a
subalgebra generated by a list of polynomials and the inverse of an
additional polynomial. The polynomials have special form, see the paragraph
above.

Dependencies:
- sage.all: Provides the PolynomialRing class used for polynomial operations.
"""

from sage.all import PolynomialRing


def is_element_of_subalgebra(gens, p, denom=1, denom_in_t=1, Pt=False,
                             denom_var=0):
    """
    Determines if a polynomial is an element of a
    subalgebra generated by a list of polynomials and the inverse of an
    additional polynomial. This is not a general purpose function, the
    generator polynomials have very specific form since they come from
    applying the method of integral curves to a triangular derivation.

    Args:
        gens (list): A list of polynomials that generate the subalgebra.
        p (Polynomial): The polynomial to check.
        denom (int, optional): The current denominator. Defaults to 1.
        denom_in_t (int, optional): The denominator in terms of t.
        Defaults to 1.
        Pt (bool or PolynomialRing, optional): The polynomial ring.
        Defaults to False.
        denom_var (int, optional): The denominator variable. Defaults to 0.

    Returns:
        tuple: A tuple containing a boolean indicating if `p` lies in the
        subalgebra and the polynomial in terms of the generators.

    Examples:
        sage: P = PolynomialRing( QQ, 10, 't', order = 'invlex' )
        sage: P.inject_variables()
        Defining t0, t1, t2, t3, t4, t5, t6, t7, t8, t9
        sage: gens = [t0, t1, t0*t3 - t1*t2, t4, t5, t6, t7, t8, t9]
        sage: p = t0*t3*t6 - t1*t2*t6 + t4*t5
        sage: is_element_of_subalgebra(gens, p)
        (True, t3*t4 + t2*t5)

    """
    # print("gens is ", gens, "\n", "p is ", p, "\n", "denom is ", denom)
    FF = p.parent().base_ring()

    # set up a polynomial ring with enough variables to refer to
    # gens and denoms as indeterminates
    if isinstance(Pt, bool):
        nr_gens = len(gens)
        names = ['t'+str(i) for i in range(nr_gens)]
        Pt = PolynomialRing(FF, nr_gens, names=names)
    else:
        nr_gens = Pt.ngens()

    if p == 1:
        return True, Pt(1)
    if p == 0:
        return True, Pt(0)

    # compute the leading monomials of the generators and the
    # denominators
    lm_gens = [x.lm() for x in gens]

    # extract info for the parent if p
    P = p.parent()
    gensP = P.gens()[::-1]

    # dict_gens is a dictionary that shows for the generators of P
    # in which element of gens it occurs in the leading term
    #
    # first list_gens contains a list which shows what is the highest
    # generator in the leading term of the k-th element of gens.
    list_gens = [max(x.variables()) for x in lm_gens]

    dict_gens = {list_gens[k]: k for k in range(len(list_gens))}
    # start the reduction. initialize tpol and den_mon
    tpol, den_exp = 0, 0
    if denom_var == 0:
        P_dn = PolynomialRing(Pt, 'dn')
        denom_var = P_dn.gens()[0]

    # we reduce p until it is zero
    while p != 0:
        # initialize the monomial to be added, the reduction pol and the
        # leading term of p.
        tmon, reduction_pol, lm_p = 1, 1, p.lm()

        # for each generator of P
        for g in gensP:
            # find the degree of lm_p in g. After several steps of reduction
            # lm_p might be a generalized monomial; that is, a monomial with
            # negative exponents.
            if lm_p.parent() is P:
                lm_p_num_degree = lm_p.degree(g)
                lm_p_den_degree = P.one().degree(g)
            else:
                lm_p_num_degree = P(lm_p.numerator()).degree(g)
                lm_p_den_degree = P(lm_p.denominator()).degree(g)
            lm_p_degree_g = lm_p_num_degree - lm_p_den_degree

            # lm_p_num_degree = P(lm_p.numerator()).degree(g)
            # lm_p_den_degree = P(lm_p.denominator()).degree(g)
            lm_p_degree_g = lm_p_num_degree - lm_p_den_degree

            # if g occurs with positive degree update the reduction pol
            # (which is a pol in the xi)
            # and tmon (monomial in the ti)
            if lm_p_degree_g > 0:
                reduction_pol *= gens[dict_gens[g]]**lm_p_degree_g
                tmon *= Pt.gens()[dict_gens[g]]**lm_p_degree_g
                # update also lm_p. this is when lm_p might become a
                # generalized monomial
                lm_p /= gens[dict_gens[g]].lm()**lm_p_degree_g

        # if p is a polynomial in gens, then then lm_p reduce to 1
        # in the previous cycle. If this is not the case, we must correct
        # by the denominators

        if lm_p != 1:
            # need to get denoms involved
            lm_p_denom = lm_p.denominator()

            # we try to reduce the denominator of lm_p with the
            # leading term of denom
            # if lm_p_denom is divisible by denom, then divide
            while lm_p.denominator() != 1:
                lm_p *= lm_p_denom
                p *= denom
                den_exp += 1

        # now compute the coefficient
        coeff = FF(-p.coefficient(p.lm()) *
                   reduction_pol.coefficient(reduction_pol.lm()))
        # modify tpol and also p
        tpol -= coeff*tmon/denom_var**den_exp
        p += coeff*reduction_pol
        # assert p.lm() < lm_old

    # return the final result, with numerator and denominator separately.
    if tpol.denominator() != 1:
        pass
        # breakpoint()
    tpol = tpol.subs({denom_var: denom_in_t})
    return True, tpol
