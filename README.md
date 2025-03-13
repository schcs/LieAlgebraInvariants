# LieAlgebraInvariants

This [SageMath](https://www.sagemath.org) module contains some functionality to compute algebraically independent generators for the field of rational invariants for nilpotent Lie algebras. 

```python
sage: from examples_lie_algs import lie_algebra_upper_triangular_matrices
sage: from rational_invariants import rational_invariant_field
sage: l = lie_algebra_upper_triangular_matrices( 6, strict = True )
sage: l.rational_invariant_field()
[x16,
 x16*x25 - x15*x26,
 x16^2*x25*x34 - x16*x15*x26*x34 - x16^2*x24*x35 + x16*x26*x14*x35 + x16*x15*x36*x24 - x16*x14*x25*x36]
sage: from examples_lie_algs import standard_filiform_lie_algebra
sage: l = standard_filiform_lie_algebra( 10, field = QQ )
sage: l.rational_invariant_field()
[y0,
 y0*y2 - 1/2*y1^2,
 y0^2*y3 - y0*y1*y2 + 1/3*y1^3,
 y0^3*y4 - y0^2*y1*y3 + 1/2*y0*y1^2*y2 - 1/8*y1^4,
 y0^4*y5 - y0^3*y1*y4 + 1/2*y0^2*y1^2*y3 - 1/6*y0*y1^3*y2 + 1/30*y1^5,
 y0^5*y6 - y0^4*y1*y5 + 1/2*y0^3*y1^2*y4 - 1/6*y0^2*y1^3*y3 + 1/24*y0*y1^4*y2 - 1/144*y1^6,
 y0^6*y7 - y0^5*y1*y6 + 1/2*y0^4*y1^2*y5 - 1/6*y0^3*y1^3*y4 + 1/24*y0^2*y1^4*y3 - 1/120*y0*y1^5*y2 + 1/840*y1^7,
 y0^7*y8 - y0^6*y1*y7 + 1/2*y0^5*y1^2*y6 - 1/6*y0^4*y1^3*y5 + 1/24*y0^3*y1^4*y4 - 1/120*y0^2*y1^5*y3 + 1/720*y0*y1^6*y2 - 1/5760*y1^8]
sage: from examples_lie_algs import lie_alg_family_6_19
sage: l = lie_alg_family_6_19( QQ )
sage: l
Lie algebra on 6 generators (x1, x2, x3, x4, x5, x6) over Rational function field in ε over Rational Field
sage: l.rational_invariant_field()
[x1, x1*x6 + (-ε)*x1*x4 + (-1/2)*x3^2 + (-1/2*ε)*x2^2]
```

In the output, we obtain three algebraically independent generators for the field of rational invariants of the Lie algebra of `6x6` strictly upper triangular matrices.