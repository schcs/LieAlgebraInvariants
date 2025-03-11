# LieAlgebraInvariants

This module is contains some functionality to compute algebraically independent generators for the field of rational invariants for nilpotent Lie algebras. 

```python
sage: from examples_lie_algs import lie_algebra_upper_triangular_matrices
sage: from characteristics import rational_invariant_field
sage: l = lie_algebra_upper_triangular_matrices( 6, strict = True )
sage: l.rational_invariant_field()
[x16,
 x16*x25 - x15*x26,
 x16^2*x25*x34 - x16*x15*x26*x34 - x16^2*x24*x35 + x16*x26*x14*x35 + x16*x15*x36*x24 - x16*x14*x25*x36]

```
