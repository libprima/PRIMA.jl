# User visible changes in `PRIMA.jl`


## Version 0.1.1

- Keywords for other constraints than bounds have been renamed as
  `nonlinear_ineq` for non-linear inequality constraints, `linear_ineq` for
  linear inequality constraints, and `linear_eq` for linear equality
  constraints.

- In `cobyla`, if the caller is not interested in the values of `c(x)`
  expressing the non-linear inequality constraints, the number of such
  constraints may be provided instead of the array to store `c(x)`.

## Version 0.1.0

First public version.
