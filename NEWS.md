# User visible changes in `PRIMA.jl`

## Version 0.1.2

- General method `prima(f, x0; kwds...)` which solves the problem with the most
  suitable Powell's algorithm (among COBYLA, LINCOA, BOBYQA, or NEWUOA)
  depending on the constraints imposed via the keywords `kwds...`.

- Objective function and non-linear constraints for `cobyla`:
  - Non-linear equality constraints can be specified by keyword `nonlinear_eq`.
  - The objective function is called as `f(x)` like in other algorithms.
  - The functions implementing the non-linear constraints are passed by
    keywords `nonlinear_eq` and `nonlinear_ineq`.

- Algorithms have a more similar interface:
  - All algorithms have the same positional input, only the available keywords
    may change.
  - All algorithms have the same convention for the objective function.
    Non-linear constraints, if any, are provided by other user-defined
    functions.
  - All algorithms yield a 2-tuple `(x,info)` with `x` an approximate solution
    (provided algorithm was successful) and `info` a structured object
    collecting other outputs from the algorithm. The following properties are
    available for all algorithms: `info.fx` is the objective function value
    `f(x)`, `info.nf` is the number of evaluations of the objective function,
    and `info.status` is the termination status of the algorithm.
  - `issuccess(info)` yields whether algorithm was successful.
  - `PRIMA.reason(info)` yields a textual description of the termination status
    of the algorithm.

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
