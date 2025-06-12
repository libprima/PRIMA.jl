# User visible changes in `PRIMA.jl`

## Version 0.2.3 (2025-06-12)

### Fixed

-  Use `task_local_storage` for storing the reference to the Julia objective function.
   This replaces the storage indexed by `Threads.threadid()` which is wrong because, since
   Julia 1.8, a given task can be executed by different threads in its life time. In
   principle, this solves issue #37.

## Version 0.2.2 (2024-10-16)

### Changed

- Update compat and continuous integration workflow.


## Version 0.2.1 (2024-09-27)

### Added

- Use `CUTEst`.
- Add tests related to issue #19.

### Changed

- Simplify `check_bounds` as suggested in issue #24.
- Revise the default `maxfun` to `500*n` and default `rhoend` to `1.0e-6*rhobeg`. This
  aligns with Powell's settings in his tests and the default settings of the Fortran code.
  Anyway, `maxfun = 100*n` is too small for a DFO solver.

### Fixed

- Fix `issuccess` (issue #23).


## Version 0.2.0 (2023-10-23)

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

- Scaling of the variables: the variables `x ∈ ℝⁿ` may be scaled by the scaling
  factors provided by the caller via keywords `scale` to re-express the problem
  in the scaled variables `u ∈ ℝⁿ` such that `u[i] = x[i]/scale[i]`. Note that
  the objective function, the constraints (linear and non-linear), and the
  bounds remain specified in the variables. Scaling the variables is useful to
  improve the conditioning of the problem, to make the scaled variables `u`
  having approximately the same magnitude, and to adapt to heterogeneous
  variables or with different units.

## Version 0.1.1 (2023-10-17)

- Keywords for other constraints than bounds have been renamed as
  `nonlinear_ineq` for non-linear inequality constraints, `linear_ineq` for
  linear inequality constraints, and `linear_eq` for linear equality
  constraints.

- In `cobyla`, if the caller is not interested in the values of `c(x)`
  expressing the non-linear inequality constraints, the number of such
  constraints may be provided instead of the array to store `c(x)`.

## Version 0.1.0

First public version.
