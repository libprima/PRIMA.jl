# PRIMA [![Build Status](https://github.com/emmt/PRIMA.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/emmt/PRIMA.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Build Status](https://ci.appveyor.com/api/projects/status/github/emmt/PRIMA.jl?svg=true)](https://ci.appveyor.com/project/emmt/PRIMA-jl) [![Coverage](https://codecov.io/gh/emmt/PRIMA.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/emmt/PRIMA.jl)

This package is a Julia interface to the
[PRIMA](https://github.com/libprima/prima) library, a **R**eference
**I**mplementation for **P**owell's methods with **M**odernization and
**A**melioration, which re-implements algorithms originally by M.J.D. Powell
for minimizing a multi-variate objective function possibly under constraints
and without derivatives.

Formally, these algorithms are designed to solve problems of the form:

``` julia
min f(x)    s.t.   x ∈ Ω ⊆ ℝⁿ
```

where `f: Ω → ℝ` is the function to minimize, `Ω ⊆ ℝⁿ` is the set of feasible
variables, and `n ≥ 1` is the number of variables. The most general feasible
set is:

``` julia
Ω = { x ∈ ℝⁿ | xl ≤ x ≤ xu, Aₑ⋅x = bₑ, Aᵢ⋅x ≤ bᵢ, and c(x) ≤ 0 }
```

where `xl ∈ ℝⁿ` and `xu ∈ ℝⁿ` are lower and upper bounds, `Aₑ` and `bₑ`
implement linear equality constraints, `Aᵢ` and `bᵢ` implement linear
inequality constraints, and `c: ℝⁿ → ℝᵐ` implements `m` non-linear constraints.

Five algorithms are provided by the `PRIMA` package:

- `uobyqa` (*Unconstrained Optimization BY Quadratic Approximations*) is for
  unconstrained optimization, that is `Ω = ℝⁿ`. According to M.J.D. Powell,
  `newuoa` is superior to `uobyqa`.

- `newuoa` is also for unconstrained optimization. According to M.J.D. Powell,
  `newuoa` is superior to `uobyqa`.

- `bobyqa` (*Bounded Optimization BY Quadratic Approximations*) is for simple
  bound constrained problems, that is `Ω = { x ∈ ℝⁿ | xl ≤ x ≤ xu }`.

- `lincoa` (*LINearly Constrained Optimization*) is for constrained
  optimization problems with bound constraints, linear equality constraints,
  and linear inequality constraints.

- `cobyla` (*Constrained Optimization BY Linear Approximations*) is for general
  constrained problems with bound constraints, non-linear constraints, linear
  equality constraints, and linear inequality constraints.

All these algorithms are based on trust region strategies and determine the
change of variables according to an affine or quadratic local approximation of
the objective function. This approximation interpolates the objective function
at a given number of points (set by keyword `npt` by some of the algorithms).
No derivatives of the objective function are needed. These algorithms are well
suited to problems with a non-analytic objective function that takes time to be
evaluated.

The table below summarizes the characteristics of the different Powell's
methods, *"linear"* constraints includes equality and inequality linear
constraints.

| Method   | Model     | Constraints                |
|:---------|:----------|:---------------------------|
| `newuoa` | quadratic | none                       |
| `uobyqa` | quadratic | none                       |
| `bobyqa` | quadratic | bounds                     |
| `lincoa` | quadratic | bounds, linear             |
| `cobyla` | affine    | bounds, linear, non-linear |

These methods are called as follows:

``` julia
using PRIMA
x, fx, nf, rc = uobyqa(f, x0; kwds...)
x, fx, nf, rc = newuoa(f, x0; kwds...)
x, fx, nf, rc = bobyqa(f, x0; kwds...)
x, fx, nf, rc, cstrv = cobyla(f, x0; kwds...)
x, fx, nf, rc, cstrv = lincoa(f, x0; kwds...)
```

where `f` is the objective function and `x0` is the initial solution.
Constraints and options may be specified by keywords `kwds...` (see below). The
result is the 4-tuple `(x, fx, nf, rc)` or the 5-tuple `(x, fx, nf, rc, cstrv)`
where `x` is the (approximate) solution found by the algorithm, `fx` is the
value of `f(x)`, `nf` is the number of calls to `f`, `rc` is a status code (an
enumeration value of type `PRIMA.Status`), and `cstrv` is the amount of
constraint violation.

For example, `rc` can be:

- `PRIMA.SMALL_TR_RADIUS` if the radius of the trust region becomes smaller or
  equal the value of keyword `rhobeg`, in other words, the algorithm has
  converged in terms of variable precision;

- `PRIMA.FTARGET_ACHIEVED` if the objective function is smaller of equal the
  value of keyword `ftarget`, in other words, the algorithm has converged in
  terms of function value.

There are other possibilities which all indicate an error. Calling:

``` julia
PRIMA.reason(rc)
```

yields a textual explanation about the reason that leads the algorithm to stop.

For all algorithms, except `cobyla`, the user-defined function takes a single
argument, the variables of the problem, and returns the value of the objective
function. It has the following signature:

``` julia
function objfun(x::Vector{Cdouble})
    return f(x)
end
```

The `cobyla` algorithm can account for `m` non-linear constraints expressed by
`c(x) ≤ 0`. For this algorithm, the user-defined function takes two arguments,
the `x` variables `x` and the values taken by the constraints `cx`, it shall
overwrite the array of constraints with `cx = c(x)` and return the value of the
objective function. It has the following signature:

``` julia
function objfuncon(x::Vector{Cdouble}, cx::Vector{Cdouble})
    copyto!(cx, c(x)) # set constraints
    return f(x)       # return value of objective function
end
```

The keywords allowed by the different algorithms are summarized by the
following table.

| Keyword      | Description                            | Algorithms                   |
|:-------------|:---------------------------------------|:-----------------------------|
| `rhobeg`     | Initial trust region radius            | all                          |
| `rhoend`     | Final trust region radius              | all                          |
| `ftarget`    | Target objective function value        | all                          |
| `maxfun`     | Maximum number of function evaluations | all                          |
| `iprint`     | Verbosity level                        | all                          |
| `npt`        | Number of points in local model        | `bobyqa`, `lincoa`, `newuoa` |
| `xl`         | Lower bound                            | `bobyqa`, `cobyla`, `lincoa` |
| `xu`         | Upper bound                            | `bobyqa`, `cobyla`, `lincoa` |
| `nlconstr`   | Non-linear constraints                 | `cobyla`                     |
| `eqconstr`   | Linear equality constraints            | `cobyla`, `lincoa`           |
| `ineqconstr` | Linear inequality constraints          | `cobyla`, `lincoa`           |

Assuming `n = length(x)` is the number of variables, then:

- `rhobeg` (default value `1.0`) is the initial radius of the trust region.

- `rhoend` (default value `rhobeg*1e-4`) is the final radius of the trust
  region used to decide whether the algorithm has converged in the variables.

- `ftarget` (default value `-Inf`) is another convergence setting. The
  algorithm is stopped as soon as `f(x) ≤ ftarget` and the status
  `PRIMA.FTARGET_ACHIEVED` is returned.

- `iprint` (default value `PRIMA.MSG_NONE`) sets the level of verbosity of the
   algorithm. Possible values are `PRIMA.MSG_EXIT`, `PRIMA.MSG_RHO`, or
   `PRIMA.MSG_FEVL`.

- `maxfun` (default `100n`) is the maximum number of function evaluations
  allowed for the algorithm. If the number of calls to `f(x)` exceeds this
  value, the algorithm is stopped and the status `PRIMA.MAXFUN_REACHED` is
  returned.

- `npt` (default value `2n + 1`) is the number of points used to approximate
  the local behavior of the objective function and such that `n + 1 ≤ npt ≤
  (n + 1)*(n + 2)/2`. The default value corresponds to the one recommended by
  M.J.D. Powell.

- `xl` and `xu` (default `fill(+Inf, n)` and `fill(-Inf, n)`) are the
  elementwise lower and upper bounds for the variables. Feasible variables are
  such that `xl ≤ x ≤ xu` (elementwise).

- `nlconstr` (default `nothing`) may be specified as a vector of `m` double
  precision floating-point values which are passed to the user-defined function
  to store `c(x)` the non-linear constraints in `x`. This keyword only exists
  for `cobyla`.

- `eqconstr` (default `nothing`) may be specified as a tuple `(A,b)`
  to represent linear equality constraints. Feasible variables are
  such that `A'⋅x = b` holds elementwise.

- `neqconstr` (default `nothing`) may be specified as a tuple `(A,b)` to
  represent linear inequality constraints. Feasible variables are such that
  `A'⋅x ≤ b` holds elementwise.
