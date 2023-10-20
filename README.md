# PRIMA [![Build Status](https://github.com/libprima/PRIMA.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/emmt/PRIMA.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Build Status](https://ci.appveyor.com/api/projects/status/github/emmt/PRIMA.jl?svg=true)](https://ci.appveyor.com/project/emmt/PRIMA-jl) [![Coverage](https://codecov.io/gh/emmt/PRIMA.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/emmt/PRIMA.jl)

This package is a Julia interface to the [PRIMA
library](https://github.com/libprima/prima), a *Reference Implementation for
Powell's Methods with Modernization and Amelioration*, by [Zaikun
Zhang](https://www.zhangzk.net/) who re-implemented and improved algorithms
originally by [M.J.D.
Powell](https://en.wikipedia.org/wiki/Michael_J._D._Powell) for minimizing a
multi-variate objective function possibly under constraints and without
derivatives.

Depending on the problem to solve, other Julia package(s) with similar
objectives may be of interest:

- [BlackBoxOptim](https://github.com/robertfeldt/BlackBoxOptim.jl) is a global
  optimization package for single- and multi-variate problems that does not
  require the objective function to be differentiable.

- [NOMAD](https://github.com/bbopt/NOMAD.jl) is an interface to the *Mesh
  Adaptive Direct Search algorithm* (MADS) designed for difficult blackbox
  optimization problems.

Formally, these algorithms are designed to solve problems of the form:

``` julia
min f(x)    subject to   x ∈ Ω ⊆ ℝⁿ
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

The five Powell's algorithms of the [PRIMA
library](https://github.com/libprima/prima) are provided by the `PRIMA`
package:

- `uobyqa` (*Unconstrained Optimization BY Quadratic Approximations*) is for
  unconstrained optimization, that is `Ω = ℝⁿ`.

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

All these algorithms are trust region methods where the variables are updated
according to an affine or a quadratic local approximation interpolating the
objective function at a given number of points (set by keyword `npt` by some of
the algorithms). No derivatives of the objective function are needed. These
algorithms are well suited to problems with a non-analytic objective function
that takes time to be evaluated.

The table below summarizes the characteristics of the different Powell's
methods (*"linear"* constraints includes equality and inequality linear
constraints).

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
x, info = uobyqa(f, x0; kwds...)
x, info = newuoa(f, x0; kwds...)
x, info = bobyqa(f, x0; kwds...)
x, info = cobyla(f, x0; kwds...)
x, info = lincoa(f, x0; kwds...)
```

where `f` is the objective function and `x0` specifies the initial values of
the variables (and is left unchanged). Constraints and options may be specified
by keywords `kwds...` (see below).

The objective function is called as `f(x)` with `x` the variables, it must
implement the following signature:

``` julia
f(x::Vector{Cdouble})::Real
```

All the algorithms return a 2-tuple `(x, info)` with `x` the variables and
`info` a structured object collecting all other information. If
`issuccess(info)` is true, then the algorithm was successful and `x` is an
approximate solution of the problem.

The output `info` has the following properties:

``` julia
info.fx       # value of the objective function f(x) on return
info.nf       # number of calls to the objective function
info.status   # final status code
info.cstrv    # amount of constraints violation, 0.0 if unconstrained
info.nl_eq    # non-linear equality constraints, empty vector if none
info.nl_ineq  # non-linear inequality constraints, empty vector if none
```

Calling one of:

``` julia
issuccess(info)
issuccess(info.status)
```

yield whether the algorithm has converged. If this is the case, `info.status`
can be one of:

- `PRIMA.SMALL_TR_RADIUS` if the radius of the trust region becomes smaller or
  equal the value of keyword `rhobeg`, in other words, the algorithm has
  converged in terms of variable precision;

- `PRIMA.FTARGET_ACHIEVED` if the objective function is smaller of equal the
  value of keyword `ftarget`, in other words, the algorithm has converged in
  terms of function value.

There are other possibilities which all indicate a failure. Calling one of:

``` julia
PRIMA.reason(info)
PRIMA.reason(info.status)
```

yield a textual explanation about the reason that leads the algorithm to stop.

The keywords allowed by the different algorithms are summarized by the
following table.

| Keyword          | Description                            | Algorithms                   |
|:-----------------|:---------------------------------------|:-----------------------------|
| `rhobeg`         | Initial trust region radius            | all                          |
| `rhoend`         | Final trust region radius              | all                          |
| `ftarget`        | Target objective function value        | all                          |
| `maxfun`         | Maximum number of function evaluations | all                          |
| `iprint`         | Verbosity level                        | all                          |
| `npt`            | Number of points in local model        | `bobyqa`, `lincoa`, `newuoa` |
| `xl`             | Lower bound                            | `bobyqa`, `cobyla`, `lincoa` |
| `xu`             | Upper bound                            | `bobyqa`, `cobyla`, `lincoa` |
| `nonlinear_eq`   | Non-linear equality constraints        | `cobyla`                     |
| `nonlinear_ineq` | Non-linear inequality constraints      | `cobyla`                     |
| `linear_eq`      | Linear equality constraints            | `cobyla`, `lincoa`           |
| `linear_ineq`    | Linear inequality constraints          | `cobyla`, `lincoa`           |

Assuming `n = length(x)` is the number of variables, then:

- `rhobeg` (default value `1.0`) is the initial radius of the trust region.

- `rhoend` (default value `1e-4*rhobeg`) is the final radius of the trust
  region. The algorithm stops when the trust region radius becomes smaller or
  equal `rhoend` and the status `PRIMA.SMALL_TR_RADIUS` is returned.

- `ftarget` (default value `-Inf`) is another convergence setting. The
  algorithm stops as soon as `f(x) ≤ ftarget` and the status
  `PRIMA.FTARGET_ACHIEVED` is returned.

- `iprint` (default value `PRIMA.MSG_NONE`) sets the level of verbosity of the
  algorithm. Possible values are `PRIMA.MSG_EXIT`, `PRIMA.MSG_RHO`, or
  `PRIMA.MSG_FEVL`.

- `maxfun` (default `100n`) is the maximum number of function evaluations
  allowed for the algorithm. If the number of calls to `f(x)` exceeds this
  value, the algorithm is stopped and the status `PRIMA.MAXFUN_REACHED` is
  returned.

- `npt` (default value `2n + 1`) is the number of points used to approximate
  the local behavior of the objective function and such that `n + 2 ≤ npt ≤
  (n + 1)*(n + 2)/2`. The default value corresponds to the one recommended by
  M.J.D. Powell.

- `xl` and `xu` (default `fill(+Inf, n)` and `fill(-Inf, n)`) are the
  element-wise lower and upper bounds for the variables. Feasible variables are
  such that `xl ≤ x ≤ xu`.

- `nonlinear_eq` (default `nothing`) may be specified with a function, say
  `c_eq`, implementing non-linear equality constraints defined by `c_eq(x) =
  0`. On return, the values of the non-linear equality constraints are given by
  `info.nl_eq` to save calling `c_eq(x)`.

- `nonlinear_ineq` (default `nothing`) may be specified with a function, say
  `c_ineq`, implementing non-linear inequality constraints defined by
  `c_ineq(x) ≤ 0`. On return, the values of the non-linear inequality
  constraints are given by `info.nl_ineq` to save calling `c_ineq(x)`.

- `linear_eq` (default `nothing`) may be specified as a tuple `(Aₑ,bₑ)` to
  impose the linear equality constraints `Aₑ⋅x = bₑ`.

- `linear_ineq` (default `nothing`) may be specified as a tuple `(Aᵢ,bᵢ)` to
  impose the linear inequality constraints `Aᵢ⋅x ≤ bᵢ`.


## References

The following 5 first references respectively describe Powell's algorithms
`cobyla`, `uobyqa`, `newuoa`, `bobyqa`, and `lincoa` while the last reference
provides a good and comprehensive introduction to these algorithms.

1. M.J.D. Powell, *"A direct search optimization method that models the
   objective and constraint functions by linear interpolation"* in Advances in
   Optimization and Numerical Analysis Mathematics and Its Applications, vol.
   **275**, pp. 51-67 (1994).

2. M.J.D. Powell, *"UOBYQA: unconstrained optimization by quadratic
   approximation"*, in Mathematical Programming, vol. **92**, pp. 555-582
   (2002) [DOI](https://doi.org/10.1007/s101070100290).

3. M.J.D. Powell, *"The NEWUOA software for unconstrained optimization without
   derivatives"*, in Nonconvex Optimization and Its Applications, pp. 255-297
   (2006).

4. M.J.D. Powell, *"The BOBYQA algorithm for bound constrained optimization
   without derivatives"*, Department of Applied Mathematics and Theoretical
   Physics, Cambridge, England, Technical Report NA2009/06 (2009).

5. M.J.D. Powell, *"On fast trust region methods for quadratic models with
   linear constraints"*, Technical Report of the Department of Applied
   Mathematics and Theoretical Physics, Cambridge University (2014).

6. T.M. Ragonneau & Z. Zhang, *"PDFO: A Cross-Platform Package for Powell's
   Derivative-Free Optimization Solvers"*,
   [arXiv](https://doi.org/10.48550/arXiv.2302.13246) (2023).
