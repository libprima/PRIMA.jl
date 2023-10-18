module PRIMA

using PRIMA_jll
const libprimac = PRIMA_jll.libprimac
include("wrappers.jl")

export bobyqa, cobyla, lincoa, newuoa, uobyqa

using TypeUtils

#------------------------------------------------------------------------------
# PUBLIC INTERFACE

"""
    PRIMA.reason(rc) -> str

yields a textual message explaining `rc`, the code returned by one of the PRIMA
optimizers.

"""
reason(status::Union{Integer,Status}) = unsafe_string(prima_get_rc_string(status))

# The high level wrappers. First the methods, then their documentation.
for func in (:bobyqa, :newuoa, :uobyqa, :lincoa, :cobyla)
    func! = Symbol(func, "!")
    @eval begin
        function $func(f, x0::AbstractVector{<:Real}; kwds...)
            x = copyto!(Vector{Cdouble}(undef, length(x0)), x0)
            return x, $func!(f, x; kwds...)...
        end
    end
end

const _doc_2_inputs_4_outputs = """
The arguments are the objective function `f` and the initial variables `x0`
which are left unchanged on output. The returned value is the 4-tuple `(x, fx,
nf, rc)` with `x` the (approximate) solution, `fx` the value of `f(x)`, `nf`
the number of calls to `f`, and `rc` a status code (see [`PRIMA.Status`](@ref)
and [`PRIMA.reason`](@ref)).
"""

const _doc_2_inputs_5_outputs = """
The arguments are the objective function `f` and the initial variables `x0`
which are left unchanged on output. The returned value is the 5-tuple `(x, fx,
nf, rc, cstrv)` with `x` the (approximate) solution, `fx` the value of `f(x)`,
`nf` the number of calls to `f`, `rc` a status code (see [`PRIMA.Status`](@ref)
and [`PRIMA.reason`](@ref)), and `cstrv` is the amount of constraint violation.
"""

const _doc_common_keywords = """
- `rhobeg` (default value `1.0`) is the initial radius of the trust region.

- `rhoend` (default value `1e-4*rhobeg`) is the final radius of the trust
  region used to decide whether the algorithm has converged in the variables.

- `ftarget` (default value `-Inf`) is another convergence setting. The
  algorithm is stopped as soon as `f(x) ≤ ftarget` and the status
  `PRIMA.FTARGET_ACHIEVED` is returned.

- `maxfun` (default `100n`) is the maximum number of function evaluations
  allowed for the algorithm. If the number of calls to `f(x)` exceeds this
  value, the algorithm is stopped and the status `PRIMA.MAXFUN_REACHED` is
  returned.

- `iprint` (default value `PRIMA.MSG_NONE`) sets the level of verbosity of the
   algorithm. Possible values are `PRIMA.MSG_EXIT`, `PRIMA.MSG_RHO`, or
   `PRIMA.MSG_FEVL`.
"""

const _doc_npt = """
- `npt` (default value `2n + 1`) is the number of points used to approximate
  the local behavior of the objective function and such that `n + 2 ≤ npt ≤
  (n + 1)*(n + 2)/2`. The default value corresponds to the one recommended by
  M.J.D. Powell.
"""

const _doc_bound_constraints = """
- `xl` and `xu` (default `fill(+Inf, n)` and `fill(-Inf, n)`) are the
  elementwise lower and upper bounds for the variables. Feasible variables are
  such that `xl ≤ x ≤ xu` (elementwise).
"""

const _doc_nonlinear_constraints = """
- `nonlinear_ineq` (default `nothing`) may be specified with the number `m` of
   non-linear inequality constraints expressed `c(x) ≤ 0`. If the caller is
   interested in the values of `c(x)` at the returned solution, the keyword may
   be set with a vector of `m` double precision floating-point values to store
   `c(x)`.
"""

const _doc_linear_constraints = """
- `linear_eq` (default `nothing`) may be specified as a tuple `(Aₑ,bₑ)`
  to represent linear equality constraints. Feasible variables are
  such that `Aₑ⋅x = bₑ` holds elementwise.

- `linear_ineq` (default `nothing`) may be specified as a tuple `(Aᵢ,bᵢ)` to
  represent linear inequality constraints. Feasible variables are such that
  `Aᵢ⋅x ≤ bᵢ` holds elementwise.
"""

"""
    uobyqa(f, x0; kwds...) -> (x, fx, nf, rc)

approximately solves the unconstrained optimization problem:

    min f(x)    s.t.   x ∈ ℝⁿ

by M.J.D. Powell's UOBYQA (for \"Unconstrained Optimization BY Quadratic
Approximations\") method. This algorithm is based on a trust region method
where variables are updated according to a quadratic local approximation
interpolating the objective function. No derivatives of the objective function
are needed.

$(_doc_2_inputs_4_outputs)

The objective function takes a single argument, the variables `x`, and returns
the function value, it shall implement the following signature:

    f(x::Vector{Cdouble})::Real

Allowed keywords are (`n = length(x)` is the number of variables):

$(_doc_common_keywords)

""" uobyqa

"""
    newuoa(f, x0; kwds...) -> (x, fx, nf, rc)

approximately solves the unconstrained optimization problem:

    min f(x)    s.t.   x ∈ ℝⁿ

by M.J.D. Powell's NEWUOA method. This algorithm is based on a trust region
method where variables are updated according to a quadratic local approximation
interpolating the objective function at a number of `npt` points. No
derivatives of the objective function are needed.

$(_doc_2_inputs_4_outputs)

The objective function takes a single argument, the variables `x`, and returns
the function value, it shall implement the following signature:

    f(x::Vector{Cdouble})::Real

Allowed keywords are (`n = length(x)` is the number of variables):

$(_doc_common_keywords)

$(_doc_npt)

""" newuoa

"""
    bobyqa(f, x0; kwds...) -> (x, fx, nf, rc)

approximately solves the bound constrained optimization problem:

    min f(x)    s.t.   x ∈ Ω ⊆ ℝⁿ

with

   Ω = { x ∈ ℝⁿ | xl ≤ x ≤ xu }

by M.J.D. Powell's BOBYQA (for \"Bounded Optimization BY Quadratic
Approximations\") method. This algorithm is based on a trust region method
where variables are updated according to a quadratic local approximation
interpolating the objective function at a number of `npt` points. No
derivatives of the objective function are needed.

$(_doc_2_inputs_4_outputs)

The objective function takes a single argument, the variables `x`, and returns
the function value, it shall implement the following signature:

    f(x::Vector{Cdouble})::Real

Allowed keywords are (`n = length(x)` is the number of variables):

$(_doc_common_keywords)

$(_doc_npt)

$(_doc_bound_constraints)

""" bobyqa

"""
    cobyla(f, x0; kwds...) -> (x, fx, nf, rc, cstrv)

approximately solves the constrained optimization problem:

    min f(x)    s.t.   x ∈ Ω ⊆ ℝⁿ

with

    Ω = { x ∈ ℝⁿ | xl ≤ x ≤ xu, Aₑ⋅x = bₑ, Aᵢ⋅x ≤ bᵢ, and c(x) ≤ 0 }

by M.J.D. Powell's COBYLA (for \"Constrained Optimization BY Linear
Approximations\") method. This algorithm is based on a trust region method
where variables are updated according to a linear local approximation of the
objective function. No derivatives of the objective function are needed.

$(_doc_2_inputs_5_outputs)

The objective function takes two arguments, the variables `x` and a vector `cx`
to store the non-linear constraints, and returns the function value, it shall
implement the following signature:

    f(x::Vector{Cdouble}, cx::Vector{Cdouble})::Real

where the e,tries of `cx` are to be overwritten by the non-linear constraints
`c(x)`.

Allowed keywords are (`n = length(x)` is the number of variables):

$(_doc_common_keywords)

$(_doc_bound_constraints)

$(_doc_linear_constraints)

$(_doc_nonlinear_constraints)

""" cobyla

"""
    lincoa(f, x0; kwds...) -> (x, fx, nf, rc, cstrv)

approximately solves the constrained optimization problem:

    min f(x)    s.t.   x ∈ Ω ⊆ ℝⁿ

with

    Ω = { x ∈ ℝⁿ | xl ≤ x ≤ xu, Aₑ⋅x = bₑ, and Aᵢ⋅x ≤ bᵢ }

by M.J.D. Powell's LINCOA (for \"LINearly Constrained Optimization\") method.
This algorithm is based on a trust region method where variables are updated
according to a quadratic local approximation of the objective function. No
derivatives of the objective function are needed.

$(_doc_2_inputs_5_outputs)

The objective function takes a single argument, the variables `x`, and returns
the function value, it shall implement the following signature:

    f(x::Vector{Cdouble})::Real

Allowed keywords are (`n = length(x)` is the number of variables):

$(_doc_common_keywords)

$(_doc_npt)

$(_doc_bound_constraints)

$(_doc_linear_constraints)

""" lincoa

"""
    PRIMA.bobyqa!(f, x; kwds...) -> (fx, nf, rc)

in-place version of [`PRIMA.bobyqa`](@ref) which to see for details. On entry,
argument `x` is a dense vector of double precision value with the initial
variables; on return, `x` is overwritten by an approximate solution.

"""
function bobyqa!(f, x::DenseVector{Cdouble};
                 xl::Union{AbstractVector{<:Real},Nothing} = nothing,
                 xu::Union{AbstractVector{<:Real},Nothing} = nothing,
                 rhobeg::Real = 1.0,
                 rhoend::Real = 1e-4*rhobeg,
                 iprint::Union{Integer,Message} = MSG_NONE,
                 ftarget::Real = -Inf,
                 maxfun::Integer = 100*length(x),
                 npt::Integer = 2*length(x) + 1)
    # Check arguments and get constraints.
    n = length(x) # number of variables
    _check_rho(rhobeg, rhoend)
    _check_npt(npt, n)
    xl = _get_lower_bound(xl, n)
    xu = _get_upper_bound(xu, n)

    fx = Ref{Cdouble}(NaN) # to store f(x) on return
    nf = Ref{Cint}(0)      # to store number of function calls
    fw = ObjFun(f, n)      # wrapper to objective function
    fp = _push_wrapper(fw) # pointer to C-callable function
    try
        # Call low-level optimizer.
        rc = prima_bobyqa(fp, n, x, fx, xl, xu, nf, rhobeg, rhoend, ftarget, maxfun, npt, iprint)
        return (fx[], Int(nf[]), rc)
    finally
        _pop_wrapper(fw)
    end
end

"""
    PRIMA.newuoa!(f, x; kwds...) -> (fx, nf, rc)

in-place version of [`PRIMA.newuoa`](@ref) which to see for details. On entry,
argument `x` is a dense vector of double precision value with the initial
variables; on return, `x` is overwritten by an approximate solution.

"""
function newuoa!(f, x::DenseVector{Cdouble};
                 rhobeg::Real = 1.0,
                 rhoend::Real = 1e-4*rhobeg,
                 ftarget::Real = -Inf,
                 maxfun::Integer = 100*length(x),
                 npt::Integer = 2*length(x) + 1,
                 iprint::Union{Integer,Message} = MSG_NONE)
    # Check arguments.
    n = length(x) # number of variables
    _check_rho(rhobeg, rhoend)
    _check_npt(npt, n)

    fx = Ref{Cdouble}(NaN) # to store f(x) on return
    nf = Ref{Cint}(0)      # to store number of function calls
    fw = ObjFun(f, n)      # wrapper to objective function
    fp = _push_wrapper(fw) # pointer to C-callable function
    try
        # Call low-level optimizer.
        rc = prima_newuoa(fp, n, x, fx, nf, rhobeg, rhoend, ftarget, maxfun, npt, iprint)
        return (fx[], Int(nf[]), rc)
    finally
        _pop_wrapper(fw)
    end
end

"""
    PRIMA.uobyqa!(f, x; kwds...) -> (fx, nf, rc)

in-place version of [`PRIMA.uobyqa`](@ref) which to see for details. On entry,
argument `x` is a dense vector of double precision value with the initial
variables; on return, `x` is overwritten by an approximate solution.

"""
function uobyqa!(f, x::DenseVector{Cdouble};
                 rhobeg::Real = 1.0,
                 rhoend::Real = 1e-4*rhobeg,
                 ftarget::Real = -Inf,
                 maxfun::Integer = 100*length(x),
                 iprint::Union{Integer,Message} = MSG_NONE)
    # Check arguments.
    n = length(x) # number of variables
    _check_rho(rhobeg, rhoend)

    fx = Ref{Cdouble}(NaN) # to store f(x) on return
    nf = Ref{Cint}(0)      # to store number of function calls
    fw = ObjFun(f, n)      # wrapper to objective function
    fp = _push_wrapper(fw) # pointer to C-callable function
    try
        # Call low-level optimizer.
        rc = prima_uobyqa(fp, n, x, fx, nf, rhobeg, rhoend, ftarget, maxfun, iprint)
        return (fx[], Int(nf[]), rc)
    finally
        _pop_wrapper(fw)
    end
end

"""
    PRIME.LinearConstraints

is the type of `(A,b)`, the 2-tuple representing linear equality constraints
`A⋅x = b` or linear inequality constraints `A⋅x ≤ b` where `A` is a matrix, `x`
is the vector of variables, and `b` is a vector.

"""
const LinearConstraints = Tuple{AbstractMatrix{<:Real},AbstractVector{<:Real}}

"""
    PRIMA.cobyla!(f, x; kwds...) -> (fx, nf, rc)

in-place version of [`PRIMA.cobyla`](@ref) which to see for details. On entry,
argument `x` is a dense vector of double precision value with the initial
variables; on return, `x` is overwritten by an approximate solution.

"""
function cobyla!(f, x::DenseVector{Cdouble};
                 nonlinear_ineq::Union{AbstractVector{<:Real},Integer,Nothing} = nothing,
                 linear_ineq::Union{LinearConstraints,Nothing} = nothing,
                 linear_eq::Union{LinearConstraints,Nothing} = nothing,
                 xl::Union{AbstractVector{<:Real},Nothing} = nothing,
                 xu::Union{AbstractVector{<:Real},Nothing} = nothing,
                 rhobeg::Real = 1.0,
                 rhoend::Real = 1e-4*rhobeg,
                 ftarget::Real = -Inf,
                 maxfun::Integer = 100*length(x),
                 iprint::Union{Integer,Message} = MSG_NONE)
    # Check arguments and get constraints.
    n = length(x) # number of variables
    _check_rho(rhobeg, rhoend)
    xl = _get_lower_bound(xl, n)
    xu = _get_upper_bound(xu, n)
    m_nlcon, nlcon = _get_nonlinear_constraints(nonlinear_ineq)
    m_eq, A_eq, b_eq = _get_linear_constraints(linear_eq, n)
    m_ineq, A_ineq, b_ineq = _get_linear_constraints(linear_ineq, n)

    cstrv = Ref{Cdouble}(NaN)     # to store constraint violation
    fx = Ref{Cdouble}(NaN)        # to store f(x) on return
    nf = Ref{Cint}(0)             # to store number of function calls
    fw = ObjFunCon(f, n, m_nlcon) # wrapper to objective function
    fp = _push_wrapper(fw)        # pointer to C-callable function

    try
        # Call low-level optimizer.
        rc = prima_cobyla(m_nlcon, fp, n, x, fx, cstrv, nlcon,
                          m_ineq, A_ineq, b_ineq, m_eq, A_eq, b_eq,
                          xl, xu, nf, rhobeg, rhoend, ftarget, maxfun, iprint)
        return (fx[], Int(nf[]), rc, cstrv[])
    finally
        _pop_wrapper(fw)
    end
end

"""
    PRIMA.lincoa!(f, x; kwds...) -> (fx, nf, rc)

in-place version of [`PRIMA.lincoa`](@ref) which to see for details. On entry,
argument `x` is a dense vector of double precision value with the initial
variables; on return, `x` is overwritten by an approximate solution.

"""
function lincoa!(f, x::DenseVector{Cdouble};
                 linear_ineq::Union{LinearConstraints,Nothing} = nothing,
                 linear_eq::Union{LinearConstraints,Nothing} = nothing,
                 xl::Union{AbstractVector{<:Real},Nothing} = nothing,
                 xu::Union{AbstractVector{<:Real},Nothing} = nothing,
                 rhobeg::Real = 1.0,
                 rhoend::Real = 1e-4*rhobeg,
                 ftarget::Real = -Inf,
                 maxfun::Integer = 100*length(x),
                 npt::Integer = 2*length(x) + 1,
                 iprint::Union{Integer,Message} = MSG_NONE)
    # Check arguments and get constraints.
    n = length(x) # number of variables
    _check_rho(rhobeg, rhoend)
    _check_npt(npt, n)
    xl = _get_lower_bound(xl, n)
    xu = _get_upper_bound(xu, n)
    m_eq, A_eq, b_eq = _get_linear_constraints(linear_eq, n)
    m_ineq, A_ineq, b_ineq = _get_linear_constraints(linear_ineq, n)

    cstrv = Ref{Cdouble}(NaN) # to store constraint violation
    fx = Ref{Cdouble}(NaN)    # to store f(x) on return
    nf = Ref{Cint}(0)         # to store number of function calls
    fw = ObjFun(f, n)         # wrapper to objective function
    fp = _push_wrapper(fw)    # pointer to C-callable function
    try
        # Call low-level optimizer.
        rc = prima_lincoa(fp, n, x, fx, cstrv,
                          m_ineq, A_ineq, b_ineq, m_eq, A_eq, b_eq, xl, xu,
                          nf, rhobeg, rhoend, ftarget, maxfun, npt, iprint)
        return (fx[], Int(nf[]), rc, cstrv[])
    finally
        _pop_wrapper(fw)
    end
end

#------------------------------------------------------------------------------
# PRIVATE METHODS AND TYPES

abstract type AbstractObjFun end

# Small structure to wrap simple objective functions for BOBYQA, NEWUOA,
# UOBYQA, and LINCOA at low level.
struct ObjFun <: AbstractObjFun
    f::Any  # user-defined objective function
    nx::Int # number of variables
end

# Small structure to wrap objective functions with constraints for COBYLA at
# low level.
struct ObjFunCon <: AbstractObjFun
    f::Any  # user-defined objective function
    nx::Int # number of variables
    nc::Int # number of non-linear constraints
end

# Global variable storing the per-thread stacks of objective functions indexed
# by thread identifier and then by execution order. On start of an
# optimization, an object linked to the user-defined objective function is
# pushed. This object is popped out of the stack on return of the optimization
# call, whatever happens. It is therefore to wrap th call to the optimization
# method in a `try-finally` clause.
const _objfun_stack = Vector{Vector{ObjFun}}(undef, 0)
const _objfuncon_stack = Vector{Vector{ObjFunCon}}(undef, 0)

# Private function `_get_stack` yields the stack for the caller thread and for
# a given type of objective function.
_get_stack(fw::ObjFun) = _get_stack(_objfun_stack)
_get_stack(fw::ObjFunCon) = _get_stack(_objfuncon_stack)
function _get_stack(stack::Vector{Vector{T}}) where {T}
    i = Threads.threadid()
    while length(stack) < i
        push!(stack, Vector{T}(undef, 0))
    end
    return stack[i]
end

# The following wrappers are intended to be C-callable functions. Their
# respective signatures must follow the prototypes `prima_obj` and
# `prima_objcon` in `prima.h` header.
function _objfun_wrapper(x_ptr::Ptr{Cdouble},  # (input) variables
                         f_ptr::Ptr{Cdouble})  # (output) function value
    fw = last(_objfun_stack[Threads.threadid()]) # retrieve the objective function
    x = unsafe_wrap(Array, x_ptr, fw.nx)
    unsafe_store!(f_ptr, as(Cdouble, fw.f(x)))
    return nothing
end
function _objfuncon_wrapper(x_ptr::Ptr{Cdouble},  # (input) variables
                            f_ptr::Ptr{Cdouble},  # (output) function value
                            c_ptr::Ptr{Cdouble})  # (output) constraints
    fw = last(_objfuncon_stack[Threads.threadid()]) # retrieve the objective function
    x = unsafe_wrap(Array, x_ptr, fw.nx)
    c = unsafe_wrap(Array, c_ptr, fw.nc)
    unsafe_store!(f_ptr, as(Cdouble, fw.f(x, c)))
    return nothing
end

# Private functions `_push_wrapper` and `_pop_wrapper` are to be used in a
# `try-finally` clause as explained above.
_c_wrapper(::ObjFun) =
    @cfunction(_objfun_wrapper, Cvoid, (Ptr{Cdouble}, Ptr{Cdouble},))
_c_wrapper(::ObjFunCon) =
    @cfunction(_objfuncon_wrapper, Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},))
function _push_wrapper(fw::AbstractObjFun)
    push!(_get_stack(fw), fw)
    return _c_wrapper(fw)
end
function _pop_wrapper(fw::AbstractObjFun)
    stack = _get_stack(fw)
    last(stack) === fw || error(
        "objective function is not the last one in the caller thread stask")
    resize!(stack, length(stack) - 1)
    return nothing
end

function _check_rho(rhobeg, rhoend)
    isfinite(rhobeg) && rhobeg > zero(rhobeg) || throw(ArgumentError(
        "`rhobeg` must be finite and positive"))
    isfinite(rhoend) && rhoend ≥ zero(rhoend) || throw(ArgumentError(
        "`rhoend` must be finite and non-negative"))
    rhoend ≤ rhobeg|| throw(ArgumentError(
        "`rhoend` must be less or equal `rhobeg`"))
    nothing
end

function _check_npt(npt::Integer, n::Integer)
    n + 2 ≤ npt ≤ (n + 1)*(n + 2)/2 || throw(ArgumentError(
        "`n + 2 ≤ npt ≤ (n + 1)*(n + 2)/2` must hold"))
    nothing
end

@static if !isdefined(Base, :Returns)
    struct Returns{T}
        value::T
    end
    (obj::Returns)(args...; kwds...) = obj.value
end

# Null-array to represent missing argument. It is sufficient to implement the abstract
# array API plus the Base.unsafe_convert method.
struct NullArray{T,N} <: AbstractArray{T,N} end
const NullVector{T} = NullArray{T,1}
const NullMatrix{T} = NullArray{T,2}
Base.length(::NullArray{T,N}) where {T,N} = 0
Base.size(::NullArray{T,N}) where {T,N} = ntuple(Returns(0), Val(N))
Base.axes(::NullArray{T,N}) where {T,N} = ntuple(Returns(Base.OneTo(0)), Val(N))
Base.unsafe_convert(::Type{Ptr{S}}, ::NullArray{T,N}) where {T,N,S<:Union{Cvoid,T}} = Ptr{S}(0)

# FIXME: Matrix A of linear constraints is in row-major order (this is imposed
# by the C interface which transposes the matrix A).
_get_linear_constraints(::Nothing, n::Integer) =
    0, NullMatrix{Cdouble}(), NullVector{Cdouble}()
function _get_linear_constraints(Ab::LinearConstraints, n::Integer)
    A, b = Ab
    Base.has_offset_axes(A) && error(
        "matrix `A` of linear constraints must have 1-based indices")
    Base.has_offset_axes(b) && error(
        "vector `b` of linear constraints must have 1-based indices")
    m = length(b) # number of constraints
    size(A) == (m,n) || throw(DimensionMismatch(
        "matrix `A` of linear constraints has incompatible dimensions"))
    T = Cdouble
    # FIXME: Like in FORTRAN, Julia matrices are in column-major storage order,
    # but we must transpose the matrix A in linear constraints because we call
    # the FORTRAN code through a C interface which consider that matrices are
    # in row-major storage. As a result, the matrix `A` will be transposed
    # twice. This isn't a big issue for a small number of variables and
    # constraints, but it's not completely satisfactory either.
    A_ = Matrix{T}(undef, n, m)
    @inbounds for i ∈ 1:m
        for j ∈ 1:n
            A_[j,i] = A[i,j]
        end
    end
    b_ = _dense_array(T, b)
    return m, A_, b_
end

_get_nonlinear_constraints(::Nothing) =
    0, NullVector{Cdouble}()
_get_nonlinear_constraints(m::Integer) =
    _get_nonlinear_constraints(Vector{Cdouble}(undef, m))
_get_nonlinear_constraints(c::AbstractVector{<:Real}) =
    length(c), _dense_array(Cdouble, c)

for (uplo, def) in ((:lower, typemin),
                    (:upper, typemax))
    func = Symbol("_get_$(uplo)_bound")
    @eval begin
        $func(::Nothing, n::Integer) = fill($def(Cdouble), n)
        function $func(b::AbstractVector, n::Integer)
            Base.has_offset_axes(b) && error(
                $("$uplo bound must have 1-based indices"))
            length(b) == n || throw(DimensionMismatch(string(
                $("$uplo bound must have "), n, " elements")))
            return _dense_array(Cdouble, b)
        end
    end
end

# Yields an array with given element type and which can be safely passed to a C
# function (i.e. it has contiguous elements in memory).
_dense_array(::Type{T}, A::DenseArray{T,N}) where {T,N} = A
_dense_array(::Type{T}, A::AbstractArray{<:Any,N}) where {N,T} = convert(Array{T,N}, A)

end
