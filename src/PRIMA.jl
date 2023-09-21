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

# The high level wrappers.
for func in (:bobyqa, :newuoa, :uobyqa, :lincoa, :cobyla)
    func! = Symbol(func, "!")
    @eval begin
        @doc """
            $($func)(f, x0; kwds...) -> (x, fx, nf, rc)

        attempts to minimize objective function `f` by the
        [`PRIMA.$($func)!`](@ref) method starting with variables `x0` (which
        are left unchanged on output). The result is the 4-tuple `(x, fx, nf,
        rc)` with `x` the (approximate) solution found, `fx` the value of
        `f(x)`, `nf` the number of calls to `f`, and `rc` a status code (see
        [`PRIMA.Status`](@ref) and [`PRIMA.reason`](@ref)).

        See the documention of [`PRIMA.$($func)!`](@ref) for a description of
        the method and of the allowed keywords `kwds...`.

        """
        function $func(f, x0::AbstractVector{<:Real}; kwds...)
            x = copyto!(Vector{Cdouble}(undef, length(x0)), x0)
            return x, $func!(f, x; kwds...)...
        end
    end
end

"""
    PRIMA.bobyqa!(f, x; kwds...) -> (fx, nf, rc)

attempts to minimize objective function `f` by M.J.D. Powell's BOBYQA method
for "Bounded Optimization BY Quadratic Approximations" and without derivatives.
On entry, argument `x` is a dense vector of double precision value with the
initial variables; on return, `x` is overwritten by an approximate solution and
the 3-tuple `(fx, nf, rc)` is returned with, `fx` the value of `f(x)`, `nf` the
number of calls to `f`, and `rc` a status code (see [`PRIMA.Status`](@ref) and
[`PRIMA.reason`](@ref)).

Call [`PRIMA.bobyqa`](@ref) for a version that preserves the initial variables.

"""
function bobyqa!(f, x::DenseVector{Cdouble};
                 xl::DenseVector{Cdouble} = fill(typemin(Cdouble), length(x)),
                 xu::DenseVector{Cdouble} = fill(typemax(Cdouble), length(x)),
                 rhobeg::Real = 1.0,
                 rhoend::Real = rhobeg*1e-4,
                 iprint::Union{Integer,Message} = MSG_NONE,
                 ftarget::Real = -Inf,
                 maxfun::Integer = 100*length(x),
                 npt::Integer = 2*length(x) + 1)
    n = length(x) # number of variables
    @assert length(xl) == n
    @assert length(xu) == n
    @assert isfinite(rhobeg) && rhobeg > zero(rhobeg)
    @assert isfinite(rhoend) && rhoend ≥ zero(rhoend) && rhoend ≤ rhobeg
    @assert n + 1 ≤ npt ≤ (n + 1)*(n + 2)/2

    fx = Ref{Cdouble}(NaN) # to store f(x) on return
    nf = Ref{Cint}(0)      # to store number of function calls
    fw = ObjFun(f, n)      # wrapper to objective function
    fp = _push_wrapper(fw) # pointer to C-callable function
    try
        # Call low-level optimizer.
        rc = prima_bobyqa(calfun, n, x, f, xl, xu, nf, rhobeg, rhoend, ftarget, maxfun, npt, iprint)
        return (fx[], Int(nf[]), rc)
    finally
        _pop_wrapper(fw)
    end
end

"""
    PRIMA.newuoa!(f, x; kwds...) -> (fx, nf, rc)

attempts to minimize objective function `f` by M.J.D. Powell's NEWUOA method
for unconstrained optimization without derivatives. On entry, argument `x` is
a dense vector of double precision value with the initial variables; on return,
`x` is overwritten by an approximate solution and the 3-tuple `(fx, nf, rc)` is
returned with, `fx` the value of `f(x)`, `nf` the number of calls to `f`, and
`rc` a status code (see [`PRIMA.Status`](@ref) and [`PRIMA.reason`](@ref)).

Call [`PRIMA.newuoa`](@ref) for a version that preserves the
initial variables.

"""
function newuoa!(f, x::DenseVector{Cdouble};
                 rhobeg::Real = 1.0,
                 rhoend::Real = rhobeg*1e-4,
                 ftarget::Real = -Inf,
                 maxfun::Integer = 100*length(x),
                 npt::Integer = 2*length(x) + 1,
                 iprint::Union{Integer,Message} = MSG_NONE)
    n = length(x) # number of variables
    @assert isfinite(rhobeg) && rhobeg > zero(rhobeg)
    @assert isfinite(rhoend) && rhoend ≥ zero(rhoend) && rhoend ≤ rhobeg
    @assert n + 1 ≤ npt ≤ (n + 1)*(n + 2)/2

    fx = Ref{Cdouble}(NaN) # to store f(x) on return
    nf = Ref{Cint}(0)      # to store number of function calls
    fw = ObjFun(f, n)      # wrapper to objective function
    fp = _push_wrapper(fw) # pointer to C-callable function
    try
        # Call low-level optimizer.
        rc = prima_newuoa(calfun, n, x, f, nf, rhobeg, rhoend, ftarget, maxfun, npt, iprint)
        return (fx[], Int(nf[]), rc)
    finally
        _pop_wrapper(fw)
    end
end

"""
    PRIMA.uobyqa!(f, x; kwds...) -> (fx, nf, rc)

attempts to minimize objective function `f` by M.J.D. Powell's UOBYQA method
for "Unconstrained Optimization BY Quadratic Approximation" and without
derivatives. On entry, argument `x` is a dense vector of double precision value
with the initial variables; on return, `x` is overwritten by an approximate
solution and the 3-tuple `(fx, nf, rc)` is returned with, `fx` the value of
`f(x)`, `nf` the number of calls to `f`, and `rc` a status code (see
[`PRIMA.Status`](@ref) and [`PRIMA.reason`](@ref)).

Call [`PRIMA.uobyqa`](@ref) for a version that preserves the initial variables.

"""
function uobyqa!(f, x::DenseVector{Cdouble};
                 rhobeg::Real = 1.0,
                 rhoend::Real = rhobeg*1e-4,
                 ftarget::Real = -Inf,
                 maxfun::Integer = 100*length(x),
                 iprint::Union{Integer,Message} = MSG_NONE)
    n = length(x) # number of variables
    @assert isfinite(rhobeg) && rhobeg > zero(rhobeg)
    @assert isfinite(rhoend) && rhoend ≥ zero(rhoend) && rhoend ≤ rhobeg

    fx = Ref{Cdouble}(NaN) # to store f(x) on return
    nf = Ref{Cint}(0)      # to store number of function calls
    fw = ObjFun(f, n)      # wrapper to objective function
    fp = _push_wrapper(fw) # pointer to C-callable function
    try
        # Call low-level optimizer.
        rc = prima_uobyqa(calfun, n, x, f, nf, rhobeg, rhoend, ftarget, maxfun, iprint)
        return (fx[], Int(nf[]), rc)
    finally
        _pop_wrapper(fw)
    end
end

const LinearConstraint = Tuple{DenseMatrix{Cdouble},DenseVector{Cdouble}}

function cobyla!(f, x::DenseVector{Cdouble};
                 nlconstr::Union{DenseVector{Cdouble},Nothing} = nothing,
                 ineqconstr::Union{LinearConstraint,Nothing} = nothing,
                 eqconstr::Union{LinearConstraint,Nothing} = nothing,
                 xl::DenseVector{Cdouble} = fill(typemin(Cdouble), length(x)),
                 xu::DenseVector{Cdouble} = fill(typemax(Cdouble), length(x)),
                 rhobeg::Real = 1.0,
                 rhoend::Real = rhobeg*1e-4,
                 ftarget::Real = -Inf,
                 maxfun::Integer = 100*length(x),
                 iprint::Union{Integer,Message} = MSG_NONE)
    n = length(x) # number of variables
    @assert length(xl) == n
    @assert length(xu) == n
    @assert isfinite(rhobeg) && rhobeg > zero(rhobeg)
    @assert isfinite(rhoend) && rhoend ≥ zero(rhoend) && rhoend ≤ rhobeg

    m_nlcon, nlcon_ptr = _get_nonlinear_constraints(nlconstr)
    m_eq, A_eq, b_eq = _get_linear_constraints(eqconstr, n)
    m_ineq, A_ineq, b_ineq = _get_linear_constraints(ineqconstr, n)

    cstrv = Ref{Cdouble}(NaN)     # to store constraint violation
    fx = Ref{Cdouble}(NaN)        # to store f(x) on return
    nf = Ref{Cint}(0)             # to store number of function calls
    fw = ObjFunCon(f, n, m_nlcon) # wrapper to objective function
    fp = _push_wrapper(fw)        # pointer to C-callable function

    try
        # Call low-level optimizer.
        rc = GC.@preserve nlconstr ineqconstr eqconstr prima_cobyla(m_nlcon, calcfc,
             n, x, f, cstrv, nlconstr, m_ineq, Aineq, bineq, m_eq, Aeq, beq, xl, xu,
             nf, rhobeg, rhoend, ftarget, maxfun, iprint)
        return (fx[], Int(nf[]), rc, cstrv[])
    finally
        _pop_wrapper(fw)
    end
end

function lincoa!(f, x::DenseVector{Cdouble};
                 ineqconstr::Union{LinearConstraint,Nothing} = nothing,
                 eqconstr::Union{LinearConstraint,Nothing} = nothing,
                 xl::DenseVector{Cdouble} = fill(typemin(Cdouble), length(x)),
                 xu::DenseVector{Cdouble} = fill(typemax(Cdouble), length(x)),
                 rhobeg::Real = 1.0,
                 rhoend::Real = rhobeg*1e-4,
                 ftarget::Real = -Inf,
                 maxfun::Integer = 100*length(x),
                 npt::Integer = 2*length(x) + 1,
                 iprint::Union{Integer,Message} = MSG_NONE)
    n = length(x) # number of variables
    @assert length(xl) == n
    @assert length(xu) == n
    @assert isfinite(rhobeg) && rhobeg > zero(rhobeg)
    @assert isfinite(rhoend) && rhoend ≥ zero(rhoend) && rhoend ≤ rhobeg
    @assert n + 1 ≤ npt ≤ (n + 1)*(n + 2)/2

    m_eq, A_eq, b_eq = _get_linear_constraints(eqconstr, n)
    m_ineq, A_ineq, b_ineq = _get_linear_constraints(ineqconstr, n)

    cstrv = Ref{Cdouble}(NaN) # to store constraint violation
    fx = Ref{Cdouble}(NaN)    # to store f(x) on return
    nf = Ref{Cint}(0)         # to store number of function calls
    fw = ObjFun(f, n)         # wrapper to objective function
    fp = _push_wrapper(fw)    # pointer to C-callable function
    try
        # Call low-level optimizer.
        rc = GC.@preserve ineqconstr eqconstr prima_lincoa(calfun, n, x, f,
             cstrv, m_ineq, Aineq, bineq, m_eq, Aeq, beq, xl, xu, nf,
             rhobeg, rhoend, ftarget, maxfun, npt, iprint)
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
    f::Any  # user-defined onjective function
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
# a given type of objectve function.
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
_c_wrapper(::ObjFun) = @cfunction(_objfun_wrapper, Cvoid, (Ptr{Cdouble}, Ptr{Cdouble},))
_c_wrapper(::ObjFunCon) = @cfunction(_objfuncon_wrapper, Cvoid, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},))
function _push_wrapper(fw::AbstractObjFun)
    push!(_get_stack(fw), fw)
    return _c_wrapper(fw)
end
function _pop_wrapper(fw::AbstractObjFun)
    stack = _get_stack(fw)
    last(stack) === fw || error
    ("objective function is not the last one in the caller tread stask")
    resize!(stack, length(stack) - 1)
    return nothing
end

# FIXME: Matrix A of linear constraints is in row-major order (this is imposed
# by the C interface which transposes the matrix A).
_get_linear_constraints(::Nothing, n::Integer) = 0, Ptr{Cdouble}(0), Ptr{Cdouble}(0)
function _get_linear_constraints(Ab::LinearConstraint, n::Integer)
    A, b = Ab
    m = length(b) # number of constraints
    @assert size(A) == (n, m)
    return m, pointer(A), pointer(b)
end

_get_nonlinear_constraints(::Nothing, n::Integer) = 0, Ptr{Cdouble}(0)
_get_nonlinear_constraints(c::DenseVector{Cdouble}) =
    length(c), pointer(c)

end
