module PRIMANLPModelsExt

if isdefined(Base, :get_extension)
    using PRIMA, NLPModels
else
    using ..PRIMA, ..NLPModels
end

# This structure is to wrap a non-linear problem model into a callable object.
struct ObjectiveFunction{F<:AbstractNLPModel} <: Function
    nlp::F
end
(f::ObjectiveFunction)(x::AbstractVector) = NLPModels.obj(f.nlp, x)

function check_variables(nlp::AbstractNLPModel, x0::AbstractVector)
    length(x0) == get_nvar(nlp) || error(
        "initial variables must have $(get_nvar(nlp)) elements")
    return x0
end

const Variables = AbstractVector{<:Real}

for func in (:uobyqa, :newuoa)
    @eval function PRIMA.$func(nlp::AbstractNLPModel, x0::Variables = get_x0(nlp); kwds...)
        has_bounds(nlp) && error("`$($func)` cannot solve problems with bound constraints")
        has_equalities(nlp) && error("`$($func)` cannot solve problems with equality constraints")
        has_inequalities(nlp) && error("`$($func)` cannot solve problems with inequality constraints")
        return uobyqa(ObjectiveFunction(nlp), check_variables(nlp, x0);  kwds...)
    end
end

function PRIMA.bobyqa(nlp::AbstractNLPModel, x0::Variables = get_x0(nlp); kwds...)
    has_equalities(nlp) && error("`bobyqa` cannot solve problems with equality constraints")
    has_inequalities(nlp) && error("`bobyqa` cannot solve problems with inequality constraints")
    return bobyqa(ObjectiveFunction(nlp), check_variables(nlp, x0); kwds...,
                  xl = get_lvar(nlp), xu = get_uvar(nlp))
end

function PRIMA.lincoa(nlp::AbstractNLPModel, x0::Variables = get_x0(nlp); kwds...)
    nlin = get_nlin(nlp) # number of linear constraints
    nnln = get_nnln(nlp) # number of non-linear constraints
    nlin == 0 || error("linear constraints not yet implemented for NLPModels in `lincoa`")
    nnln == 0 || error("`lincoa` cannot solve problems with non-linear constraints")
    return lincoa(ObjectiveFunction(nlp), check_variables(nlp, x0); kwds...,
                  xl = get_lvar(nlp), xu = get_uvar(nlp),
                  linear_eq = nothing, linear_ineq = nothing)
end

function PRIMA.cobyla(nlp::AbstractNLPModel, x0::Variables = get_x0(nlp); kwds...)
    nlin = get_nlin(nlp) # number of linear constraints
    nnln = get_nnln(nlp) # number of non-linear constraints
    nlin == 0 || error("linear constraints not yet implemented for NLPModels in `cobyla`")
    nnln == 0 || error("non-linear constraints not yet implemented for NLPModels in `cobyla`")
    return cobyla(ObjectiveFunction(nlp), check_variables(nlp, x0); kwds...,
                  xl = get_lvar(nlp), xu = get_uvar(nlp),
                  linear_eq = nothing, linear_ineq = nothing,
                  nonlinear_eq = nothing, nonlinear_ineq = nothing)
end

function PRIMA.prima(nlp::AbstractNLPModel, x0::Variables = get_x0(nlp); kwds...)
    nlin = get_nlin(nlp) # number of linear constraints
    nnln = get_nnln(nlp) # number of non-linear constraints
    if nnln > 0
        # Only COBYLA can deal with non-linear constraints.
        error("solving problem with non-linear constraints by COBYLA not yet implemented")
        return cobyla(nlp, x0; kwds...)
    elseif nlin > 0
        # LINCOA can deal with bounds and linear constraints.
        error("solving problem with linear constraints by LINCOA not yet implemented")
        return lincoa(nlp, x0; kwds...)
    elseif has_bounds(nlp)
        # BOBYQA can deal with bounds.
        return bobyqa(nlp, x0; kwds...)
    else
        # Use NEWUOA for unconstrained problems.
        return newuoa(nlp, x0; kwds...)
    end
end

end # module
