# typedef enum prima_message_t;
@enum Message::Cint begin
    MSG_NONE = 0
    MSG_EXIT = 1
    MSG_RHO = 2
    MSG_FEVL = 3
end

# typedef enum prima_rc_t;
@enum Status::Cint begin
    SMALL_TR_RADIUS = 0
    FTARGET_ACHIEVED = 1
    TRSUBP_FAILED = 2
    MAXFUN_REACHED = 3
    MAXTR_REACHED = 20
    NAN_INF_X = -1
    NAN_INF_F = -2
    NAN_INF_MODEL = -3
    NO_SPACE_BETWEEN_BOUNDS = 6
    DAMAGING_ROUNDING = 7
    ZERO_LINEAR_CONSTRAINT = 8
    CALLBACK_TERMINATE = 30
    INVALID_INPUT = 100
    ASSERTION_FAILS = 101
    VALIDATION_FAILS = 102
    MEMORY_ALLOCATION_FAILS = 103
    NULL_OPTIONS = 110
    NULL_PROBLEM = 111
    NULL_X0 = 112
    NULL_RESULT = 113
    NULL_FUNCTION = 114
    PROBLEM_SOLVER_MISMATCH_NONLINEAR_CONSTRAINTS = 115
    PROBLEM_SOLVER_MISMATCH_LINEAR_CONSTRAINTS = 116
    PROBLEM_SOLVER_MISMATCH_BOUNDS = 117
end

function prima_get_rc_string(rc)
    @ccall libprimac.prima_get_rc_string(rc::Cint)::Cstring
end

"""
    PrimaObj

is the Julia type for:

```c
typedef void (*prima_obj_t)(const double x[], double *const f,
                            const void *data);
```

the objective function required by UOBYQA, NEWUOA, BOBYQA, and LINCOA with
arguments:

* `x`: on input, the vector of variables (should not be modified);

* `f`: on output, the value of the function; a `NaN` value can be passed to
  signal an evaluation error;

* `data`: user data.

"""
const PrimaObj = Ptr{Cvoid}

"""
    PrimaObjcon

is the Julia type for:

```c
typedef void (*prima_objcon_t)(const double x[], double *const f,
                               double constr[], const void *data);
```

the objective & constraint function required by COBYLA with arguments:

* `x`: on input, the vector of variables (should not be modified)

* `f`: on output, the value of the function a `NaN` value can be passed to
  signal an evaluation error

* `constr`: on output, the value of the constraints (of size `m_nlcon`), with
  the constraints being `constr ≤ 0` `NaN` values can be passed to signal
  evaluation errors;

* `data`: user data.

"""
const PrimaObjcon = Ptr{Cvoid}

"""
    PrimaCallback

is the Julia type for:

```c
typedef void (*prima_callback_t)(const int n, const double x[], const double f,
                                 const int nf, const int tr, const double cstrv,
                                 const int m_nlcon, const double nlconstr[],
                                 bool *const terminate);
```

the callback function to report algorithm progress with arguments:

* `n`: number of variables;

* `x`: the current best point;

* `f`: the function value of the current best point;

* `nf`: number of function evaluations;

* `tr`: number of trust-region iterations;

* `cstrv`: the constraint violation of the current best point (LINCOA and
  COBYLA only);

* `m_nlcon`: number of nonlinear constraints (COBYLA only);

* `nlconstr`: nonlinear constraint values of the current best point (COBYLA
  only);

* `terminate`: a boolean to ask for termination.

"""
const PrimaCallback = Ptr{Cvoid}

const UserData = Ptr{Cvoid}

"""
```julia
prima_bobyqa(calfun, data, n, x, f, xl, xu, nf, rhobeg, rhoend, ftarget,
             maxfun, npt, iprint, callback, info) -> status
```

Julia wrapper for the C function:

```c
int bobyqa_c(prima_obj_t calfun, const void *data, const int n, double x[],
             double *const f, const double xl[], const double xu[],
             int *const nf, const double rhobeg, const double rhoend,
             const double ftarget, const int maxfun, const int npt, const int iprint,
             const prima_callback_t callback, int *const info);
```

"""
function prima_bobyqa(calfun, data, n, x, f, xl, xu, nf, rhobeg, rhoend, ftarget,
                      maxfun, npt, iprint, callback, info)
    @ccall libprimac.bobyqa_c(calfun::PrimaObj, data::UserData, n::Cint, x::Ptr{Cdouble},
                              f::Ptr{Cdouble}, xl::Ptr{Cdouble}, xu::Ptr{Cdouble},
                              nf::Ptr{Cint}, rhobeg::Cdouble, rhoend::Cdouble,
                              ftarget::Cdouble, maxfun::Cint, npt::Cint, iprint::Cint,
                              callback::PrimaCallback, info::Ptr{Cint})::Status
end

"""
```julia
prima_newuoa(calfun, data, n, x, f, nf, rhobeg, rhoend, ftarget, maxfun, npt,
             iprint, callback, info) -> status
```

is the Julia wrapper for the C function:

```c
int newuoa_c(prima_obj_t calfun, const void *data, const int n, double x[],
             double *const f, int *const nf, const double rhobeg, const double rhoend,
             const double ftarget, const int maxfun, const int npt, const int iprint,
             const prima_callback_t callback, int *const info);
```

"""
function prima_newuoa(calfun, data, n, x, f, nf, rhobeg, rhoend, ftarget, maxfun, npt,
                      iprint, callback, info)
    @ccall libprimac.newuoa_c(calfun::PrimaObj, data::UserData, n::Cint, x::Ptr{Cdouble},
                              f::Ptr{Cdouble}, nf::Ptr{Cint}, rhobeg::Cdouble,
                              rhoend::Cdouble, ftarget::Cdouble, maxfun::Cint,
                              npt::Cint, iprint::Cint,
                              callback::PrimaCallback, info::Ptr{Cint})::Status
end

"""
```julia
prima_uobyqa(calfun, data, n, x, f, nf, rhobeg, rhoend, ftarget, maxfun, iprint,
             callback, info) -> status
```

is the Julia wrapper for the C function:

```c
int uobyqa_c(prima_obj_t calfun, const void *data, const int n, double x[],
             double *const f, int *const nf, const double rhobeg, const double rhoend,
             const double ftarget, const int maxfun, const int iprint,
             const prima_callback_t callback, int *const info);
```

"""
function prima_uobyqa(calfun, data, n, x, f, nf, rhobeg, rhoend, ftarget, maxfun, iprint,
                      callback, info)
    @ccall libprimac.uobyqa_c(calfun::PrimaObj, data::UserData, n::Cint, x::Ptr{Cdouble},
                              f::Ptr{Cdouble}, nf::Ptr{Cint}, rhobeg::Cdouble,
                              rhoend::Cdouble, ftarget::Cdouble, maxfun::Cint, iprint::Cint,
                              callback::PrimaCallback, info::Ptr{Cint})::Status
end


"""
```julia
prima_cobyla(m_nlcon, calcfc, data, n, x, f, cstrv, nlconstr, m_ineq, Aineq, bineq,
             m_eq, Aeq, beq, xl, xu, f0, nlconstr0, nf, rhobeg, rhoend, ftarget,
             maxfun, iprint, ctol, callback, info) -> status
```

is the Julia wrapper for the C function:

```c
int cobyla_c(const int m_nlcon, const prima_objcon_t calcfc, const void *data, const int n,
             double x[], double *const f, double *const cstrv, double nlconstr[],
             const int m_ineq, const double Aineq[], const double bineq[],
             const int m_eq, const double Aeq[], const double beq[],
             const double xl[], const double xu[],
             const double f0, const double nlconstr0[],
             int *const nf, const double rhobeg, const double rhoend, const double ftarget,
             const int maxfun, const int iprint, const double ctol,
             const prima_callback_t callback, int *const info);
```

"""
function prima_cobyla(m_nlcon, calcfc, data, n, x, f, cstrv, nlconstr, m_ineq, Aineq, bineq,
                      m_eq, Aeq, beq, xl, xu, f0, nlconstr0, nf, rhobeg, rhoend, ftarget,
                      maxfun, iprint, ctol, callback, info)
    @ccall libprimac.cobyla_c(m_nlcon::Cint, calcfc::PrimaObjcon, data::UserData, n::Cint,
                              x::Ptr{Cdouble}, f::Ptr{Cdouble}, cstrv::Ptr{Cdouble},
                              nlconstr::Ptr{Cdouble}, m_ineq::Cint, Aineq::Ptr{Cdouble},
                              bineq::Ptr{Cdouble}, m_eq::Cint, Aeq::Ptr{Cdouble},
                              beq::Ptr{Cdouble}, xl::Ptr{Cdouble}, xu::Ptr{Cdouble},
                              f0::Cdouble, nlconstr0::Ptr{Cdouble},
                              nf::Ptr{Cint}, rhobeg::Cdouble, rhoend::Cdouble,
                              ftarget::Cdouble, maxfun::Cint, iprint::Cint, ctol::Cdouble,
                              callback::PrimaCallback, info::Ptr{Cint})::Status
end


"""
```julia
prima_lincoa(calfun, data, n, x, f, cstrv, m_ineq, Aineq, bineq, m_eq, Aeq, beq,
             xl, xu, nf, rhobeg, rhoend, ftarget, maxfun, npt, iprint, ctol,
             callback, info)
```

is the Julia wrapper for the C function:

```c
int lincoa_c(prima_obj_t calfun, const void *data, const int n, double x[],
             double *const f, double *const cstrv,
             const int m_ineq, const double Aineq[], const double bineq[],
             const int m_eq, const double Aeq[], const double beq[],
             const double xl[], const double xu[],
             int *const nf, const double rhobeg, const double rhoend,
             const double ftarget, const int maxfun, const int npt,
             const int iprint, const double ctol,
             const prima_callback_t callback, int *const info);
```

"""
function prima_lincoa(calfun, data, n, x, f, cstrv, m_ineq, Aineq, bineq, m_eq, Aeq, beq,
                      xl, xu, nf, rhobeg, rhoend, ftarget, maxfun, npt, iprint, ctol,
                      callback, info)
    @ccall libprimac.lincoa_c(calfun::PrimaObj, data::UserData, n::Cint, x::Ptr{Cdouble},
                              f::Ptr{Cdouble}, cstrv::Ptr{Cdouble},
                              m_ineq::Cint, Aineq::Ptr{Cdouble}, bineq::Ptr{Cdouble},
                              m_eq::Cint, Aeq::Ptr{Cdouble}, beq::Ptr{Cdouble},
                              xl::Ptr{Cdouble}, xu::Ptr{Cdouble},
                              nf::Ptr{Cint}, rhobeg::Cdouble, rhoend::Cdouble,
                              ftarget::Cdouble, maxfun::Cint, npt::Cint,
                              iprint::Cint, ctol::Cdouble,
                              callback::PrimaCallback, info::Ptr{Cint})::Status
end
