@enum prima_algorithm_t::Cint begin
    UOBYQA = 0
    NEWUOA = 1
    BOBYQA = 2
    LINCOA = 3
    COBYLA = 4
end

@enum prima_message_t::Cint begin
    MSG_NONE = 0
    MSG_EXIT = 1
    MSG_RHO = 2
    MSG_FEVL = 3
end

@enum prima_rc_t::Cint begin
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
    @ccall libprimac.prima_get_rc_string(rc::prima_rc_t)::Cstring
end

# typedef void ( * prima_obj_t ) ( const double x [ ] , double * const f , const void * data )
const prima_obj_t = Ptr{Cvoid}

# typedef void ( * prima_objcon_t ) ( const double x [ ] , double * const f , double constr [ ] , const void * data )
const prima_objcon_t = Ptr{Cvoid}

# typedef void ( * prima_callback_t ) ( const int n , const double x [ ] , const double f , const int nf , const int tr , const double cstrv , const int m_nlcon , const double nlconstr [ ] , bool * const terminate )
const prima_callback_t = Ptr{Cvoid}

struct prima_problem_t
    n::Cint
    calfun::prima_obj_t
    calcfc::prima_objcon_t
    x0::Ptr{Cdouble}
    xl::Ptr{Cdouble}
    xu::Ptr{Cdouble}
    m_ineq::Cint
    Aineq::Ptr{Cdouble}
    bineq::Ptr{Cdouble}
    m_eq::Cint
    Aeq::Ptr{Cdouble}
    beq::Ptr{Cdouble}
    m_nlcon::Cint
    f0::Cdouble
    nlconstr0::Ptr{Cdouble}
end

function prima_init_problem(problem, n)
    @ccall libprimac.prima_init_problem(problem::Ptr{prima_problem_t}, n::Cint)::Status
end

struct prima_options_t
    rhobeg::Cdouble
    rhoend::Cdouble
    maxfun::Cint
    iprint::Cint
    ftarget::Cdouble
    npt::Cint
    ctol::Cdouble
    data::Ptr{Cvoid}
    callback::prima_callback_t
end

function prima_init_options(options)
    @ccall libprimac.prima_init_options(options::Ptr{prima_options_t})::Status
end

struct prima_result_t
    x::Ptr{Cdouble}
    f::Cdouble
    cstrv::Cdouble
    nlconstr::Ptr{Cdouble}
    nf::Cint
    status::Cint
    message::Cstring
end

function prima_free_result(result)
    @ccall libprimac.prima_free_result(result::Ptr{prima_result_t})::Status
end

function prima_minimize(algorithm, problem, options, result)
    @ccall libprimac.prima_minimize(algorithm::prima_algorithm_t, problem::prima_problem_t,
                                    options::prima_options_t,
                                    result::Ptr{prima_result_t})::Status
end
