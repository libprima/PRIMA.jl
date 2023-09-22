@enum Message::Int32 begin
    MSG_NONE = 0
    MSG_EXIT = 1
    MSG_RHO = 2
    MSG_FEVL = 3
end

@enum Status::Int32 begin
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
    INVALID_INPUT = 100
    ASSERTION_FAILS = 101
    VALIDATION_FAILS = 102
    MEMORY_ALLOCATION_FAILS = 103
end

function prima_get_rc_string(rc)
    @ccall libprimac.prima_get_rc_string(rc::Cint)::Cstring
end

# typedef void ( * prima_obj ) ( const double x [ ] , double * f )
const prima_obj = Ptr{Cvoid}

# typedef void ( * prima_objcon ) ( const double x [ ] , double * f , double constr [ ] )
const prima_objcon = Ptr{Cvoid}

function prima_bobyqa(calfun, n, x, f, xl, xu, nf, rhobeg, rhoend, ftarget, maxfun, npt,
                      iprint)
    @ccall libprimac.prima_bobyqa(calfun::prima_obj, n::Cint, x::Ptr{Cdouble},
                                  f::Ptr{Cdouble}, xl::Ptr{Cdouble}, xu::Ptr{Cdouble},
                                  nf::Ptr{Cint}, rhobeg::Cdouble, rhoend::Cdouble,
                                  ftarget::Cdouble, maxfun::Cint, npt::Cint,
                                  iprint::Cint)::Cint
end

function prima_newuoa(calfun, n, x, f, nf, rhobeg, rhoend, ftarget, maxfun, npt, iprint)
    @ccall libprimac.prima_newuoa(calfun::prima_obj, n::Cint, x::Ptr{Cdouble},
                                  f::Ptr{Cdouble}, nf::Ptr{Cint}, rhobeg::Cdouble,
                                  rhoend::Cdouble, ftarget::Cdouble, maxfun::Cint,
                                  npt::Cint, iprint::Cint)::Cint
end

function prima_uobyqa(calfun, n, x, f, nf, rhobeg, rhoend, ftarget, maxfun, iprint)
    @ccall libprimac.prima_uobyqa(calfun::prima_obj, n::Cint, x::Ptr{Cdouble},
                                  f::Ptr{Cdouble}, nf::Ptr{Cint}, rhobeg::Cdouble,
                                  rhoend::Cdouble, ftarget::Cdouble, maxfun::Cint,
                                  iprint::Cint)::Cint
end

function prima_cobyla(m_nlcon, calcfc, n, x, f, cstrv, nlconstr, m_ineq, Aineq, bineq, m_eq,
                      Aeq, beq, xl, xu, nf, rhobeg, rhoend, ftarget, maxfun, iprint)
    @ccall libprimac.prima_cobyla(m_nlcon::Cint, calcfc::prima_objcon, n::Cint,
                                  x::Ptr{Cdouble}, f::Ptr{Cdouble}, cstrv::Ptr{Cdouble},
                                  nlconstr::Ptr{Cdouble}, m_ineq::Cint, Aineq::Ptr{Cdouble},
                                  bineq::Ptr{Cdouble}, m_eq::Cint, Aeq::Ptr{Cdouble},
                                  beq::Ptr{Cdouble}, xl::Ptr{Cdouble}, xu::Ptr{Cdouble},
                                  nf::Ptr{Cint}, rhobeg::Cdouble, rhoend::Cdouble,
                                  ftarget::Cdouble, maxfun::Cint, iprint::Cint)::Cint
end

function prima_lincoa(calfun, n, x, f, cstrv, m_ineq, Aineq, bineq, m_eq, Aeq, beq, xl, xu,
                      nf, rhobeg, rhoend, ftarget, maxfun, npt, iprint)
    @ccall libprimac.prima_lincoa(calfun::prima_obj, n::Cint, x::Ptr{Cdouble},
                                  f::Ptr{Cdouble}, cstrv::Ptr{Cdouble}, m_ineq::Cint,
                                  Aineq::Ptr{Cdouble}, bineq::Ptr{Cdouble}, m_eq::Cint,
                                  Aeq::Ptr{Cdouble}, beq::Ptr{Cdouble}, xl::Ptr{Cdouble},
                                  xu::Ptr{Cdouble}, nf::Ptr{Cint}, rhobeg::Cdouble,
                                  rhoend::Cdouble, ftarget::Cdouble, maxfun::Cint,
                                  npt::Cint, iprint::Cint)::Cint
end
