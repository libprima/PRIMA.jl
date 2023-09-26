using PRIMA
using Test, TypeUtils

@testset "PRIMA.jl" begin
    @testset "Utils " begin
        for status in instances(PRIMA.Status)
            @test PRIMA.reason(status) isa String
        end
    end

    # User defined objective function:
    function f(x::AbstractVector{T}) where {T<:AbstractFloat}
        @assert length(x) == 2
        x1, x2 = x
        return as(T, 5*(x1 - 3)*(x1 - 3) + 7*(x2 - 2)*(x2 - 2) + as(T, 1//10)*(x1 + x2) - 10)
    end

    # Objective function for COBYLA (same name but different signature to implement non-linear constraint).
    function f(x::AbstractVector{T}, c::AbstractVector{T}) where {T<:AbstractFloat}
        nlconstr!(x, c)
        return f(x)
    end

    # Non-linear constraint.
    function nlconstr!(x::AbstractVector{T}, c::AbstractVector{T}) where {T<:AbstractFloat}
        @assert length(x) == 2
        @assert length(c) == 1
        x1, x2 = x
        c[firstindex(c)] = x1*x1 + x2*x2 - 13 # ‖x‖² ≤ 13
        return nothing
    end

    # Array to store non-linear constraints.
    nlconstr = Array{Cdouble}(undef, 1)

    # Initial solution.
    x0 = [0.0, 0.0]
    n = length(x0)
    x0_sav = copy(x0) # for testing that x0 does not get overwritten

    # Inequality constraints: x₁ ≤ 4, x₂ ≤ 3, x₁ + x₂ ≤ 10
    A_ineq = [1.0 0.0;
              0.0 1.0;
              1.0 1.0]
    b_ineq = [4.0, 3.0, 10.0]

    # Bounds.
    xl = [-6.0, -6.0]
    xu = [ 6.0,  6.0]

    function check_bounds(xl, x, xu)
        flag = true
        for (lᵢ, xᵢ, uᵢ) in zip(xl, x, xu)
            flag &= lᵢ ≤ xᵢ ≤ uᵢ
        end
        return flag
    end

    @testset "NEWUOA" begin
        x, fx, nf, rc = @inferred PRIMA.newuoa(f, x0;
                                               rhobeg = 1.0, rhoend = 1e-3,
                                               ftarget = -Inf,
                                               maxfun = 200n,
                                               npt = 2n + 1,
                                               iprint = PRIMA.MSG_EXIT)
        msg = PRIMA.reason(rc)
        println("x = $x, f(x) = $fx, rc = $rc, msg = '$msg', evals = $nf")
        @test abs(x[1] - 3) ≤ 2e-2 && abs(x[2] - 2) ≤ 2e-2
        @test f(x) ≈ fx
        @test x0 == x0_sav
    end

    @testset "UOBYQA" begin
        x, fx, nf, rc = @inferred PRIMA.uobyqa(f, x0;
                                               rhobeg = 1.0, rhoend = 1e-3,
                                               ftarget = -Inf,
                                               maxfun = 200n,
                                               iprint = PRIMA.MSG_EXIT)
        msg = PRIMA.reason(rc)
        println("x = $x, f(x) = $fx, rc = $rc, msg = '$msg', evals = $nf")
        @test abs(x[1] - 3) ≤ 2e-2 && abs(x[2] - 2) ≤ 2e-2
        @test f(x) ≈ fx
        @test x0 == x0_sav
    end

    @testset "BOBYQA" begin
        x, fx, nf, rc = @inferred PRIMA.bobyqa(f, x0;
                                               xl, xu,
                                               rhobeg = 1.0, rhoend = 1e-3,
                                               ftarget = -Inf,
                                               maxfun = 200n,
                                               npt = 2n + 1,
                                               iprint = PRIMA.MSG_EXIT)
        msg = PRIMA.reason(rc)
        println("x = $x, f(x) = $fx, rc = $rc, msg = '$msg', evals = $nf")
        @test abs(x[1] - 3) ≤ 2e-2 && abs(x[2] - 2) ≤ 2e-2
        @test f(x) ≈ fx
        @test x0 == x0_sav
        @test check_bounds(xl, x, xu)
    end

    @testset "COBYLA" begin
        x, fx, nf, rc, cstrv = @inferred PRIMA.cobyla(f, x0;
                                                      nlconstr,
                                                      ineqconstr = (A_ineq, b_ineq),
                                                      xl, xu,
                                                      rhobeg = 1.0, rhoend = 1e-3,
                                                      ftarget = -Inf,
                                                      maxfun = 200*n,
                                                      iprint = PRIMA.MSG_EXIT)
        msg = PRIMA.reason(rc)
        println("x = $x, f(x) = $fx, cstrv = $cstrv, nlconstr = $nlconstr, rc = $rc, msg = '$msg', evals = $nf")
        @test abs(x[1] - 3) ≤ 2e-2 && abs(x[2] - 2) ≤ 2e-2
        @test f(x) ≈ fx
        @test x0 == x0_sav
        @test check_bounds(xl, x, xu)
    end


    @testset "LINCOA" begin
        x, fx, nf, rc, cstrv = @inferred PRIMA.lincoa(f, x0;
                                                      ineqconstr = (A_ineq, b_ineq),
                                                      xl, xu,
                                                      rhobeg = 1.0, rhoend = 1e-3,
                                                      ftarget = -Inf,
                                                      maxfun = 200*n,
                                                      npt = 2n + 1,
                                                      iprint = PRIMA.MSG_EXIT)
        msg = PRIMA.reason(rc)
        println("x = $x, f(x) = $fx, cstrv = $cstrv, rc = $rc, msg = '$msg', evals = $nf")
        @test abs(x[1] - 3) ≤ 2e-2 && abs(x[2] - 2) ≤ 2e-2
        @test f(x) ≈ fx
        @test x0 == x0_sav
        @test check_bounds(xl, x, xu)
    end

end

nothing
