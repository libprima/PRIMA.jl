using PRIMA
using Test, TypeUtils
Sys.WORD_SIZE > 32 && using CUTEst

optimizer_name(::typeof(PRIMA.uobyqa)) = "UOBYQA"
optimizer_name(::typeof(PRIMA.newuoa)) = "NEWUOA"
optimizer_name(::typeof(PRIMA.bobyqa)) = "BOBYQA"
optimizer_name(::typeof(PRIMA.cobyla)) = "COBYLA"
optimizer_name(::typeof(PRIMA.lincoa)) = "LINCOA"
optimizer_name(::typeof(PRIMA.prima))  = "PRIMA"
optimizer_name(algo::Symbol) = optimizer_name(optimizer(algo))

optimizer(algo::Symbol) =
    algo === :uobyqa ? PRIMA.uobyqa :
    algo === :newuoa ? PRIMA.newuoa :
    algo === :bobyqa ? PRIMA.bobyqa :
    algo === :cobyla ? PRIMA.cobyla :
    algo === :lincoa ? PRIMA.lincoa :
    algo === :prima  ? PRIMA.prima  :
    error("unknown optimizer `:$algo`")

print_1(x::AbstractVector, info::PRIMA.Info) = print_1(stdout, x, info)
function print_1(io::IO, x::AbstractVector, info::PRIMA.Info)
    msg = PRIMA.reason(info)
    println(io, "x = $x, f(x) = $(info.fx), status = $(info.status), msg = '$msg', evals = $(info.nf)")
end

print_2(x::AbstractVector, info::PRIMA.Info) = print_2(stdout, x, info)
function print_2(io::IO, x::AbstractVector, info::PRIMA.Info)
    msg = PRIMA.reason(info)
    println(io, "x = $x, f(x) = $(info.fx), cstrv = $(info.cstrv), status = $(info.status), msg = '$msg', evals = $(info.nf)")
end

print_3(x::AbstractVector, info::PRIMA.Info) = print_3(stdout, x, info)
function print_3(io::IO, x::AbstractVector, info::PRIMA.Info)
    msg = PRIMA.reason(info)
    println("x = $x, f(x) = $(info.fx), cstrv = $(info.cstrv), c(x) = $(info.nl_ineq), status = $(info.status), msg = '$msg', evals = $(info.nf)")
end

@testset "PRIMA.jl" begin
    @testset "Utils " begin
        for status in instances(PRIMA.Status)
            @test PRIMA.reason(status) isa String
        end
        # Check NULL-array API.
        let NullArray = PRIMA.NullArray
            A = NullArray{Int16,3}()
            @test eltype(A) === Int16
            @test size(A) == (0,0,0)
            @test axes(A) == (1:0,1:0,1:0)
            @test length(A) == 0
            @test pointer(A) === Ptr{Int16}(0)
        end
    end

    @testset "Simple objective function" begin
        # User defined objective function:
        function f(x::AbstractVector{T}) where {T<:AbstractFloat}
            @assert length(x) == 2
            x1, x2 = x
            return as(T, 5*(x1 - 3)*(x1 - 3) + 7*(x2 - 2)*(x2 - 2) + as(T, 1//10)*(x1 + x2) - 10)
        end

        # Non-linear inequality constraints.
        function c_ineq(x::AbstractVector{T}) where {T<:AbstractFloat}
            @assert length(x) == 2
            x1, x2 = x
            return x1*x1 + x2*x2 - 13 # ‖x‖² ≤ 13
        end

        # Initial solution.
        x0 = [0.0, 0.0]
        n = length(x0)
        x0_sav = copy(x0) # for testing that x0 does not get overwritten
        scl = 2.0 # scaling factor, a power of two should not change results
        scale = [scl, scl]

        # Inequality constraints: x₁ ≤ 4, x₂ ≤ 3, x₁ + x₂ ≤ 10
        A_ineq = [1.0 0.0;
                  0.0 1.0;
                  1.0 1.0]
        b_ineq = [4.0, 3.0, 10.0]

        # Bounds.
        xl = [-6.0, -6.0]
        xu = [ 6.0,  6.0]

        check_bounds(xl, x, xu) = all(xl .≤ x .≤ xu)

        @testset "NEWUOA" begin
            println("\nNEWUOA:")
            kwds = (rhobeg = 1.0, rhoend = 1e-6, ftarget = -Inf,
                    maxfun = 500n, npt = 2n + 1, iprint = PRIMA.MSG_EXIT)
            x, info = @inferred PRIMA.newuoa(f, x0; kwds...)
            print_1(x, info)
            @test issuccess(info)
            @test x ≈ [3,2] atol=2e-2 rtol=0
            @test f(x) ≈ info.fx
            @test x0 == x0_sav
            # Solve problem with general driver.
            x1, info1 = @inferred PRIMA.prima(f, x0; kwds...)
            @test x1 == x
            @test info1 == info
            # Solve problem with scaling factors.
            kwds = (scale, rhobeg = 1.0/scl, rhoend = 1e-6/scl, ftarget = -Inf,
                    maxfun = 500n, npt = 2n + 1, iprint = PRIMA.MSG_EXIT)
            x1, info1 = @inferred PRIMA.newuoa(f, x0; kwds...)
            @test x1 ≈ x
        end

        @testset "UOBYQA" begin
            println("\nUOBYQA:")
            kwds = (rhobeg = 1.0, rhoend = 1e-6, ftarget = -Inf,
                    maxfun = 500n, iprint = PRIMA.MSG_EXIT)
            x, info = @inferred PRIMA.uobyqa(f, x0; kwds...)
            print_1(x, info)
            @test issuccess(info)
            @test x ≈ [3,2] atol=2e-2 rtol=0
            @test f(x) ≈ info.fx
            @test x0 == x0_sav
            # Solve problem with scaling factors.
            kwds = (scale, rhobeg = 1.0/scl, rhoend = 1e-6/scl, ftarget = -Inf,
                    maxfun = 500n, iprint = PRIMA.MSG_EXIT)
            x1, info1 = @inferred PRIMA.uobyqa(f, x0; kwds...)
            @test x1 ≈ x
        end

        @testset "BOBYQA" begin
            println("\nBOBYQA:")
            kwds = (xl = xl, xu = xu,
                    rhobeg = 1.0, rhoend = 1e-6, ftarget = -Inf,
                    maxfun = 500n, npt = 2n + 1, iprint = PRIMA.MSG_EXIT)
            x, info = @inferred PRIMA.bobyqa(f, x0; kwds...)
            print_1(x, info)
            @test issuccess(info)
            @test x ≈ [3,2] atol=2e-2 rtol=0
            @test f(x) ≈ info.fx
            @test x0 == x0_sav
            @test check_bounds(xl, x, xu)
            # Solve problem with general driver.
            x1, info1 = @inferred PRIMA.prima(f, x0; kwds...)
            @test x1 == x
            @test info1 == info
            # Solve problem with scaling factors.
            kwds = (scale, rhobeg = 1.0/scl, rhoend = 1e-6/scl, ftarget = -Inf,
                    maxfun = 500n, npt = 2n + 1, iprint = PRIMA.MSG_EXIT)
            x1, info1 = @inferred PRIMA.bobyqa(f, x0; kwds...)
            @test x1 ≈ x
        end

        @testset "COBYLA" begin
            println("\nCOBYLA:")
            kwds = (xl = xl, xu = xu, linear_ineq = (A_ineq, b_ineq),
                    rhobeg = 1.0, rhoend = 1e-6, ftarget = -Inf,
                    maxfun = 500n, iprint = PRIMA.MSG_EXIT)
            # First call with just the number of non-linear inequality constraints.
            x, info = @inferred PRIMA.cobyla(f, x0; kwds...,
                                             nonlinear_ineq = c_ineq)
            print_3(x, info)
            @test issuccess(info)
            @test x ≈ [3,2] atol=2e-2 rtol=0
            @test f(x) ≈ info.fx
            @test x0 == x0_sav
            @test check_bounds(xl, x, xu)
            # Call with given number of non-linear inequality constraints.
            x1, info1 = @inferred PRIMA.cobyla(f, x0; kwds...,
                                               nonlinear_ineq = (length(c_ineq(x0)), c_ineq))
            @test x1 == x
            @test info1 == info
            # Solve problem with general driver.
            x1, info1 = @inferred PRIMA.prima(f, x0; kwds...,
                                              nonlinear_ineq = c_ineq)
            @test x1 == x
            @test info1 == info
            # Solve problem with scaling factors.
            kwds = (xl = xl, xu = xu, linear_ineq = (A_ineq, b_ineq),
                    scale, rhobeg = 1.0/scl, rhoend = 1e-6/scl, ftarget = -Inf,
                    maxfun = 500n, iprint = PRIMA.MSG_EXIT)
            x1, info1 = @inferred PRIMA.cobyla(f, x0; kwds...,
                                               nonlinear_ineq = c_ineq)
            @test x1 ≈ x
        end

        @testset "LINCOA" begin
            println("\nLINCOA:")
            kwds = (xl = xl, xu = xu, linear_ineq = (A_ineq, b_ineq),
                    rhobeg = 1.0, rhoend = 1e-6, ftarget = -Inf,
                    maxfun = 500n, npt = 2n + 1, iprint = PRIMA.MSG_EXIT)
            x, info = @inferred PRIMA.lincoa(f, x0; kwds...)
            print_2(x, info)
            @test issuccess(info)
            @test x ≈ [3,2] atol=2e-2 rtol=0
            @test f(x) ≈ info.fx
            @test x0 == x0_sav
            @test check_bounds(xl, x, xu)
            # Solve problem with general driver.
            x1, info1 = @inferred PRIMA.prima(f, x0; kwds...)
            @test x1 == x
            @test info1 == info
            # Solve problem with scaling factors.
            kwds = (xl = xl, xu = xu, linear_ineq = (A_ineq, b_ineq),
                    scale, rhobeg = 1.0/scl, rhoend = 1e-6/scl, ftarget = -Inf,
                    maxfun = 500n, npt = 2n + 1, iprint = PRIMA.MSG_EXIT)
            x1, info1 = @inferred PRIMA.lincoa(f, x0; kwds...)
            @test x1 ≈ x
        end
    end

    @testset "Rosenbrock" begin
        """
            Rosenbrock(x, y; a=1, b=100, hard=false) -> (a - x)² + b⋅(y - x²)²

        yields the value of the Rosenbrock function at coordinates `(x,y)`. The
        Rosenbrock function is a non-convex function with a global minimum at `(a,a²)`.

        If `x` is a vector, then:

            Rosenbrock(x; a=1, b=100, hard=false)

        yields the generalized Rosenbrock function.

        """
        function Rosenbrock(x1::Real, x2::Real; a::Real=1, b::Real=100)
            T = float(promote_type(typeof(x1), typeof(x2)))
            x1 = T(x1)
            t1 = T(a) - x1
            t2 = T(x2) - x1*x1
            return t1*t1 + T(b)*t2*t2
        end
        function Rosenbrock(x::AbstractVector{T};
                            a::Real=1, b::Real=100, hard::Bool=false) where {T<:AbstractFloat}
            hard || @assert iseven(length(x))
            stp = (hard ? 1 : 2)
            a, b = T(a), T(b)
            fx = zero(T)
            @inbounds for i ∈ firstindex(x):stp:lastindex(x)-1
                fx += Rosenbrock(x[i], x[i+1]; a, b)
            end
            return fx
        end

        # Initial solution and settings. NOTE: The problem is rather difficult
        # for COBYLA, hence `maxfun` is very large and `rhoend` very small.
        f(x) = Rosenbrock(x)
        f(x, cx) = f(x)
        f2(x) = Rosenbrock(x; a=2) # unconstrained global min. at (2,4)
        f2(x, cx) = f2(x)
        n = 2
        rhobeg = 1.0
        rhoend = 1e-6
        ftarget = -Inf
        maxfun = 3000n
        npt = 2n + 1
        iprint = PRIMA.MSG_EXIT

        # Linear inequality constraints
        # defining a closed convex region delimited by:
        #
        #     2x - 3y ≤  6
        #     4x + 3y ≤ 12
        #     -x ≤ 1
        A_ineq = [ 2 -3; 4  3; -1  0];
        b_ineq = [6, 12, 1]

        @testset "$(optimizer_name(optim))" for optim in (:uobyqa, :newuoa, :bobyqa, :cobyla, :lincoa)

            println("\nUnconstrained minimization of Rosenbrock function by $(optimizer_name(optim)):")
            x0 = [-1, 2]
            if optim === :uobyqa
                x, info = @inferred uobyqa(f, x0; rhobeg, rhoend, ftarget, maxfun, iprint)
            elseif optim === :newuoa
                x, info = @inferred newuoa(f, x0; rhobeg, rhoend, ftarget, maxfun, iprint, npt)
            elseif optim === :bobyqa
                x, info = @inferred bobyqa(f, x0; rhobeg, rhoend, ftarget, maxfun, iprint, npt)
            elseif optim === :cobyla
                x, info = @inferred cobyla(f, x0; rhobeg, rhoend, ftarget, maxfun, iprint)
            elseif optim === :lincoa
                x, info = @inferred lincoa(f, x0; rhobeg, rhoend, ftarget, maxfun, iprint, npt)
            else
                continue
            end
            print_1(x, info)
            @test issuccess(info)
            @test x ≈ [1,1] rtol=0 atol=(optim == :cobyla ? 3e-2 : 2e-2)
            @test f(x) ≈ info.fx

            # Bound constrained optimization.
            if optim ∈ (:bobyqa, :cobyla, :lincoa)
                println("\nBound constrained minimization of Rosenbrock function by $(optimizer_name(optim)):")
                x0 = [0.5, 2.5]
                xl = [-Inf, 1.2]
                xu = [+Inf, +Inf]
                if optim === :bobyqa
                    x, info = @inferred bobyqa(f, x0; xl, xu,
                                               rhobeg, rhoend, ftarget, maxfun, iprint, npt)
                elseif optim === :cobyla
                    x, info = @inferred cobyla(f, x0; xl, xu,
                                               rhobeg, rhoend, ftarget, maxfun, iprint)
                elseif optim === :lincoa
                    x, info = @inferred lincoa(f, x0; xl, xu,
                                               rhobeg, rhoend, ftarget, maxfun, iprint, npt)
                else
                    continue
                end
                print_1(x, info)
                @test issuccess(info)
                @test x ≈ [1.095247,1.2] rtol=0 atol=2e-2
                @test f(x) ≈ info.fx
            end

            # Linearly constrained optimization.
            if optim ∈ (:cobyla, :lincoa)
                println("\nConstrained minimization of Rosenbrock function by $(optimizer_name(optim)):")
                x0 = [-1, 2] # starting point
                linear_ineq = (A_ineq, b_ineq)
                if optim === :cobyla
                    x, info = @inferred cobyla(f, x0; linear_ineq,
                                               rhobeg, rhoend, ftarget, maxfun, iprint)
                elseif optim === :lincoa
                    x, info = @inferred lincoa(f, x0; linear_ineq,
                                               rhobeg, rhoend, ftarget, maxfun, iprint, npt)
                else
                    continue
                end
                print_1(x, info)
                @test issuccess(info)
                @test x ≈ [1.0,1.0] rtol=0 atol=(optim == :cobyla ? 3e-2 : 2e-2)
                @test f(x) ≈ info.fx

                println("\nIdem but with one linear inequality constraint replaced by a bound constraint:")
                x0 = [1, 2] # starting point
                linear_ineq = (A_ineq[1:2,:], b_ineq[1:2])
                xl = [-1,-Inf]
                if optim === :cobyla
                    x, info = @inferred cobyla(f, x0; xl, linear_ineq,
                                               rhobeg, rhoend, ftarget, maxfun, iprint)
                elseif optim === :lincoa
                    x, info = @inferred lincoa(f, x0; xl, linear_ineq,
                                               rhobeg, rhoend, ftarget, maxfun, iprint, npt)
                else
                    continue
                end
                print_1(x, info)
                @test issuccess(info)
                @test x ≈ [1.0,1.0] rtol=0 atol=(optim == :cobyla ? 3e-2 : 2e-2)
                @test f(x) ≈ info.fx

                # The solution is on the line: 4x + 3y = 12 (the boundary of the first constraint).
                println("\nIdem but one linear constraint is active at the solution:")
                x0 = [1, 2] # starting point
                linear_ineq = (A_ineq, b_ineq)
                xl = [-1,-Inf]
                if optim === :cobyla
                    x, info = @inferred cobyla(f2, x0; linear_ineq,
                                               rhobeg, rhoend, ftarget, maxfun, iprint)
                elseif optim === :lincoa
                    x, info = @inferred lincoa(f2, x0; linear_ineq,
                                               rhobeg, rhoend, ftarget, maxfun, iprint, npt)
                else
                    continue
                end
                print_1(x, info)
                @test issuccess(info)
                @test x ≈ [1.441832,2.077557] rtol=0 atol=(optim == :cobyla ? 3e-2 : 2e-2)
                @test f2(x) ≈ info.fx
            end
        end
    end

    @testset "User examples" begin
        # See https://github.com/libprima/PRIMA.jl/issues/19
        let cost_func(x) = sum(abs2, x), x0 = randn(4),
            opts = (rhobeg=0.1, rhoend=1e-8)
            @testset "$algo" for algo in (:uobyqa, :newuoa, :bobyqa,
                                          :cobyla, :lincoa, :prima)
                optim = optimizer(algo)
                x1, res = @inferred optim(cost_func, x0; opts...)
                @test issuccess(res)
                @test maximum(abs.(x1)) ≤ 1e-8
            end
        end
    end

    if Sys.WORD_SIZE > 32
        @testset "Unconstrained CUTEst problem $name" for name in ("TOINTQOR", "OSBORNEB", "LANCZOS1LS",)
            x1, res1 = @inferred PRIMA.prima_CUTEst(name; maxfun=5000)
            @test issuccess(res1)
            x2, res2 = @inferred PRIMA.newuoa_CUTEst(name; maxfun=5000)
            @test issuccess(res2)
            @test x1 ≈ x2
        end
    end

end

nothing
