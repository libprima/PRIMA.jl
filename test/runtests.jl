using PRIMA
using Test, TypeUtils

optimizer_name(::typeof(PRIMA.uobyqa)) = "UOBYQA"
optimizer_name(::typeof(PRIMA.newuoa)) = "NEWUOA"
optimizer_name(::typeof(PRIMA.bobyqa)) = "BOBYQA"
optimizer_name(::typeof(PRIMA.cobyla)) = "COBYLA"
optimizer_name(::typeof(PRIMA.lincoa)) = "LINCOA"
optimizer_name(algo::Symbol) = optimizer_name(optimizer(algo))

optimizer(algo::Symbol) =
    algo === :uobyqa ? PRIMA.uobyqa :
    algo === :newuoa ? PRIMA.newuoa :
    algo === :bobyqa ? PRIMA.bobyqa :
    algo === :cobyla ? PRIMA.cobyla :
    algo === :lincoa ? PRIMA.lincoa :
    error("unknown optimizer `:$algo`")

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

        # Array to store non-linear constraints.
        c = Array{Cdouble}(undef, 1)

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
            println("\nNEWUOA:")
            x, fx, nf, rc = @inferred PRIMA.newuoa(f, x0;
                                                   rhobeg = 1.0, rhoend = 1e-3,
                                                   ftarget = -Inf,
                                                   maxfun = 200n,
                                                   npt = 2n + 1,
                                                   iprint = PRIMA.MSG_EXIT)
            msg = PRIMA.reason(rc)
            println("x = $x, f(x) = $fx, rc = $rc, msg = '$msg', evals = $nf")
            @test x ≈ [3,2] atol=2e-2 rtol=0
            @test f(x) ≈ fx
            @test x0 == x0_sav
        end

        @testset "UOBYQA" begin
            println("\nUOBYQA:")
            x, fx, nf, rc = @inferred PRIMA.uobyqa(f, x0;
                                                   rhobeg = 1.0, rhoend = 1e-3,
                                                   ftarget = -Inf,
                                                   maxfun = 200n,
                                                   iprint = PRIMA.MSG_EXIT)
            msg = PRIMA.reason(rc)
            println("x = $x, f(x) = $fx, rc = $rc, msg = '$msg', evals = $nf")
            @test x ≈ [3,2] atol=2e-2 rtol=0
            @test f(x) ≈ fx
            @test x0 == x0_sav
        end

        @testset "BOBYQA" begin
            println("\nBOBYQA:")
            x, fx, nf, rc = @inferred PRIMA.bobyqa(f, x0;
                                                   xl, xu,
                                                   rhobeg = 1.0, rhoend = 1e-3,
                                                   ftarget = -Inf,
                                                   maxfun = 200n,
                                                   npt = 2n + 1,
                                                   iprint = PRIMA.MSG_EXIT)
            msg = PRIMA.reason(rc)
            println("x = $x, f(x) = $fx, rc = $rc, msg = '$msg', evals = $nf")
            @test x ≈ [3,2] atol=2e-2 rtol=0
            @test f(x) ≈ fx
            @test x0 == x0_sav
            @test check_bounds(xl, x, xu)
        end

        @testset "COBYLA" begin
            println("\nCOBYLA:")
            # First call with just the number of non-linear inequality constraints.
            x, fx, nf, rc, cstrv = @inferred PRIMA.cobyla(f, x0;
                                                          nonlinear_ineq = c_ineq,
                                                          linear_ineq = (A_ineq, b_ineq),
                                                          xl, xu,
                                                          rhobeg = 1.0, rhoend = 1e-3,
                                                          ftarget = -Inf,
                                                          maxfun = 200*n,
                                                          iprint = PRIMA.MSG_EXIT)
            msg = PRIMA.reason(rc)
            println("x = $x, f(x) = $fx, cstrv = $cstrv, rc = $rc, msg = '$msg', evals = $nf")
            @test x ≈ [3,2] atol=2e-2 rtol=0
            @test f(x) ≈ fx
            @test x0 == x0_sav
            @test check_bounds(xl, x, xu)
            # Second call with an array to store the non-linear inequality constraints.
            x, fx, nf, rc, cstrv = @inferred PRIMA.cobyla(f, x0;
                                                          nonlinear_ineq = (c, c_ineq),
                                                          linear_ineq = (A_ineq, b_ineq),
                                                          xl, xu,
                                                          rhobeg = 1.0, rhoend = 1e-3,
                                                          ftarget = -Inf,
                                                          maxfun = 200*n,
                                                          iprint = PRIMA.MSG_EXIT)
            msg = PRIMA.reason(rc)
            println("x = $x, f(x) = $fx, cstrv = $cstrv, c(x) = $c, rc = $rc, msg = '$msg', evals = $nf")
            @test x ≈ [3,2] atol=2e-2 rtol=0
            @test f(x) ≈ fx
            @test x0 == x0_sav
            @test check_bounds(xl, x, xu)
        end

        @testset "LINCOA" begin
            println("\nLINCOA:")
            x, fx, nf, rc, cstrv = @inferred PRIMA.lincoa(f, x0;
                                                          linear_ineq = (A_ineq, b_ineq),
                                                          xl, xu,
                                                          rhobeg = 1.0, rhoend = 1e-3,
                                                          ftarget = -Inf,
                                                          maxfun = 200*n,
                                                          npt = 2n + 1,
                                                          iprint = PRIMA.MSG_EXIT)
            msg = PRIMA.reason(rc)
            println("x = $x, f(x) = $fx, cstrv = $cstrv, rc = $rc, msg = '$msg', evals = $nf")
            @test x ≈ [3,2] atol=2e-2 rtol=0
            @test f(x) ≈ fx
            @test x0 == x0_sav
            @test check_bounds(xl, x, xu)
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
        rhoend = 1e-5
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
                x, fx, nf, rc = @inferred uobyqa(f, x0;
                                                 rhobeg, rhoend, ftarget, maxfun, iprint)
            elseif optim === :newuoa
                x, fx, nf, rc = @inferred newuoa(f, x0;
                                                 rhobeg, rhoend, ftarget, maxfun, iprint, npt)
            elseif optim === :bobyqa
                x, fx, nf, rc = @inferred bobyqa(f, x0;
                                                 rhobeg, rhoend, ftarget, maxfun, iprint, npt)
            elseif optim === :cobyla
                x, fx, nf, rc = @inferred cobyla(f, x0;
                                                 rhobeg, rhoend, ftarget, maxfun, iprint)
            elseif optim === :lincoa
                x, fx, nf, rc = @inferred lincoa(f, x0;
                                                 rhobeg, rhoend, ftarget, maxfun, iprint, npt)
            else
                continue
            end
            msg = PRIMA.reason(rc)
            println("x = $x, f(x) = $fx, rc = $rc, msg = '$msg', evals = $nf")
            @test x ≈ [1,1] rtol=0 atol=(optim == :cobyla ? 3e-2 : 2e-2)
            @test f(x) ≈ fx

            # Bound constrained optimization.
            if optim ∈ (:bobyqa, :cobyla, :lincoa)
                println("\nBound constrained minimization of Rosenbrock function by $(optimizer_name(optim)):")
                x0 = [0.5, 2.5]
                xl = [-Inf, 1.2]
                xu = [+Inf, +Inf]
                if optim === :bobyqa
                    x, fx, nf, rc = @inferred bobyqa(f, x0;
                                                     xl, xu,
                                                     rhobeg, rhoend, ftarget, maxfun, iprint, npt)
                elseif optim === :cobyla
                    x, fx, nf, rc = @inferred cobyla(f, x0;
                                                     xl, xu,
                                                     rhobeg, rhoend, ftarget, maxfun, iprint)
                elseif optim === :lincoa
                    x, fx, nf, rc = @inferred lincoa(f, x0;
                                                     xl, xu,
                                                     rhobeg, rhoend, ftarget, maxfun, iprint, npt)
                else
                    continue
                end
                msg = PRIMA.reason(rc)
                println("x = $x, f(x) = $fx, rc = $rc, msg = '$msg', evals = $nf")
                @test x ≈ [1.095247,1.2] rtol=0 atol=2e-2
                @test f(x) ≈ fx
            end

            # Linearly constrained optimization.
            if optim ∈ (:cobyla, :lincoa)
                println("\nConstrained minimization of Rosenbrock function by $(optimizer_name(optim)):")
                x0 = [-1, 2] # starting point
                linear_ineq = (A_ineq, b_ineq)
                if optim === :cobyla
                    x, fx, nf, rc = @inferred cobyla(f, x0;
                                                     linear_ineq,
                                                     rhobeg, rhoend, ftarget, maxfun, iprint)
                elseif optim === :lincoa
                    x, fx, nf, rc = @inferred lincoa(f, x0;
                                                     linear_ineq,
                                                     rhobeg, rhoend, ftarget, maxfun, iprint, npt)
                else
                    continue
                end
                msg = PRIMA.reason(rc)
                println("x = $x, f(x) = $fx, rc = $rc, msg = '$msg', evals = $nf")
                @test x ≈ [1.0,1.0] rtol=0 atol=(optim == :cobyla ? 3e-2 : 2e-2)
                @test f(x) ≈ fx

                println("\nIdem but with one linear inequality constraint replaced by a bound constraint:")
                x0 = [1, 2] # starting point
                linear_ineq = (A_ineq[1:2,:], b_ineq[1:2])
                xl = [-1,-Inf]
                if optim === :cobyla
                    x, fx, nf, rc = @inferred cobyla(f, x0;
                                                     xl, linear_ineq,
                                                     rhobeg, rhoend, ftarget, maxfun, iprint)
                elseif optim === :lincoa
                    x, fx, nf, rc = @inferred lincoa(f, x0;
                                                     xl, linear_ineq,
                                                     rhobeg, rhoend, ftarget, maxfun, iprint, npt)
                else
                    continue
                end
                msg = PRIMA.reason(rc)
                println("x = $x, f(x) = $fx, rc = $rc, msg = '$msg', evals = $nf")
                @test x ≈ [1.0,1.0] rtol=0 atol=(optim == :cobyla ? 3e-2 : 2e-2)
                @test f(x) ≈ fx

                # The solution is on the line: 4x + 3y = 12 (the boundary of the first constraint).
                println("\nIdem but one linear constraint is active at the solution:")
                x0 = [1, 2] # starting point
                linear_ineq = (A_ineq, b_ineq)
                xl = [-1,-Inf]
                if optim === :cobyla
                    x, fx, nf, rc = @inferred cobyla(f2, x0;
                                                     linear_ineq,
                                                     rhobeg, rhoend, ftarget, maxfun, iprint)
                elseif optim === :lincoa
                    x, fx, nf, rc = @inferred lincoa(f2, x0;
                                                     linear_ineq,
                                                     rhobeg, rhoend, ftarget, maxfun, iprint, npt)
                else
                    continue
                end
                msg = PRIMA.reason(rc)
                println("x = $x, f(x) = $fx, rc = $rc, msg = '$msg', evals = $nf")
                @test x ≈ [1.441832,2.077557] rtol=0 atol=(optim == :cobyla ? 3e-2 : 2e-2)
                @test f2(x) ≈ fx
            end
        end
    end
end

nothing
