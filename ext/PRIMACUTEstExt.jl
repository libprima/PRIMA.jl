module PRIMACUTEstExt

using PRIMA, CUTEst

for func in (:uobyqa, :newuoa, :bobyqa, :lincoa, :cobyla, :prima)
    @eval function PRIMA.$(Symbol(func,"_CUTEst"))(name::AbstractString; kwargs...)
        nlp = CUTEstModel{Float64}(name)
        try
            return $func(nlp; kwargs...)
        finally
            finalize(nlp)
        end
    end
end

end # module
