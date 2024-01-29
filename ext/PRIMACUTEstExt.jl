module PRIMACUTEstExt

if isdefined(Base, :get_extension)
    using NLPModels
    using CUTEst
    using PRIMA
else
    using ..NLPModels
    using ..CUTEst
    using ..PRIMA
end

for func in (:uobyqa, :newuoa, :bobyqa, :lincoa, :cobyla, :prima)
    @eval function PRIMA.$(Symbol(func,"_CUTEst"))(name::AbstractString; kwds...)
        nlp = CUTEstModel(name)
        try
            return $func(nlp; kwds...)
        finally
            finalize(nlp)
        end
    end
end

end # module
