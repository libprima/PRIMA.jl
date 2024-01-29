module PRIMACUTEstExt

if isdefined(Base, :get_extension)
    using PRIMA, CUTEst
else
    using ..PRIMA, ..CUTEst
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
