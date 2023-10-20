# Script to parse PRIMA headers and generate Julia wrappers.
using PRIMA_jll
using Clang
using Clang.Generators
using JuliaFormatter

function main()

    cd(@__DIR__)
    include = joinpath(PRIMA_jll.artifact_dir, "include", "prima")
    headers = [joinpath(include, "prima.h")]

    options = load_options(joinpath(@__DIR__, "prima.toml"))
    options["general"]["output_file_path"] = joinpath("..", "src", "wrappers.jl")

    args = get_default_args()
    push!(args, "-I$include")

    ctx = create_context(headers, args, options)
    build!(ctx)

    path = options["general"]["output_file_path"]
    code = readlines(path)
    for repl in [
        # Simplify the name of non-exported but public symbols (they will be
        # prefixed by the module name).
        r"\bprima_message\b" => "Message",
        r"\bprima_rc\b" => "Status",
        r"\bPRIMA_" => "",
        # Force enums to have signed value of type Cint (see PRIMA library doc.
        # about negative `iprint` values).
        r"^ *@enum +(\w+) *:: *U?Int\d+ +begin *$" => s"@enum \1::Cint begin",
        # All algorithms return a Status.
        r"\) *:: *Cint *$" => s")::Status",
        # Remove some useless code.
        r"^\s*const\s+PRIMAC_API\s*=.*$" => "",
        ]
        for i in eachindex(code)
            code[i] = replace(code[i], repl)
        end
    end
    open(path, "w") do io
        foreach(line -> println(io, line), code)
    end
    format_file(path, YASStyle())
    return nothing
end

# If we want to use the file as a script with `julia wrapper.jl`
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
