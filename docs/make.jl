push!(LOAD_PATH,"../src/")
using Documenter, DynACof

makedocs(sitename="DynACof")
deploydocs(;
    repo="github.com/VEZY/DynACof.jl.git",
)
