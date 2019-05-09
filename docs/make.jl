using Documenter, DynACof

makedocs(;
    modules=[DynACof],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/VEZY/DynACof.jl/blob/{commit}{path}#L{line}",
    sitename="DynACof.jl",
    authors="remi.vezy",
    assets=[],
)

deploydocs(;
    repo="github.com/VEZY/DynACof.jl",
)
