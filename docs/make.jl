using GeomOpt
using Documenter

DocMeta.setdocmeta!(GeomOpt, :DocTestSetup, :(using GeomOpt); recursive=true)

makedocs(;
    modules=[GeomOpt],
    authors="Christoph Ortner <christohortner@gmail.com> and contributors",
    sitename="GeomOpt.jl",
    format=Documenter.HTML(;
        canonical="https://ACEsuit.github.io/GeomOpt.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/ACEsuit/GeomOpt.jl",
    devbranch="main",
)
