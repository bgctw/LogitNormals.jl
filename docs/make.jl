using LogitNormals
using Documenter

DocMeta.setdocmeta!(LogitNormals, :DocTestSetup, :(using LogitNormals); recursive=true)

makedocs(;
    modules=[LogitNormals],
    authors="Thomas Wutzler <twutz@bgc-jena.mpg.de> and contributors",
    repo="https://github.com/bgctw/LogitNormals.jl/blob/{commit}{path}#{line}",
    sitename="LogitNormals.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://bgctw.github.io/LogitNormals.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/bgctw/LogitNormals.jl",
    devbranch="main",
)
