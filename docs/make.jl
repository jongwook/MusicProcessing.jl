
push!(LOAD_PATH,"../src/")
using MusicProcessing
using Documenter

DocMeta.setdocmeta!(MusicProcessing, :DocTestSetup, :(using MusicProcessing); recursive=true)

makedocs(;
    modules=[MusicProcessing],
    authors="JuliaMusic",
    repo="github.com/JuliaMusic/MusicProcessing.jl/blob/{commit}{path}#{line}",
    sitename="MusicProcessing.jl",
    format=Documenter.HTML(;
        prettyurls=Base.get(ENV, "CI", "false") == "true",
        canonical="github.com/JuliaMusic/MusicProcessing.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaMusic/MusicProcessing.jl",
    devbranch="main",
    push_preview = true
) 