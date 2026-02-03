using Documenter
using AttenuatedTotalReflectance

makedocs(
    sitename="AttenuatedTotalReflectance.jl",
    modules=[AttenuatedTotalReflectance],
    pages=[
        "Home" => "index.md",
        # Add more pages here as you grow
    ]
)

deploydocs(
    repo="github.com/Xyarc/AttenuatedTotalReflectance.jl.git",
    devbranch="main"
)