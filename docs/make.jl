using Documenter
using AttenuatedTotalReflectance

makedocs(
    sitename="AttenuatedTotalReflectance.jl",
    modules=[AttenuatedTotalReflectance],
    pages=[
       
        "Home" => "index.md",
        "Manual" => [
            "Getting Started" => "guide.md",
            "Physics Background" => "physics.md"
        ],
        "API Reference" => "api.md"  # You can put all your @docs blocks here
    ]
)
    


deploydocs(
    repo="github.com/Xyarc/AttenuatedTotalReflectance.jl.git",
    devbranch="main"
)