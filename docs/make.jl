using Documenter
using StarStats

makedocs(
    sitename = "StarStats",
     authors="Pablo Marchant <pablo.marchant@ugent.be>, Alina Istrate, Polina Smirnova <polina.smirnova@ugent.be>",
         repo="https://github.com/orlox/StarStats.jl/blob/{commit}{path}#{line}",
    format = Documenter.HTML(),
    modules = [StarStats],
    pages=[
        "Home" => "index.md",
        "Examples" => [],
        "Reference" => [
            "EquivalentEvolutionaryPoint.md",
            "SimulationData.md",
            "StellarModelGrid.md",
            "SimplexInterpolation.md",
            "StellarModelSet.md",
        ]
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/orlox/StarStats.jl"
)
