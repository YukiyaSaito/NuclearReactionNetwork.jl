using Documenter

include("../src/LinearInterpolations.jl")
include("../src/LinearSolvers.jl")
include("../src/Astro.jl")
include("../src/Network.jl")
include("../src/ReactionTypes.jl")
include("../src/NetworkDatas.jl")
include("../src/InOut.jl")
include("../src/NetworkSolve.jl")
include("../src/main.jl")

makedocs(sitename="NuclearReactionNetwork.jl")