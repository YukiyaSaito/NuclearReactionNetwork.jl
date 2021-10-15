module NuclearReactionNetwork

# Write your package code here.
include("./LinearInterpolations.jl")
include("./LinearSolvers.jl")
include("./Astro.jl")
include("./Network.jl")
include("./ReactionTypes.jl")
include("./NetworkDatas.jl")
include("./InOut.jl")
include("./NetworkSolve.jl")

include("./main.jl")

import .LinearInterpolations
import .LinearSolvers
import .ReactionTypes
import .Network
import .NetworkDatas
import .NetworkSolve
import .InOut
import .main

# test
import .read_test

# test
export main
export read_test

end
