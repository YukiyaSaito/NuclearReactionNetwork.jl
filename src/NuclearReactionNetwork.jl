module NuclearReactionNetwork

# Write your package code here.
include("./ReactionTypes.jl")
include("./Network.jl")

include("./NetworkSolve.jl")

# include("./InOut.jl")

include("./main.jl")

import .ReactionTypes
import .Network
import .NetworkSolve
# import .InOut
import .main

# test
import .read_probdecay!
import .read_test



export initialize_reactions
export read_boundary
export ProbDecay
export initialize_jacobian

export main

# test
export get_networksize
export read_probdecay!
export read_test

end
