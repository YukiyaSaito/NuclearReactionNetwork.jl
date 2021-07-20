using NuclearReactionNetwork
using Test
using Profile

@testset "NuclearReactionNetwork.jl" begin
    # Write your tests here.
    # read_test("/Users/yukiya/Documents/Physics/PhD/ReactionNetwork.jl/moller2003.dat")
    @time networksize, timestep = read_test("/Users/yukiya/Documents/Physics/PhD/prism/input/extent/default.dat")
    # Profile.print()
    # @test networksize == 13396
    # @test typeof(timestep)==Float64
end

